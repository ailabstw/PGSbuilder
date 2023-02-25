#!/usr/bin/python3

import sys, os, random, json, pyreadr, joblib, copy
import numpy as np
import pandas as pd
import sklearn.metrics as metrics
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.base import clone
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from collections import OrderedDict
from glob import glob


# build the beta dataframe from all PRS algorithms except GenEpi (self.ss_df)
class Weights():
    def __init__(self, bfile, ss_file, prs_dir):
        self.basename = bfile.split('/')[-1]
        self.prs_dir = prs_dir

        # bim_df
        self.bim_df = pd.read_csv('{}.bim'.format(bfile), sep='\s+', names=['CHR', 'ID', 'CM', 'POS', 'ALT', 'REF'])
        self.bim_df = self.bim_df[['CHR', 'POS', 'ID', 'REF', 'ALT']]

        # ss_df
        self.ss_df = pd.read_csv(ss_file, sep='\s+')
        self.ss_df = self.ss_df[['#CHROM','POS','ID','REF','ALT','A1','P','LOG10_P', 'BETA']]
        self.ss_df = self.bim_df.merge(self.ss_df[['ID', 'A1', 'P', 'LOG10_P', 'BETA']], on='ID', how='left')

        # fill NA
        self.ss_df.fillna(0, inplace=True)


    def __call__(self):
        print('\n\n###### Building Beta (Weight) Dataframe ######\n\n')
        ### C+T
        # threshold
        try:
            with open('{}/CandT/best_pvalue_range'.format(self.prs_dir), 'r') as f:
                threshold = f.readlines()[0].split()[0]
        except:
            threshold = 0
            pass

        # beta
        snp_file = '{}/CandT/{}.valid.snp'.format(self.prs_dir, self.basename)
        if os.path.isfile(snp_file):
            print('Loading C+T ...')
            with open(snp_file, 'r') as f:
                snp_list = ''.join(f.readlines()).split()
            self._get_ct_beta('CandT', snp_list, float(threshold))


        ### PRSice2
        # threshold
        try:
            with open('{}/PRSice2/best_pvalue_range'.format(self.prs_dir), 'r') as f:
                threshold = f.readlines()[0].split()[0]
        except:
            threshold = 0
            pass

        # beta
        snp_file = '{}/PRSice2/{}.valid.snp'.format(self.prs_dir, self.basename)
        if os.path.isfile(snp_file):
            print('Loading PRSice2 ...')
            with open(snp_file, 'r') as f:
                snp_list = ''.join(f.readlines()).split()
            self._get_ct_beta('PRSice2', snp_list, float(threshold))
        

        ### Lassosum
        beta_file = '{}/Lassosum/{}.beta'.format(self.prs_dir, self.basename)
        if os.path.isfile(beta_file):
            print('Loading Lassosum ...')
            with open(beta_file, 'r') as f:
                beta_list = ''.join(f.readlines()).split()
                beta_list = list(map(float, beta_list))
            self.ss_df['Lassosum'] = beta_list
            self.ss_df.loc[self.ss_df['A1']!=self.ss_df['ALT'], 'Lassosum'] *= -1
        

        ### LDpred2
        beta_prefix = '{}/LDpred2/{}'.format(self.prs_dir, self.basename)
        if os.path.isfile('{}.beta.rds'.format(beta_prefix)):
            print('Loading LDpred2 ...')
            self._get_ldpred2_beta(beta_prefix)
            self.ss_df.loc[self.ss_df['A1']!=self.ss_df['ALT'], 'LDpred2'] *= -1
        

        ### PRScs
        beta_file = '{}/PRScs/effect_size.txt'.format(self.prs_dir)
        if os.path.isfile(beta_file):
            print('Loading PRScs ...')
            self._get_prscs_beta(beta_file)
            self.ss_df.loc[self.ss_df['A1']!=self.ss_df['ALT'], 'PRScs'] *= -1


        ### fill NA
        self.ss_df.fillna(0, inplace=True)

        print('\n\n###### Complete ######\n\n')
        

    def _get_ct_beta(self, col_name, snp_list, threshold=1):
        self.ss_df[col_name] = self.ss_df['BETA']
        self.ss_df.loc[~self.ss_df['ID'].isin(snp_list), col_name] = 0
        self.ss_df.loc[self.ss_df['P'] > threshold, col_name] = 0

    
    def _get_ldpred2_beta(self, prefix):
        df = pyreadr.read_r('{}.snp.rds'.format(prefix))[None]
        df['beta_modify'] = pyreadr.read_r('{}.beta.rds'.format(prefix))[None]
        df['_NUM_ID_'] = df['_NUM_ID_'] - 1
        self.ss_df['LDpred2'] = 0
        idx = self.ss_df.iloc[df['_NUM_ID_']].index
        self.ss_df.loc[idx, 'LDpred2'] = df['beta_modify'].to_numpy()


    def _get_prscs_beta(self, file):
        df = pd.read_csv(file, sep='\s+', names=['CHR', 'ID', 'POS', 'A1', 'A2', 'PRScs'])
        self.ss_df = self.ss_df.merge(df[['ID', 'PRScs']], on='ID', how='left')
        self.ss_df['PRScs'].fillna(0, inplace=True)



# build the dataframe for PRS results
class PRSResults():
    def __init__(self, bfile_prefix, pred_prefix, method):
        self.bfile_prefix = bfile_prefix
        self.pred_prefix = pred_prefix
        self.method = method
        self.df = pd.read_csv('{}.fam'.format(bfile_prefix), sep='\s+', names=['FID', 'IID', 'father', 'mother', 'sex', 'phenotype'])
        self.df = self.df[['FID', 'IID', 'phenotype']]
        if self.method == 'clf':
            self.df['phenotype'] = self.df['phenotype'].apply(lambda x: x - 1 if (x == 1) | (x == 2) else x) # control = 0, case = 1
        self.df['phenotype'].replace(-9, np.nan, inplace=True)


    def __call__(self):
        print('\n\n###### Building PRS Prediction Dataframe ######\n\n')
        ### prediction dataframe
        print('Loading prediction dataframe ...')
        for file in glob(f"{self.pred_prefix}.*.profile"):
            #file = '{}.{}.profile'.format(self.pred_prefix, tool)
            tool = file.replace(f"{self.pred_prefix}.","").replace(".profile","")
            self.df = self._append_score(self.df, file, tool)
        # GenEpi
        file = '{}.GenEpi.csv'.format(self.pred_prefix)
        self.df = self._append_score(self.df, file, 'GenEpi', sep=',', prev_col='score')
        print('Available algorithms:', ', '.join(self.df.columns[3:]))

        ### checking NA
        # fill NA with minimum
        # if NA > 10%, drop the algorithm
        print('Checking NA ...')
        num = self.df.shape[0]
        for tool in self.df.columns[3:]:
            na_count = self.df[tool].isna().sum()
            print('NA count of {} = {} / {}, {:.2f}%'.format(tool, na_count, num, na_count/num*100))
            if na_count/num > 0.1:
                self.df = self.df.drop(columns=[tool])
                print('Drop {} because more than 10\% are NA'.format(tool))
            else:
                minimum = self.df[tool].min()
                self.df[tool] = self.df[tool].fillna(minimum)
                print('Fill {} with minimum ({})'.format(tool, minimum))

        return self.df
    

    def _append_score(self, df, file, post_col, sep='\s+', prev_col='SCORESUM'):
        if not os.path.isfile(file):
            return df
        pred_df = pd.read_csv(file, sep=sep)
        pred_df = pred_df[['FID', 'IID', prev_col]]
        pred_df = pred_df.rename(columns={prev_col: post_col})
        df = df.merge(pred_df, on=['FID', 'IID'], how='left')
        return df



# build the dataframe for PRS+cov results
class CovResults():
    def __init__(self, prs_df, cov_df, method):
        print('\n\n###### Building PRS+Covariates Prediction Dataframe ######\n\n')

        # checking method
        if method not in ['clf', 'reg']:
            print('Method must be clf or reg')
            sys.exit()
        self.method = method

        # getting predictors and covariates
        self.tools = list(prs_df.columns)[3:]
        self.covs = list(cov_df.columns)[2:]
        print('{} PRS algorithms'.format(len(self.tools)))
        print('PRS algorithms:', ', '.join(self.tools))
        print('{} covariates'.format(len(self.covs)))
        print('Covariates:', ', '.join(self.covs))

        # building dataframe
        self.df = prs_df.merge(cov_df, on=['FID', 'IID'], how='inner')


    def Train(self):
        print('\n\n###### Training Regression Models with Covariates ######\n\n')

        # define functions
        scaler = StandardScaler()
        if self.method == 'clf':
            reg = LogisticRegression()
        else:
            reg = LinearRegression()

        # prediction dataframe
        self.df = self.df.dropna(subset=['phenotype'], axis=0) # drop NA for building regression model
        pred_df = self.df.loc[:, ['FID', 'IID', 'phenotype']]
        y = self.df['phenotype'].to_numpy()

        # apply regression
        model_dict = dict()
        tools = ['cov'] + self.tools
        for tool in tools:
            model_dict[tool] = dict() # columns, numerical_columns, scaler_mean, scaler_scale, reg_coef, reg_intercept
            
            # columns
            if tool == 'cov':
                print('Training covariates only ...')
                cols = self.covs
            else:
                print('Training {} with covariates ...'.format(tool))
                cols = [tool] + self.covs
            model_dict[tool]['columns'] = cols

            # numerical columns
            numerical_cols = list()
            for i in range(len(cols)):
                if self.df[cols[i]].dropna().unique().shape[0] > 2:
                    numerical_cols.append(i)
            model_dict[tool]['numerical_columns'] = numerical_cols

            # scaler: only applied on numerical features
            x = self.df[cols].to_numpy()
            scaler_mean = 0 
            scaler_scale = 0
            if len(numerical_cols) > 0:
                x[:, numerical_cols] = scaler.fit_transform(x[:, numerical_cols])
                scaler_mean = scaler.mean_.tolist()
                scaler_scale = scaler.scale_.tolist()
                
            model_dict[tool]['scaler_mean'] = scaler_mean
            model_dict[tool]['scaler_scale'] = scaler_scale
            
            # regression
            x = np.nan_to_num(x) # fill NA with 0 (mean of scaler = 0)
            reg.fit(x, y)
            if self.method == 'clf':
                pred = reg.predict_proba(x)[:,1]
                model_dict[tool]['reg_coef'] = reg.coef_.tolist()[0]
                model_dict[tool]['reg_intercept'] = reg.intercept_[0]
            else:
                pred = reg.predict(x)
                model_dict[tool]['reg_coef'] = reg.coef_.tolist()
                model_dict[tool]['reg_intercept'] = reg.intercept_
            

            # prediction
            pred_df[tool] = pred
        
        print('\n\n###### Complete ######\n\n')
        return pred_df, model_dict


    def Test(self, model_path):
        print('\n\n###### Predicting Risk Scores with Covariates ######\n\n')

        # load model
        model_dict = json.load(open(model_path, 'r'))
        model_tools = list(model_dict.keys())

        # predict
        pred_df = self.df.loc[:, ['FID', 'IID', 'phenotype']]
        y = self.df['phenotype'].to_numpy()
        tools = ['cov'] + self.tools
        for tool in tools:
            if not tool in model_tools:
                print('Skip {} because of no available model'.format(tool))
                continue
            if tool == 'cov':
                print('Testing covariates only ...')
            else:
                print('Testing {} with covariates ...'.format(tool))

            # model
            cols = model_dict[tool]['columns']
            numerical_cols = model_dict[tool]['numerical_columns']
            scaler_mean = np.array(model_dict[tool]['scaler_mean'])
            scaler_scale = np.array(model_dict[tool]['scaler_scale'])
            reg_coef = np.array(model_dict[tool]['reg_coef'])
            reg_intercept = model_dict[tool]['reg_intercept']
            
            # x
            x = self.df[cols].to_numpy()
            if len(numerical_cols) > 0:
                x[:, numerical_cols] = (x[:, numerical_cols] - scaler_mean) / scaler_scale

            # regression
            x = np.nan_to_num(x)
            if self.method == 'clf':
                pred = 1 / (1 + np.exp(-(np.matmul(x, reg_coef) + reg_intercept)))
            else:
                pred = np.matmul(x, reg_coef) + reg_intercept

            # record
            pred_df[tool] = pred
        
        print('\n\n###### Complete ######\n\n')
        return pred_df



# build population reference, including rank and histogram
class CohortRef():
    def __init__(self, pred_df, rank_ref_file=None):
        print('\n\n###### Building Cohort Reference ######\n\n')
        ### prediction dataframe
        self.df = pred_df

        ### ranking dataframe
        if rank_ref_file:
            print('Loading provided ranking reference ...')
            self.rank_ref_df = pd.read_csv(rank_ref_file, index_col=0)
            if self.rank_ref_df.index.name is not None:
                self.rank_ref_df = pd.read_csv(rank_ref_file)
        else:
            print('Building ranking reference ...')
            self.rank_ref_df = self._build_rank_ref(self.df)
            print('Building histogram...')
            self.hist_ref_df = self._build_histogram_ref(self.df, self.rank_ref_df)


    def __call__(self):
        print('Mapping ranking score ...')
        self.rank_df = self._map_rank(self.df, self.rank_ref_df)
        return self.rank_df


    def _build_rank_ref(self, df):
        rank_list = list(range(100)) + [99.5, 99.7, 99.9, 100]
        rank_df = pd.DataFrame(index=rank_list)
        for tool in df.columns[3:]:
            rank_score = np.percentile(df[tool], rank_list)
            rank_df[tool] = rank_score
        return rank_df


    def _map_rank(self, df, rank_ref_df):
        rank_df = df.loc[:, df.columns[:3]]
        for tool in df.columns[3:].tolist():
            if tool not in rank_ref_df.columns.tolist():
                continue
            interp_rank = np.interp(df[tool], rank_ref_df[tool], rank_ref_df.index.tolist())
            rank_df[tool] = interp_rank
        return rank_df


    def _build_histogram_ref(self, df, rank_ref_df):
        hist_df = pd.DataFrame()
        for tool in rank_ref_df.columns:
            score_min = rank_ref_df.loc[0, tool]
            score_max = rank_ref_df.loc[100, tool]
            hist, bin_edges = np.histogram(df[tool], bins=100, range=(score_min, score_max), density=True)
            hist_df[tool] = hist
        return hist_df



# analyze the prediction
class Analysis():
    def __init__(self, df, method, outdir):
        print('\n\n###### Preparing for Analysis ######\n\n')
        # preprocessing dataframe
        self.df = df
        self.tools = list(df.columns)[3:] # ['FID', 'IID', 'phenotype', ALGO_1, ALGO_2, ...]
        total_num = self.df.shape[0]
        self.df = self.df.dropna(subset=['phenotype'], axis=0)
        dropna_num = self.df.shape[0]
        print("Drop {} samples with no phenotype".format(total_num - dropna_num))

        # checking method
        if method not in ['clf', 'reg']:
            print('Method must be clf or reg')
            sys.exit()
        self.method = method

        # checking outdir
        self.outdir = outdir
        if not os.path.isdir(self.outdir):
            os.mkdir(self.outdir)


    def __call__(self, percentile_num=10, fontsize=8, linewidth=1, figsize=3, dpi=200):
        print('\n\n###### Analyzing PRS Performance ######\n\n')
        self.percentile_num = percentile_num
        self.fontsize = fontsize
        self.linewidth = linewidth
        self.figsize = figsize
        self.dpi = dpi

        plt.rcParams.update({'font.size': fontsize})

        if self.method == 'clf':
            print('Analyzing Classification ...')
            self.AnaCLF()
        elif self.method == 'reg':
            print('Analyzing Regression ...')
            self.AnaREG()
        
        print('\n\n###### Complete ######\n\n')


    def AnaCLF(self):
        ##### ROC, PRC
        print('Calculating metrics auROC and auPRC ...')
        self.results = dict()
        self.results['ROC'] = self._roc_plot(self.df, self.tools, title='ROC', figfile='{}/ROC.png'.format(self.outdir))
        self.results['PRC'] = self._prc_plot(self.df, self.tools, title='PRC', figfile='{}/PRC.png'.format(self.outdir))
        json.dump(self.results, open('{}/performance.json'.format(self.outdir), 'w'))

        ##### distribution
        print('Plotting prediction distribution ...')
        self._dist_plot(self.df, self.tools, figfile='{}/distribution.png'.format(self.outdir))

        ##### percentile
        print('Calculating percentile distribution ...')
        self.percentile_df = self._percentile(self.df, self.tools, n=self.percentile_num, clf=True)
        self.percentile_df.to_csv('{}/percentile.csv'.format(self.outdir), index=False)
        self._or_percentile_plot(self.percentile_df, title='percentile of OR', figfile='{}/ORpercentile.png'.format(self.outdir))


    def AnaREG(self):
        ##### Pearson and Spearman
        print('Calculating metrics Pearson and Spearman correlation ...')
        self.results = dict()
        y = self.df['phenotype']
        ##self.results['R2'] = dict()
        self.results['Pearson'] = dict()
        self.results['Spearman'] = dict()

        for tool in self.tools:
            pred = self.df[tool]
            # R square
            ##r2_score = metrics.r2_score(y, pred)
            ##self.results['R2'][tool] = r2_score
            # Pearson correlation
            pearson, pvalue = stats.pearsonr(y, pred)
            self.results['Pearson'][tool] = [pearson, pvalue]
            # Spearman correlation
            spearman, pvalue = stats.spearmanr(y, pred)
            self.results['Spearman'][tool] = [spearman, pvalue]
        json.dump(self.results, open('{}/performance.json'.format(self.outdir), 'w'))
        
        # plot
        self._bar_plot(self.results, 'Pearson', figfile='{}/pearson.png'.format(self.outdir))
        self._bar_plot(self.results, 'Spearman', figfile='{}/spearman.png'.format(self.outdir))

        ##### percentile
        print('Calculating percentile distribution ...')
        self.percentile_df = self._percentile(self.df, self.tools, n=self.percentile_num, clf=False)
        self.percentile_df.to_csv('{}/percentile.csv'.format(self.outdir), index=False)
        self._percentile_plot(self.percentile_df, title='percentile', figfile='{}/percentile.png'.format(self.outdir))


    def _roc_plot(self, df, tools, title=None, figfile=None):
        # figure
        fig, ax = plt.subplots(1, 1, figsize=(self.figsize*1.25, self.figsize), dpi=self.dpi)

        # calculate AUC
        roc_res={}
        results = dict()
        for tool in tools:
            pred = df[tool]
            fpr, tpr, _ = metrics.roc_curve(df['phenotype'], pred)
            auc = metrics.auc(fpr, tpr)
            results[tool] = auc
            ax.plot(fpr, tpr, label='%s,AUC=%.3f'%(tool, auc), linewidth=self.linewidth)
            roc_res[tool] = {'fpr':fpr.tolist(),'tpr':tpr.tolist(),'auc':auc}

        # output fpr tpr to roc.json
        json.dump(roc_res, open('{}/roc.json'.format(self.outdir), 'w'))
        # adjust label and legend
        ax.legend(loc='lower right', handletextpad=0.2, borderpad=0.2)
        ax.plot([0, 1], [0, 1], '--', color='black')
        ax.set_xlim([0, 1])
        ax.set_ylim([0, 1])
        ax.set_title(title, fontsize=self.fontsize+2)
        ax.set_ylabel('True Positive Rate')
        ax.set_xlabel('False Positive Rate')
        ax.tick_params(axis='both', which='major')
        
        # savefig
        fig.tight_layout()
        if figfile:
            fig.savefig(figfile)
        plt.close()

        return results


    def _prc_plot(self, df, tools, title=None, figfile=None):
        # figure
        fig, ax = plt.subplots(1, 1, figsize=(self.figsize*1.25, self.figsize), dpi=self.dpi)

        # calculate AUC
        results = dict()
        for tool in tools:
            pred = df[tool]
            precision, recall, _ = metrics.precision_recall_curve(df['phenotype'], pred)
            auc = metrics.average_precision_score(df['phenotype'], pred)
            results[tool] = auc
            ax.plot(recall, precision, label='%s,AP=%.3f'%(tool, auc), linewidth=self.linewidth)
        
        # adjust label and legend
        ax.legend(loc='best', handletextpad=0.2, borderpad=0.2)
        ax.set_xlim([0, 1])
        ax.set_ylim([0, 1])
        ax.set_title(title, fontsize=self.fontsize+2)
        ax.set_xlabel('Recall')
        ax.set_ylabel('Precision')
        ax.tick_params(axis='both', which='major')

        # savefig
        fig.tight_layout()
        if figfile:
            fig.savefig(figfile)
        plt.close()
        
        return results

    
    def _dist_plot(self, df, tools, figfile=None):
        # figure
        if (len(tools)==4) or (len(tools)==2):
            ncol = 2
        else:
            ncol = 3
        nrow = int(np.ceil(len(tools)/ncol))
        fig, ax = plt.subplots(nrow, ncol, figsize=(ncol*self.figsize, (nrow/1.5)*self.figsize), dpi=self.dpi)

        # distribution
        for i in range(len(tools)):
            if nrow == 1:
                cur_ax = ax[i]
            else:
                cur_ax = ax[i//ncol][i%ncol]
            tool = tools[i]
            sns.histplot(data=df, x=tool, hue='phenotype', stat='density',
                         common_norm=False, element='step', ax=cur_ax, bins=30, legend=False)
            #_ = cur_ax.set_title(tool)
            cur_ax.axvline(df.loc[df['phenotype']==0, tool].mean(), color='#3A76AD', ls='--', lw=0.7)
            cur_ax.axvline(df.loc[df['phenotype']==1, tool].mean(), color='#EF8637', ls='--', lw=0.7)

        # legend
        patches = [
            mpatches.Patch(color='#CADBEA', label='control'),
            mpatches.Patch(color='#FADFC5', label='case')
        ]
        fig.legend(handles=patches, loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.05))

        # savefig
        fig.tight_layout()
        if figfile:
            fig.savefig(figfile, bbox_inches='tight')
        plt.close()

    
    def _percentile(self, df, tools, n=10, clf=True):
        # percentile list
        percentile_li = list()
        interval = int(df.shape[0] // n)
        idx = [(i+1)*(100/n) for i in range(n)]
        for i in idx:
            percentile_li += [i,] * interval
            if i == idx[-1]:
                percentile_li += [i,] * (df.shape[0] - interval*n)
        percentile_li = percentile_li[::-1]
    
        # assign percentile to each tool
        result_df = pd.DataFrame()
        for tool in tools:
            df = df.sort_values(by=tool, ascending=False)
            temp_df = df.loc[:, ['phenotype', tool]]
            temp_df['tool'] = tool
            temp_df['percentile'] = percentile_li
            temp_df = temp_df.rename(columns={tool: 'PRS'})
            result_df = pd.concat([result_df, temp_df])
    
        # regression
        if not clf:
            return result_df
    
        # calculate OR for classification
        or_df = pd.DataFrame(columns=['tool', 'percentile', 'OR', 'ci_upper', 'ci_lower', 'pos_num', 'neg_num'])
        percentile_li = sorted(set(percentile_li))
        base_neg, base_pos = [df[df['phenotype']==i].shape[0]//n for i in [0, 1]]
        
        for tool in tools:
            # OR of each percentile
            for percentile in percentile_li:
                temp_df = result_df[(result_df['tool']==tool) & (result_df['percentile']==percentile)]
                neg, pos = [temp_df[temp_df['phenotype']==i].shape[0] for i in [0, 1]]
                OR, ci_upper, ci_lower = self._cal_or(pos, neg, base_pos, base_neg)
                or_df.loc[len(or_df)] = [tool, percentile, OR, ci_upper, ci_lower, pos, neg]
        
        return or_df


    def _cal_or(self, a, b, c, d):
        a += 0.5
        b += 0.5
        c += 0.5
        d += 0.5
        OR = (a*d) / (b*c)
        se = np.sqrt(1/a + 1/b + 1/c + 1/d) # SE(log(OR))
        ci_upper = np.exp(np.log(OR) + 1.96*se)
        ci_lower = np.exp(np.log(OR) - 1.96*se)
        return OR, ci_upper, ci_lower


    def _or_percentile_plot(self, percentile_df, title=None, figfile=None):
        # figure
        fig, ax = plt.subplots(1, 1, figsize=(2*self.figsize, self.figsize), dpi=self.dpi)

        # point plot
        sns.pointplot(x='percentile', y='OR', hue='tool', data=percentile_df,
                      join=False, dodge=0.3, scale=0.4, ci=None, ax=ax)

        # error bar
        ## get colors
        h, l = ax.get_legend_handles_labels()
        color_map = {l[i]: h[i].get_facecolor() for i in range(len(h))}
        ecolor = list(map(lambda x: color_map[x], percentile_df['tool']))
        ecolor = np.reshape(np.array(ecolor), (-1,4))

        ## get point position
        x_coords = []
        y_coords = []
        for point_pair in ax.collections:
            for x, y in point_pair.get_offsets():
                x_coords.append(x)
                y_coords.append(y)

        ## error range
        bound_arr = percentile_df[['ci_lower', 'ci_upper']].to_numpy()
        base_arr = np.reshape(percentile_df['OR'].to_numpy(), (-1,1))
        yerr = np.abs(bound_arr - base_arr).T

        ## plot
        ax.errorbar(x=x_coords, y=y_coords, yerr=yerr, fmt=' ', ecolor=ecolor, elinewidth=self.linewidth/4)

        # savefig
        ax.legend(loc='upper left', ncol=3)
        ax.set_title(title)
        fig.tight_layout()
        if figfile:
            fig.savefig(figfile)
        plt.close()


    def _percentile_plot(self, percentile_df, title=None, figfile=None):
        # plot
        fig, ax = plt.subplots(1, 1, figsize=(2*self.figsize, self.figsize), dpi=self.dpi)
        sns.pointplot(x='percentile', y='phenotype', hue='tool', data=percentile_df,
                      join=False, dodge=0.3, scale=0.4, errwidth=self.linewidth/4, ax=ax)
        
        # savefig
        ax.legend(loc='upper left', ncol=3)
        ax.set_title(title)
        fig.tight_layout()
        if figfile:
            fig.savefig(figfile)
        plt.close()


    def _bar_plot(self, result_dict, metric, title=None, figfile=None):
        # df
        df = pd.DataFrame(result_dict[metric]).iloc[0]
        
        # plot
        fig, ax = plt.subplots(1, 1, figsize=(self.figsize, self.figsize), dpi=self.dpi)
        sns.barplot(x=df.index, y=df, ax=ax)
        ax.set_title(title, fontsize=self.fontsize+2)
        ax.set_ylabel(metric)
        
        # savefig
        fig.tight_layout()
        if figfile:
            fig.savefig(figfile)
        plt.close()
