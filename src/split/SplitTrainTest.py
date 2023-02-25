#!/usr/bin/python3
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
import seaborn as sns

def ArgumentParser():
    parser = argparse.ArgumentParser(prog='SplitTrainTest')
    parser.add_argument('--fam_file', required=True, help='the fam file')
    parser.add_argument('--testing_ratio', required=False, help='the testing ratio', type=float, default=0.1)
    parser.add_argument('--random', required=False, default=False, action='store_true', help='randomly split')
    parser.add_argument('--dropna', required=False, default=False, action='store_true', help='drop missing phenotype')
    parser.add_argument('--method', required=False, default='clf', help='clf (classification) or reg (regression), default=clf')
    parser.add_argument('--out_basename', required=True, help='the output basename')
    parser.add_argument('--train_suffix', required=False, default='train', help='the suffix of training')
    parser.add_argument('--test_suffix', required=False, default='test', help='the suffix of testing')
    return parser


def GetTestList(fam_df, testing_ratio):
    test_idx = []
    for group,da in fam_df.groupby("GROUP"):
        ll = int(len(da.index)*testing_ratio)

        if ll > 0:
            da = da.sample(n = ll)
            test_idx = test_idx + da.index.values.tolist()
    test_idx.sort()
    return test_idx


if __name__=='__main__':
    # arguments
    args = ArgumentParser().parse_args()
    fam_file = args.fam_file
    testing_ratio = args.testing_ratio
    random_mode = args.random
    method = args.method
    dropna = args.dropna
    out_basename = args.out_basename
    train_suffix = args.train_suffix
    test_suffix = args.test_suffix

    print('Method: {}, DropNA: {}, Random: {}'.format(method, dropna, random_mode))

    # laod fam file
    fam_df = pd.read_csv(fam_file, sep='\s+', names=['FID', 'IID', 'father', 'mother', 'sex', 'phenotype'])
    na_num = fam_df[fam_df['phenotype'] == -9].shape[0]
    print('The ratio of missing phenotype: {0:.2f}'.format(na_num / fam_df.shape[0]))

    # drop NA or not
    if dropna:
        fam_df = fam_df[(fam_df['phenotype'] != -9) & ~fam_df.phenotype.isna() ]
        print('Drop NA')
        
        
    np.random.seed(0)
    #fam_df = fam_df.sample(n=len(fam_df.index)).reset_index(drop=True)
    

    # method: clf
    if method == 'clf':
        fam_df["GROUP"] = fam_df.phenotype
        
        test_idx = GetTestList(fam_df, testing_ratio)
        
        # add tags of train and test
        fam_df['tag'] = 'train'
        fam_df.loc[test_idx, 'tag'] = 'test'
        fam_df = fam_df.drop(columns = ["GROUP"])

        print(pd.crosstab(fam_df.phenotype, fam_df.tag))
    # method: reg
    else:
        # Add split by rank
        # group by rank
        windows_size = int(1/testing_ratio)
        fam_df["RANK"] = fam_df.phenotype.rank()
        fam_df["GROUP"] = fam_df.RANK % windows_size 

        # get test_ratio from each group
        test_idx = GetTestList(fam_df, testing_ratio)

        # add tags of train and test
        fam_df['tag'] = 'train'
        fam_df.loc[test_idx, 'tag'] = 'test'
        fam_df = fam_df.drop(columns = ["GROUP", "RANK"])

        # distribution
        fig, ax = plt.subplots(1, 1, figsize=(6, 4), dpi=200)
        sns.histplot(data=fam_df[fam_df['phenotype'] != -9], hue='tag', x='phenotype', ax=ax,
                     stat='probability', common_norm=False, element='step')
        fig.tight_layout()
        fig.savefig('{}.dist.png'.format(out_basename))
        print('Distribution plot of training and testing at {}.dist.png'.format(out_basename))
        

    # save
    fam_df[fam_df['tag']=='train'][['FID', 'IID']].to_csv('{}.{}.list'.format(out_basename, train_suffix), header=False, index=False, sep='\t')
    fam_df[fam_df['tag']=='test'][['FID', 'IID']].to_csv('{}.{}.list'.format(out_basename, test_suffix), header=False, index=False, sep='\t')
    
    
'''
test code
python3 /yilun/prs-algo/split/SplitTrainTest.py \
    --fam_file /mnt/prsdata/Test/Data/TWB2_LDL_extSumStats/TWB2_LDL.fam \
    --testing_ratio 0.1 \
    --random \
    --dropna \
    --method "reg" \
    --out_basename /tmp/TWB2_LDL \
    --train_suffix "train" \
    --test_suffix "test"    


python3 /yilun/prs-algo/split/SplitTrainTest.py \
    --fam_file /mnt/prsdata/Test/Data/DEMO_CLF/demo_hg38.fam \
    --testing_ratio 0.1 \
    --random \
    --dropna \
    --method "clf" \
    --out_basename /tmp/demo_clf \
    --train_suffix "train" \
    --test_suffix "test"    



'''