# -*- coding: utf-8 -*-
"""
Created on Nov 2020

@author: Chester (Yu-Chuan Chang)
"""

""""""""""""""""""""""""""""""
# import libraries
""""""""""""""""""""""""""""""
import argparse
import os
import sys
import numpy as np
import joblib
import re

import numpy as np
import scipy as sp
import sklearn.metrics as skMetric
import scipy.stats as stats
from scipy.optimize import curve_fit
from scipy.stats import norm

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

""""""""""""""""""""""""""""""
# define functions 
""""""""""""""""""""""""""""""
def ArgumentsParser():
    ### define arguments
    str_description = ''
    'This script is a tool for GenEpi prediction'
    parser = argparse.ArgumentParser(prog='predictor', description=str_description)
    
    ### define arguments for I/O
    parser.add_argument("-g", required=True, help="filename of the input .gen file")
    parser.add_argument("-p", required=False, help="filename of the input phenotype")
    parser.add_argument("-m", required=True, help="filename of the predicting model")
    parser.add_argument("-f", required=True, help="filename of the feature file")    
    parser.add_argument("-o", required=False, help="output file path")
    
    return parser

def FeatureGenerator(str_inputFileName_genotype, str_inputFileName_feature, str_outputFilePath = ""):
    ### set default output path
    if str_outputFilePath == "":
        str_outputFilePath = os.path.dirname(str_inputFileName_genotype) + "/predictedResult/"
    ### if output folder doesn't exist then create it
    if not os.path.exists(str_outputFilePath):
        os.makedirs(str_outputFilePath)

    ### get all selected snp ids
    list_feature_rsid_all = []
    with open(str_inputFileName_feature, "r") as file_inputFile:
        ### grep the header
        list_rsids = file_inputFile.readline().strip().split(",")
        for rsid in list_rsids:
            list_feature_rsid_all.append(rsid)
    ### get unique selected snp ids
    dict_feature_rsid_unique = {}
    for item in list_feature_rsid_all:
        for subitem in item.split("*"): # modify
            subitem_rsid=re.sub(r"\_\S*\.?\S*", "", subitem)
            if subitem_rsid not in dict_feature_rsid_unique:
                dict_feature_rsid_unique[subitem_rsid] = 1
    
    ### extract selected snp from genotype file
    list_inputFile_genotype = []
    with open(str_inputFileName_genotype, 'r') as file_inputFile:
        for line in file_inputFile:
            list_thisSnp = line.strip().split(" ")
            if list_thisSnp[1] in dict_feature_rsid_unique.keys():
                list_inputFile_genotype.append(line)
    
    ### count lines of input files
    int_num_genotype = len(list_inputFile_genotype)
    int_num_phenotype = 0
    
    ### get phenotype number
    with open(str_inputFileName_genotype, 'r') as file_inputFile:
        list_thisSnp = file_inputFile.readline().strip().split(" ")
        int_num_phenotype = int((len(list_thisSnp) - 5) / 3)

    ### get genotype file
    list_genotype = [[] for x in range(int_num_phenotype)]
    list_genotype_rsid = []
    with open(str_inputFileName_genotype, 'r') as file_inputFile:
        for line in file_inputFile:
            list_thisSnp = line.strip().split(" ")
            if list_thisSnp[1] not in dict_feature_rsid_unique.keys():
                continue
            np_this_genotype = np.empty([int_num_phenotype, 3], dtype='int8')
            for idx_subject in range(0, int_num_phenotype):
                list_allelType = [0, 0, 0]
                if list_thisSnp[idx_subject * 3 + 5 : idx_subject * 3 + 8] != ['0', '0', '0']:
                    list_allelType[np.argmax(list_thisSnp[idx_subject * 3 + 5 : idx_subject * 3 + 8])] = 1
                list_genotype[idx_subject].extend(list_allelType)
            list_genotype_rsid.append(list_thisSnp[1] + "_" + list_thisSnp[3] + "." + list_thisSnp[3])
            list_genotype_rsid.append(list_thisSnp[1] + "_" + list_thisSnp[3] + "." + list_thisSnp[4])
            list_genotype_rsid.append(list_thisSnp[1] + "_" + list_thisSnp[4] + "." + list_thisSnp[4])
    np_genotype = np.array(list_genotype, dtype=np.int8)
    np_genotype_rsid = np.array(list_genotype_rsid)
    
    ### generate feature
    np_feature = np.empty([int_num_phenotype, len(list_feature_rsid_all)], dtype='int')
    for idx_feature in range(len(list_feature_rsid_all)):
        list_feature_rsid = list_feature_rsid_all[idx_feature].split("*")
        if len(list_feature_rsid) == 1:
            #np_feature[:, idx_feature] = np_genotype[:, int(np.argwhere(np_genotype_rsid == list_feature_rsid[0]))]
            list_feature_rsid_0 = list_feature_rsid[0]
            if len(np.argwhere(np_genotype_rsid == list_feature_rsid[0])) == 0:
                split_snp=re.split('_|\.',list_feature_rsid_0)
                list_feature_rsid_0=split_snp[0]+'_'+split_snp[2]+'.'+split_snp[1]
            np_feature[:, idx_feature] = np_genotype[:, int(np.argwhere(np_genotype_rsid == "".join(list_feature_rsid_0)))]
        else:
            #np_feature[:, idx_feature] = np.multiply(np_genotype[:, int(np.argwhere(np_genotype_rsid == list_feature_rsid[0]))], np_genotype[:, int(np.argwhere(np_genotype_rsid == list_feature_rsid[1]))])
            list_feature_rsid_0 = list_feature_rsid[0]
            if len(np.argwhere(np_genotype_rsid == list_feature_rsid[0])) == 0:
                split_snp=re.split('_|\.',list_feature_rsid_0)
                list_feature_rsid_0=split_snp[0]+'_'+split_snp[2]+'.'+split_snp[1]
            list_feature_rsid_1 = list_feature_rsid[1]
            if len(np.argwhere(np_genotype_rsid == list_feature_rsid[1])) == 0:
                split_snp=re.split('_|\.',list_feature_rsid_1)
                list_feature_rsid_1=split_snp[0]+'_'+split_snp[2]+'.'+split_snp[1]
            np_feature[:, idx_feature] = np.multiply(np_genotype[:, int(np.argwhere(np_genotype_rsid == "".join(list_feature_rsid_0)))], np_genotype[:, int(np.argwhere(np_genotype_rsid == "".join(list_feature_rsid_1)))])
    
    ### output feature
    with open(os.path.join(str_outputFilePath, "Feature.csv"), "w") as file_outputFile:
        file_outputFile.writelines(",".join(list_feature_rsid_all) + "\n")
        for idx_subject in range(int_num_phenotype):
            file_outputFile.writelines(",".join(np_feature[idx_subject, :].astype(str)) + "\n")
    
    return np_feature

def IsolatedDataPredictor(str_inputFileName_genotype, str_inputFileName_model, str_inputFileName_feature, str_outputFilePath = "", str_mode = "c"):
    ### set default output path
    if str_outputFilePath == "":
        str_outputFilePath = os.path.dirname(str_inputFileName_genotype) + "/"
    ### if output folder doesn't exist then create it
    if not os.path.exists(str_outputFilePath):
        os.makedirs(str_outputFilePath)
    
    estimator = joblib.load(str_inputFileName_model)
    np_genotype = FeatureGenerator(str_inputFileName_genotype, str_inputFileName_feature, str_outputFilePath)
    
    if str_mode == "c":
        list_predict = []
        list_predict_proba = []
        this_predict = estimator.predict(np_genotype)
        this_predict_proba = estimator.predict_proba(np_genotype)
        with open(os.path.join(str_outputFilePath, "Prediction.csv"), "w") as file_outputFile:
            file_outputFile.writelines("Prediction,PRS" + "\n")
            for idx_predict, idx_predict_proba in zip(this_predict, this_predict_proba):
                file_outputFile.writelines(str(idx_predict) + "," + "{:0.4f}".format(idx_predict_proba[1]) + "\n")
                #print(str(idx_predict) + " - " + "{:0.4f}".format(idx_predict_proba[1])) # modify: mute
                list_predict.append(idx_predict)
                list_predict_proba.append(idx_predict_proba)

        return list_predict, list_predict_proba
    
    else:
        list_predict = []
        this_predict = estimator.predict(np_genotype)
        with open(os.path.join(str_outputFilePath, "Prediction.csv"), "w") as file_outputFile:
            file_outputFile.writelines("Prediction" + "\n")
            for idx_predict in zip(this_predict):
                file_outputFile.writelines(str(idx_predict[0]) + "\n")
                #print(str(idx_predict[0])) # modify: mute
                list_predict.append(idx_predict[0])

        return list_predict, None

def gaussian(x, mean, amplitude, standard_deviation):
    return amplitude * np.exp( - ((x - mean) / standard_deviation) ** 2)

def fsigmoid(x, a, b):
    float_return = 1.0/(1.0+np.exp(-a*(x-b)))
    return float_return

def PlotPolygenicScore(list_target, list_predict, list_proba, str_outputFilePath="", str_label=""):
    """

    Plot figure for polygenic score, including group distribution and prevalence to PGS

    Args:
        list_target (list): A list containing the target of each samples
        list_predict (list): A list containing the predition value of each samples
        list_proba (list): A list containing the predition probability of each samples
        str_outputFilePath (str): File path of output file
        str_label (str): The label of the output plots

    Returns:
        None
    
    """

    float_f1Score = skMetric.f1_score(list_target, list_predict)

    #-------------------------
    # group distribution
    #-------------------------
    pd_pgs = pd.concat([pd.DataFrame(list_target), pd.DataFrame(list_predict), pd.DataFrame(np.array(list_proba)[:,1])], axis=1)
    pd_pgs.columns = ['target', 'predict', 'proba']

    int_bin = 25
    plt.figure(figsize=(5,5))
    
    # plot case
    pd_case = pd_pgs[pd_pgs.target == 1.0]
    plt.hist(pd_case['proba'], bins=int_bin, label='Case', color="#e68fac", weights=np.ones_like(pd_case['proba'])/float(len(pd_case['proba'])))
    bin_heights, bin_borders = np.histogram(pd_case['proba'], bins=int_bin)
    bin_heights = bin_heights / float(len(pd_case['proba']))
    bin_widths = np.diff(bin_borders)
    bin_centers = bin_borders[:-1] + bin_widths / 2
    popt, _ = curve_fit(gaussian, bin_centers, bin_heights, maxfev=100000000)
    x_interval_for_fit = np.linspace(bin_borders[0], bin_borders[-1], 10000)
    plt.plot(x_interval_for_fit, gaussian(x_interval_for_fit, *popt), c="#b3446c")

    # plot control
    pd_control = pd_pgs[pd_pgs.target == 0.0]
    plt.hist(pd_control['proba'], bins=int_bin, label='Control', color="#4997d0", weights=np.ones_like(pd_control['proba'])/float(len(pd_control['proba'])))
    bin_heights, bin_borders = np.histogram(pd_control['proba'], bins=int_bin)
    bin_heights = bin_heights / float(len(pd_control['proba']))
    bin_widths = np.diff(bin_borders)
    bin_centers = bin_borders[:-1] + bin_widths / 2
    popt, _ = curve_fit(gaussian, bin_centers, bin_heights, maxfev=100000000)
    x_interval_for_fit = np.linspace(bin_borders[0], bin_borders[-1], 10000)
    plt.plot(x_interval_for_fit, gaussian(x_interval_for_fit, *popt), c="#00416a")

    # plot formatting
    str_method = "GenEpi"
    plt.legend(prop={'size': 12})
    plt.title(str_method + ' Predicting F1 Score: ' + "%.4f" % float_f1Score + ' ')
    plt.xlim(0, 1)
    plt.ylim(0, 0.5)
    plt.xlabel('Polygenic Score')
    plt.ylabel('Fraction of samples by group')
    plt.savefig(os.path.join(str_outputFilePath, str("GenEpi_PGS_" + str_label + ".png")), dpi=300)
    plt.close('all')

    #-------------------------
    # prevalence to PGS
    #-------------------------
    int_step = 5
    np_percentile = np.percentile(pd_pgs['proba'], q=list(range(0, 100 + int_step, int_step)))
    pd_pgs['percentile'] = np.searchsorted(np_percentile, pd_pgs['proba'], side='left') * int_step

    pd_prevalence_obs = pd_pgs.groupby('percentile').sum()[['target']]
    pd_prevalence_obs_count = pd_pgs.groupby('percentile').count()[['target']]
    pd_prevalence_obs = pd_prevalence_obs/pd_prevalence_obs_count
    pd_prevalence_obs.columns = ['obs']
    pd_prevalence_pre = pd_pgs.groupby('percentile').mean()[['proba']]
    pd_prevalence_pre.columns = ['pre']

    pd_rr_ingroup_case = pd_pgs.groupby('percentile').sum()[['target']]
    pd_rr_ingroup_sum = pd_pgs.groupby('percentile').count()[['target']]
    pd_rr = (pd_rr_ingroup_case / pd_rr_ingroup_sum) / (pd_pgs.sum()[['target']] / pd_pgs.count()[['target']])
    pd_rr.columns = ['Relative Risk']

    plt.figure(figsize=(5,5))
    sns.scatterplot(x=pd_prevalence_obs.index, y=pd_prevalence_obs['obs'], hue=pd_rr['Relative Risk'], palette=sns.cubehelix_palette(8, start=.5, rot=-.75, as_cmap=True))
    popt, pcov = sp.optimize.curve_fit(fsigmoid, pd_prevalence_pre.index, pd_prevalence_pre['pre'], method='dogbox', bounds=([0., 0.],[1., 100.]))
    sns.lineplot(x=pd_prevalence_pre.index, y=fsigmoid(pd_prevalence_pre.index, *popt), color="black")

    plt.legend(prop={'size': 12}, loc='upper left')
    plt.title(str_method + ' Predicting F1 Score: ' + "%.4f" % float_f1Score + ' ')
    plt.xlim(0, 100)
    plt.ylim(0, 1)
    plt.xlabel('Polygenic Score Percentile')
    plt.ylabel('Prevalence of Percentile Group')
    plt.savefig(os.path.join(str_outputFilePath, str("GenEpi_Prevalence_" + str_label + ".png")), dpi=300)
    plt.close('all')

    #-------------------------
    # plot ROC
    #-------------------------
    fpr, tpr, _ = skMetric.roc_curve(list_target, np.array(list_proba)[:,1])
    float_auc = skMetric.auc(fpr, tpr)

    plt.figure(figsize=(5,5))
    plt.plot(fpr, tpr, color='#e68fac', lw=2, label='Class 1 ROC curve (area = %0.2f)' % float_auc)
    plt.plot([0, 1], [0, 1], color='black', lw=2, linestyle='--')

    plt.legend(loc="lower right")
    plt.title('Receiver operating characteristic example')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.savefig(os.path.join(str_outputFilePath, str("GenEpi_ROC_" + str_label + ".png")), dpi=300)
    plt.close('all')

""""""""""""""""""""""""""""""
# main function
""""""""""""""""""""""""""""""
def main(args=None):
    ### obtain arguments from argument parser
    args = ArgumentsParser().parse_args(args)
    
    ### get arguments for I/O
    str_file_genotype = args.g
    str_path_model = args.m
    str_file_feature = args.f
    str_path_output = args.o

    if "Classifier" in str_path_model:
        list_predict, list_proba = IsolatedDataPredictor(str_file_genotype, str_path_model, str_file_feature, str_path_output, "c")
    else:
        list_predict, list_proba = IsolatedDataPredictor(str_file_genotype, str_path_model, str_file_feature, str_path_output, "r")

    ### plot prs
    if args.p is not None:
        str_file_phenotype = args.p
    
        list_target = []
        with open(str_file_phenotype, 'r') as file_inputFile:
            for line in file_inputFile:
                try: # modify (deal with NA)
                    list_target.append(float(line.strip().split(',')[-1]))
                except:
                    list_target.append(np.nan)

        # remove NA (modify)
        list_target = np.array(list_target)
        list_predict = np.array(list_predict)
        list_na = np.isnan(np.array(list_target))
        list_target = list_target[~list_na]
        list_predict = list_predict[~list_na]
        if list_proba: list_proba = np.array(list_proba)[~list_na]

        if "Classifier" in str_path_model:
            PlotPolygenicScore(list_target, list_predict, list_proba, str_path_output)
        else:
            print(sp.stats.pearsonr(list_target, list_predict))

if __name__ == "__main__":
    main()
