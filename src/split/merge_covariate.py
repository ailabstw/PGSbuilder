# categorical features should be 1 as control and 2 as case; NA: 0, -9, na, nan

import os, sys, argparse
import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description='Merge PCA and covariate files.')
    parser.add_argument('--fam', required=True, help='the fam file')
    parser.add_argument('--pca', required=False, default='', help='the eigenvector file')
    parser.add_argument('--cov', required=False, default='', help='the covariate file')
    parser.add_argument('--out', required=True, help='the output file')
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    
    # read fam file
    df = pd.read_csv(args.fam, sep='\s+', names=['FID', 'IID', 'father', 'mother', 'sex', 'phenotype'])
    df = df[['FID', 'IID']]
    df.FID = df.FID.astype(str)
    df.IID = df.IID.astype(str)

    # read PCA file
    if args.pca != '':
        temp_df = pd.read_csv(args.pca, sep='\s+')
        temp_df = temp_df.rename(columns={'#FID': 'FID'})
        temp_df.drop_duplicates(subset=['FID', 'IID'], inplace=True)
        temp_df.FID = temp_df.FID.astype(str)
        temp_df.IID = temp_df.IID.astype(str)
        df = df.merge(temp_df, on=['FID', 'IID'], how='left')

    # read covariate file
    if args.cov != '':
        temp_df = pd.read_csv(args.cov)
        if len(temp_df.columns) < 3:
            temp_df = pd.read_csv(args.cov, sep = "\s+")
        temp_df.FID = temp_df.FID.astype(str)
        temp_df.IID = temp_df.IID.astype(str)

        temp_df.drop_duplicates(subset=['FID', 'IID'], inplace=True)
        df = df.merge(temp_df, on=['FID', 'IID'], how='left')

    # replace -9 with NA
    df = df.replace(-9, np.nan)

    # save
    df.to_csv(args.out, sep='\t', index=False, na_rep='NaN')


if __name__=='__main__':
    main()
