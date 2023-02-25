#!/usr/bin/python3
import os, sys
import numpy as np
import pandas as pd


def parse_argumnet():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input_glm", type = str,  required = True, help="cov csv file")
    parser.add_argument("-o","--out_path", type = str,  default = None, help="cov csv file")
    parser.add_argument("--chr", type = str, default="chr", help="plink fam file")
    parser.add_argument("--pos", type=str, default="position")
    parser.add_argument("--id", type = str, default = "SNPs", help="output_path")
    parser.add_argument("--a1", type = str, default = "a1", help="output_path")
    parser.add_argument("--a2", type = str, default = "a2", help="output_path")
    parser.add_argument("--obs", type = str, default = "SIZE", help="output_path")
    parser.add_argument("--se", type = str, default = "PARA_SE", help="output_path")
    parser.add_argument("--stat", type = str, default = "TZ", help="output_path")
    parser.add_argument("--freq", type = str, default = "A1-FREQ", help="output_path")
    parser.add_argument("--pvalue", type = str, default ="LOG10-P", help="output_path")
    parser.add_argument("--logp", type = str, default = "P-value", help="output_path")
    args = parser.parse_args()
    
    return args


def ColNorm(df, args):
    synonym_dict = {
        '#CHROM': [ 'CHR', 'CHROM', 'CHROMOSOME', args.chr ],
        'POS': [ 'POSITION', 'BP', args.pos ],
        'ID': [ 'SNP', 'RSID', args.id ],
        'A1': [ 'ALLELE1', 'ALLELE_1', args.a1 ],
        'A2': [ 'ALLELE2', 'ALLELE_2', 'A0', 'ALLELE0', 'ALLELE_0', 'AX', args.a2 ],
        'A1_FREQ': [ 'AF', 'A1FREQ', 'FREQUENCY', 'ALLELE_FREQ', 'ALLELE_FREQUENCY', args.freq ],
        'OBS_CT': [ 'COUNT', 'N_EFF', 'NEFF', 'NUMBER', args.obs ],
        'BETA_SE': [ 'SE', 'LOG(OR)_SE', args.se ],
        'STAT': [ 'T_STAT', 'Z_STAT', "STAT", args.stat ],
        'LOG10_P': [ 'LOG_P', 'LOGP', 'LOG10P', args.logp ],
        'P': [ 'pvalue', 'Pvalue', args.pvalue ],
    }
    for name, synonym in synonym_dict.items():
        df = df.rename(columns={i:name for i in synonym})
    return df


# #CHROM, POS, ID, REF, ALT, A1, OBS_CT, BETA / OR, BETA_SE, P / LOG10_P
def CheckRequired(cols):
    check = True
    for col in ['#CHROM', 'POS', 'ID', 'REF', 'ALT']:
        if col not in cols:
            print(f'Columns related to SNP information is required, but "{col}" is missed.')
            check = False
    if 'A1' not in cols:
        print('Column A1 is required for recognizing the effective allele.')
        check = False
    if 'OBS_CT' not in cols:
        print('Column related to the sample size is required, e.g. OBS_CT, NUM, N_EFF.')
        check = False
    if ('OR' not in cols) & ('BETA' not in cols):
        print('Column of odds ratio or beta is required to represent the effect size.')
        check = False
    if 'BETA_SE' not in cols:
        print('Column of the standard error of beta (BETA_SE) is required.')
        check = False
    if ('P' not in cols) & ('LOG10_P' not in cols):
        print('Column of P-value or -log10(P) is required.')
        check = False
    return check


def Main(args):
    filename = args.input_glm
    if filename.endswith('.csv'):
        df = pd.read_csv(filename, header=0)
    elif filename.endswith('.tsv'):
        df = pd.read_csv(filename, header=0, sep='\t')
    else:
        df = pd.read_csv(filename, header=0, sep='\s+')
    df.columns = df.columns.str.upper()
    # normalize column name
    df = ColNorm(df, args)
    cols = df.columns.tolist()
    # check required columns
    if not CheckRequired(cols):
        print('Some required columns are missed.')
        sys.exit(1)
    else:
        print('Pass the check of required columns.')

    # add OR or BETA
    try:
        if ('OR' in cols) & ('BETA' not in cols):
            df['BETA'] = np.log(df['OR'])
        elif ('OR' not in cols) & ('BETA' in cols):
            df['OR'] = np.exp(df['BETA'])
    except:
        print('Fail to add OR or BETA column')
        sys.exit(1)

    # add P or LOG10_P
    try:
        if ('P' in cols) & ('LOG10_P' not in cols):
            df['LOG10_P'] = -np.log10(df['P'])
        elif ('P' not in cols) & ('LOG10_P' in cols):
            df['P'] = 10**(-df['LOG10_P'])
    except:
        print('Fail to add P-value or -log10(P) column')
        sys.exit(1)

    # add A2
    try:
        if 'A2' not in cols:
            df['A2'] = df['REF']
            df.loc[df['A1']==df['REF'], 'A2'] = df.loc[df['A1']==df['REF'], 'ALT']
    except:
        print('Fail to add A2 column')
        sys.exit(1)

    # add other columns
    for col in ['A1_FREQ', 'STAT']:
        if col not in cols:
            df[col] = np.nan

    # ordered columns
    df = df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'A1', 'A2', 'A1_FREQ', 'OBS_CT', 'OR', 'BETA', 'BETA_SE', 'STAT', 'P', 'LOG10_P']]
    if args.out_path is None:
        args.out_path = filename
    df.to_csv(args.out_path, sep='\t', index=False, na_rep='NaN')


if __name__ == '__main__':
    args = parse_argumnet()
    Main(args)




'''
python3 /yilun/prs-algo/assoc/modify_sumstats.py \
    -i /mnt/prsdata/Test/OUTPUT/DEMO_REG/gwas/yilun-reg-cov_gwas/demo_hg38.QC.PHENO1.glm.linear \
    -o /tmp/TWB2_LDL_extSumStats

'''
