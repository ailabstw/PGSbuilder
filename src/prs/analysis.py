#!/usr/bin/python3

import os, sys, argparse
from utils import *


def ArgumentParser():
    parser = argparse.ArgumentParser(prog='PRS Analysis', description='cohort reference, covariates, and performance')
    parser.add_argument('--pred_file', required=True, help='the prediction file (.csv)')
    parser.add_argument('--method', required=True, help='clf or reg')
    parser.add_argument('--mode', required=True, help='target or test')
    parser.add_argument('--out_dir', required=True, help='the output directory')
    parser.add_argument('--rank_ref_file', required=False, default="", help='the rank reference file (.csv)')
    parser.add_argument('--cov', required=False, default="", help='the covariate file (.tsv)')
    parser.add_argument('--cov_ref_dir', required=False, default="", help='the directory of covariate regression model, including scaler.joblib, reg.joblib')
    parser.add_argument('--run_performance', action='store_true', help='whether to calculate the model performance')
    parser.add_argument('--percentile_num', required=False, default=10, type=int, help='the number of percentile groups, default=10')
    return parser


def main(args=None):
    ### arguments
    args = ArgumentParser().parse_args(args)

    # prediction dataframe
    pred_df = pd.read_csv(args.pred_file)

    # output dir
    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)

    # check mode
    if args.mode == 'test':
        if not os.path.isfile(args.rank_ref_file):
            print('Rank reference file is required')
            sys.exit()
        if (os.path.isfile(args.cov)) & (not os.path.isdir(args.cov_ref_dir)):
            print('If the covariate file is provied, covariate models are required in test mode')
            sys.exit()


    ### cohort reference
    if args.mode == 'target':
        cohort_ref = CohortRef(pred_df)
        rank_df = cohort_ref()
        cohort_ref.rank_ref_df.to_csv(f'{args.out_dir}/rank_ref.csv')
        cohort_ref.hist_ref_df.to_csv(f'{args.out_dir}/hist_ref.csv')
    else:
        cohort_ref = CohortRef(pred_df, args.rank_ref_file)
        rank_df = cohort_ref()
    rank_df.to_csv(f'{args.out_dir}/rank.csv', index=False)


    ### performance
    if args.run_performance:
        analysis = Analysis(pred_df, args.method, args.out_dir)
        analysis(percentile_num=args.percentile_num)

    ### covariates
    if os.path.isfile(args.cov):
        try:
            cov_df = pd.read_csv(args.cov, sep='\t')
            cov = CovResults(pred_df, cov_df, args.method)
        except:
            cov_df = pd.read_csv(args.cov, sep=',')
            cov = CovResults(pred_df, cov_df, args.method)

        # outdir
        if not os.path.isdir(f'{args.out_dir}/cov'):
            os.mkdir(f'{args.out_dir}/cov')
        
        # target
        if args.mode == 'target':
            # train regression model
            cov_pred_df, cov_model_dict = cov.Train()
            json.dump(cov_model_dict, open(f'{args.out_dir}/cov/models.json', 'w'))

            # save weight
            weight_dict = dict()
            for k, v in cov_model_dict.items():
                weight_dict[k] = dict(zip(v['columns'], v['reg_coef']))
            json.dump(weight_dict, open(f'{args.out_dir}/cov/weight.json', 'w'))
            
            # build cohort reference
            cov_cohort_ref = CohortRef(cov_pred_df)
            cov_rank_df = cov_cohort_ref()
            cov_cohort_ref.rank_ref_df.to_csv(f'{args.out_dir}/cov/rank_ref.csv')
            cov_cohort_ref.hist_ref_df.to_csv(f'{args.out_dir}/cov/hist_ref.csv')
        
        # test
        else:
            # test regression model
            cov_pred_df = cov.Test(f'{args.cov_ref_dir}/models.json')
            
            # calculate rank
            cov_cohort_ref = CohortRef(cov_pred_df, f'{args.cov_ref_dir}/rank_ref.csv')
            cov_rank_df = cov_cohort_ref()

        # save files
        cov_pred_df.to_csv(f'{args.out_dir}/cov/prediction.csv', index=False)
        cov_rank_df.to_csv(f'{args.out_dir}/cov/rank.csv', index=False)
        
        # cov performance
        if args.run_performance:
            cov_analysis = Analysis(cov_pred_df, args.method, f'{args.out_dir}/cov')
            cov_analysis(percentile_num=args.percentile_num)


if __name__ == '__main__':
    main()
