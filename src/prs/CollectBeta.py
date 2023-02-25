import json
import os
from unittest import main

from utils import *

def ArgumentsParser():
    ### define arguments
    import argparse
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-t", "--TARGET", type = str, required=True, help="TARGET")
    parser.add_argument("-s", "--SS_STR", type = str, required=True, help="SS_STR")
    parser.add_argument("-o", "--OUTDIR", type = str, required=True, help="OUTDIR")
    parser.add_argument("-a", "--ALGO", type = str, required=True, help="ALGO")
    args = parser.parse_args()

    args.ALGO = [ i.strip() for i in args.ALGO.split(",") ]
    assert len(args.ALGO) > 0
    assert os.path.exists(args.OUTDIR)
    assert os.path.exists(f"{args.OUTDIR}/LOG")
    assert os.path.exists(args.SS_STR)

    return args


def GetWeights(TARGET, SS_STR, OUTDIR):
    weights = Weights(TARGET, SS_STR, OUTDIR)
    weights()
    
    return weights
    
    
def CheckWeights(weights, algo_list):
    fail_list_error = []
    fail_list_nobeta = []
    for i in algo_list:
        if i not in weights.ss_df:
            fail_list_error.append(i)
            continue
        
        if weights.ss_df[i].isnull().all() or (weights.ss_df[i] == 0).all():
            fail_list_nobeta.append(i)
            weights.ss_df = weights.ss_df.drop(columns = i)

    return weights, fail_list_error, fail_list_nobeta


def SaveFile(OUTDIR, weights, algo_list, fail_list_error, fail_list_nobeta):
    # save beta
    weights.ss_df.to_csv(f"{OUTDIR}/beta.tsv", sep='\t', index=False)
    
    fail_list = [ i for i in algo_list if i in fail_list_error + fail_list_nobeta ]
    DD = {
        "TOOLS": algo_list, 
        "SUCCESS": [ i for i in algo_list if i not in fail_list ],
        "FAIL": fail_list, 
        "FAIL_LOG": {},
    }
    # Add log
    for algo in fail_list:
        log_file = f"{OUTDIR}/LOG/{algo}.log"
        if os.path.exists(log_file):
            with open(log_file, "r") as FF:
                log_lines = FF.readlines()
                
            if algo in fail_list_nobeta:
                log_lines = log_lines + [f"Beta estimated by {algo} are all 0.0 or NaN\n",
                                        f"There may be too few significant SNPs to compute polygenic risk score using {algo}\n"]
            DD["FAIL_LOG"][algo] = log_lines
        else:
            DD["FAIL_LOG"][algo] = []

    with open(f"{OUTDIR}/LOG/algo_status.json", "w") as FF:
        json.dump(DD, FF, indent = 4, separators = (',',': '))

    # fail_list for removal
    with open(f"{OUTDIR}/LOG/fail_list", "w") as FF:
        for i in fail_list:
            FF.write(f"{i}\n")


if __name__ == "__main__":
    args = ArgumentsParser()
    weights = GetWeights(args.TARGET, args.SS_STR, args.OUTDIR)
    weights, fail_list_error, fail_list_nobeta = CheckWeights(weights, args.ALGO)
    SaveFile(args.OUTDIR, weights, args.ALGO, fail_list_error, fail_list_nobeta)
    