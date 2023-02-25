import argparse
import re
from datetime import datetime
import os
import pandas as pd

from lib import Utils, Misc
from lib import lib_snp as lib

'''
Example

python3 /yilun/prs-algo/docker/qc/Report/SNP_QC_report.py \
    -i /volume/prsdata/Users/yilun/RCVS/GT/QC/GT.QC \
    -C /yilun/prs-algo/docker/config.sh \
    -o /volume/prsdata/Users/yilun/Test/QC/RCVS.QC

python3 /yilun/prs-algo/docker/qc/Report/SNP_QC_report.py \
    -i /volume/prsdata/Users/chester/data/ADNI/output/gwas/qc/QualityControl/ADNI_merged_All_364_nomissid.QC \
    -C /yilun/prs-algo/docker/config.sh \
    -o /volume/prsdata/Users/yilun/Test/QC/ADNI_merged_All_364_Dedup.QC

python3 /yilun/prs-algo/docker/qc/Report/SNP_QC_report.py \
    -i /volume/prsdata/Users/yilun/Test/QC/QualityControl/TWB2_HEIGHT.base.QC \
    -C /yilun/TWB2/HEIGHT/config.sh \
    -o /volume/prsdata/Users/yilun/Test/QC/TWB2_HEIGHT.base.QC

'''

def parse_args():
    """
    Returns:
        arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input_prefix", type = str, required=True, help="input basename")
    parser.add_argument("-flipscan","--flipscan", type = str, default=None, help="inputname.flipscan")
    parser.add_argument("-smiss","--smiss", type = str, default=None, help="inputname.smiss")
    parser.add_argument("-vmiss","--vmiss", type = str, default=None, help="inputname.vmiss")
    parser.add_argument("-afreq","--afreq", type = str, default=None, help="inputname.afreq")
    parser.add_argument("-hwe","--hwe", type = str, default=None, help="inputname.hwe")
    parser.add_argument("-C","--config", type = str, required=True, help="Your gwas config")
    parser.add_argument("-p","--plot", type = int, default=1, help="Plot or not (1/0)")
    parser.add_argument("-o","--output_prefix", default="Report", help="Your gwas config")
    args = parser.parse_args()

    if args.afreq is None:
        args.afreq = f"{args.input_prefix}.afreq"

    if args.smiss is None:
        args.smiss = f"{args.input_prefix}.smiss"

    if args.vmiss is None:
        args.vmiss = f"{args.input_prefix}.vmiss"

    if args.hwe is None:
        args.hwe = f"{args.input_prefix}.hwe"

    if args.flipscan is None:
        args.flipscan = f"{args.input_prefix}.flipscan"

    return args

required_key = ["MAF_QC","GENO_QC","HWE_QC","HWE_QC_CONTROL","HWE_QC_CASE",]

FilterMAF = Utils.MakeFilter(lib.MAF_Filter)
FilterGENO = Utils.MakeFilter(lib.GENO_Filter)
FilterHWECase = Utils.MakeFilter(lib.HWE_Filter_CASE)
FilterHWEControl = Utils.MakeFilter(lib.HWE_Filter_CONTROL)
FilterHWEmissing = Utils.MakeFilter(lib.HWE_Filter_missing)
FilterHWEAll = Utils.MakeFilter(lib.HWE_Filter_ALL)
FilterFlip = Utils.MakeFilter(lib.FLIP_Filter)
FilterHWEmissing2 = Utils.MakeFilter(lib.HWE_Filter_missing2)

if __name__ == "__main__":
    args = parse_args()
    Misc.Logger(f"{args.output_prefix}.snp_qc.log","DEBUG")

    myConfig, info_list = Utils.ReadConfig(args.config, required_key = required_key)

    da_list = []
    threshod_dict_list = []
    da_list, threshod_dict_list, info_list = FilterMAF(args.afreq, myConfig, args.output_prefix, da_list, threshod_dict_list, info_list)    
    da_list, threshod_dict_list, info_list = FilterGENO(args.vmiss, myConfig, args.output_prefix, da_list, threshod_dict_list, info_list)
    da_list, threshod_dict_list, info_list = FilterFlip(args.flipscan, myConfig, args.output_prefix, da_list, threshod_dict_list, info_list)

    #HWE
    DA = pd.read_csv(args.hwe, sep = "\s+", low_memory=True)
    for fn in [FilterHWECase, FilterHWEControl, FilterHWEAll, ]:
        da_list, threshod_dict_list, info_list = fn(DA, myConfig, args.output_prefix, da_list, threshod_dict_list, info_list)
    del DA

    (DA, SNP_LIST), info_list = Utils.MergeAll(da_list, f"{args.output_prefix}.snp_qc", info_list)
    _, info_list = Misc.SaveITEM_LIST(SNP_LIST, output = f"{args.output_prefix}.snp_qc.snp_list", info_list = info_list)

    DD = {
        "Date": datetime.now().strftime('%y-%m-%d %H:%M:%S'),
        "Program": __file__,
        "Args": vars(args),
        "Filters": threshod_dict_list,
        "Config": myConfig,
        "Final_num": len(SNP_LIST), 
    }
    _, info_list = Misc.SaveJson(DD, output = f"{args.output_prefix}.snp_qc.json", info_list = info_list)

