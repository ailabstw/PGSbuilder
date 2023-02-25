import argparse
import re
from datetime import datetime
import os
from lib import Utils, Misc
from lib import lib_id as lib
''' 
Example

python3 /yilun/prs-algo/gwas/Report/IND_QC_report.py \
    -i /volume/prsdata/GWAS/TWB2/train_qc/QualityControl/TWB2.train.QC \
    -f /volume/prsdata/GWAS/TWB2/train_qc/QualityControl/TWB2.train.QC.fam \
    -C /yilun/prs-algo/docker/config.sh \
    -o /volume/prsdata/Users/yilun/Test/QC/TWB2.train.QC

python3 /yilun/prs-algo/gwas/Report/IND_QC_report.py \
    -i /volume/prsdata/Users/chester/data/ADNI/output/gwas/qc/QualityControl/ADNI_merged_All_364_nomissid.QC \
    -f /volume/prsdata/Users/chester/data/ADNI/ADNI_merged_All_364_Dedup.fam \
    -C /yilun/prs-algo/docker/config.sh \
    -o /volume/prsdata/Users/yilun/Test/QC/ADNI_merged_All_364_Dedup.QC


'''


def parse_args():
    """
    Returns:
        arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input_prefix", type = str, required=True, help="input basename")
    parser.add_argument("-f","--fam", type = str, default=None, help="inputname.het")
    parser.add_argument("-het","--het", type = str, default=None, help="inputname.het")
    parser.add_argument("-smiss","--smiss", type = str, default=None, help="inputname.smiss")
    parser.add_argument("-sexcheck","--sexcheck", type = str, default=None, help="inputname.sexcheck")
    parser.add_argument("-king","--king", type = str, default=None, help="inputname.king.cutoff.in.id")
    parser.add_argument("-C","--config", type = str, required=True, help="Your gwas config")
    parser.add_argument("-o","--output_prefix", default="Report", help="Your gwas config")
    args = parser.parse_args()

    if args.het is None:
        args.het = f"{args.input_prefix}.het"

    if args.fam is None:
        args.fam = f"{args.input_prefix}.fam"

    if args.smiss is None:
        args.smiss = f"{args.input_prefix}.smiss"

    if args.sexcheck is None:
        args.sexcheck = f"{args.input_prefix}.sexcheck"

    if args.king is None:
        args.king = f"{args.input_prefix}.king.cutoff.in.id"

    return args

required_key = ["MAF_QC","GENO_QC","HWE_QC","HWE_QC_CONTROL","HWE_QC_CASE",]

FilterSEX = Utils.MakeFilter(lib.SEXCHECK_Filter)
FilterHET = Utils.MakeFilter(lib.HET_Filter)
FilterMIND = Utils.MakeFilter(lib.MIND_Filter)
ReaderKING = Utils.MakeFilter(lib.KING_Reader)


if __name__ == "__main__":
    args = parse_args()
    Misc.Logger(f"{args.output_prefix}.ind_qc.log","DEBUG")

    myConfig, info_list = Utils.ReadConfig(args.config, required_key = required_key)

    da_list = []
    threshod_dict_list = []
    da_list, threshod_dict_list, info_list = FilterSEX(args.sexcheck, myConfig, args.output_prefix, da_list, threshod_dict_list, info_list)    
    da_list, threshod_dict_list, info_list = FilterHET(args.het, myConfig, args.output_prefix, da_list, threshod_dict_list, info_list)    
    da_list, threshod_dict_list, info_list = FilterMIND(args.smiss, myConfig, args.output_prefix, da_list, threshod_dict_list, info_list)    
    da_list, threshod_dict_list, info_list = ReaderKING(args.king, myConfig, args.output_prefix, da_list, threshod_dict_list, info_list)    

    (DA, IND_LIST), info_list = Utils.MergeAll(da_list, f"{args.output_prefix}.ind_qc", info_list)
    IND, info_list = Utils.SaveIndList(args.fam, IND_LIST, output = f"{args.output_prefix}.ind_qc.ind_list", info_list = info_list)

    DD = {
        "Date": datetime.now().strftime('%y-%m-%d %H:%M:%S'),
        "Program": __file__,
        "Args": vars(args),
        "Filters": threshod_dict_list,
        "Config": myConfig,
        "Final_num": len(IND_LIST), 
    }
    _, info_list = Misc.SaveJson(DD, output = f"{args.output_prefix}.ind_qc.json", info_list = info_list)

    p, info_list = Utils.PlotHetMiss(DA, threshod_dict_list, f"{args.output_prefix}.hetxmiss.png", info_list)
