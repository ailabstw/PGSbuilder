#!/bin/bash

# trap error
set -eo pipefail
err_report(){
  echo "Fail at $1: $2"
}
trap 'err_report ${LINENO} "${BASH_COMMAND}"' ERR

# arguments
helpFunction()
{
    echo ""
    echo "Pipeline for PRSice training and testing"
    echo "options:"
    echo "    -r, the basename of training bfile"
    echo "    -a, the filename of summary statistics"
    echo "    -m, the method (reg or clf), default=clf"
    echo "    -d, the directory of working and output"
    echo "    -o, the basename of output bfile"
    echo ""
    echo "example:"
    echo "sh PRSice.sh -r /volume/prsdata/gwas_test/TWB1_T2D/train_qc/TWB1_T2D.train.QC.PS -a /volume/prsdata/gwas_test/TWB1_T2D/train_assoc/TWB1_T2D.train.PHENO1.glm.logistic -d /volume/prsdata/Users/chester/prsice2/output/TWB1_T2D -o TWB1_T2D"
    exit 1 # Exit script after printing help
}

while getopts 'r:a:m:d:o:' flag; do
    case $flag in
        r ) TRAIN_BFILE=$OPTARG;;
        a ) ASSOC_FILE=$OPTARG;;
        m ) METHOD=$OPTARG;;
        d ) WORK_DIR=$OPTARG;;
        o ) BASENAME=$OPTARG;;
        ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
    esac
done

if [ ! -f "$TRAIN_BFILE.bed" ] || [ ! -f "$TRAIN_BFILE.bim" ] || [ ! -f "$TRAIN_BFILE.fam" ]; then
    echo "-tr missing or designating error"
    helpFunction
fi

if [ ! -f "$ASSOC_FILE" ]; then
    echo "-a missing or designating error"
    helpFunction
fi

if [ -z "$METHOD" ]; then
    METHOD='clf'
fi

if [ -z "$WORK_DIR" ]; then
    echo "-d missing"
    helpFunction
fi

if [ -z "$BASENAME" ]; then
    echo "-o missing"
    helpFunction
fi

[ $METHOD == 'clf' ] && BINARY_TARGET=T || BINARY_TARGET=F
PRSice.R \
    --prsice /usr/local/bin/PRSice_linux \
    --base "$ASSOC_FILE" \
    --snp ID --chr CHROM --bp POS --A1 A1 --A2 A2 --pvalue P \
    --target "$TRAIN_BFILE" \
    --binary-target ${BINARY_TARGET} --base-maf MAF:0.01 --base-info INFO:0.8 --stat BETA --beta \
    --perm 100 --print-snp \
    --out "${WORK_DIR}/${BASENAME}"

PVALUE=$(awk 'NR==2 {print $3}' "${WORK_DIR}/${BASENAME}.summary")
echo "${PVALUE} 0 ${PVALUE}" > "${WORK_DIR}/best_pvalue_range"

P_COLNUM=$(sed -n '1s/\s/\n/gp' $ASSOC_FILE | grep -nx 'P' | cut -d: -f1)
awk -F" " -v p_colnum=$P_COLNUM '{ print $3,$p_colnum }' "${ASSOC_FILE}" > "${WORK_DIR}/${BASENAME}.pvalue"

awk -F" " '{ print $2 }' "${WORK_DIR}/${BASENAME}.snp" > "${WORK_DIR}/${BASENAME}.valid.snp"