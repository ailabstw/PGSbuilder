#!/bin/bash

# trap error
set -eo pipefail
err_report(){
  echo "Fail at $1: $2"
}
trap 'err_report ${LINENO} "${BASH_COMMAND}"' ERR

# arguments
while getopts 'hi:c:a:m:p:d:o:' flag; do
    case $flag in
        h)
            echo "Clumping and thresholding for PRS - training (absolute paths are required)"
            echo "options:"
            echo "-i, the input bfile prefix"
            echo "-c, the covariates file (optional)"
            echo "-a, the filename of summary statistics"
            echo "-m, method (classification(clf) or regression(reg), default='clf'"
            echo "-p, the directory of src (clump_threshold_best_fit.R)"
            echo "-d, the directory of working and output"
            echo "-o, the basename of output file"
            ;;
        i) BFILE=$OPTARG;;
        c) COVAR_FILE=$OPTARG;;
        a) ASSOC_FILE=$OPTARG;;
        m) METHOD=$OPTARG;;
        p) OPT_DIR=$OPTARG;;
        d) WORK_DIR=$OPTARG;;
        o) BASENAME=$OPTARG;;
        *) echo "usage: $0 [-i] [-c] [-a] [-m] [-p] [-d] [-o]"; exit 1;;
    esac
done

if [ ! -f "$BFILE.bed" ] || [ ! -f "$BFILE.bim" ] || [ ! -f "$BFILE.fam" ]; then
    echo "-i missing or designating error"
    exit 1
fi

if [ ! -f "$ASSOC_FILE" ]; then
    echo "-a missing or designating error"
    exit 1
fi

if [ -z "$METHOD" ]; then
    METHOD='clf'
fi

if [ ! -f "$OPT_DIR/clump_threshold_best_fit.R" ]; then
    echo "-p missing"
    exit 1
fi

if [ ! -d "$WORK_DIR" ]; then
    echo "-d no such directory"
    exit 1
fi

if [ -z "$BASENAME" ]; then
    echo "-o missing"
    exit 1
fi

# clumping
plink1.9 \
    --bfile "$BFILE" \
    --clump-p1 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump "$ASSOC_FILE" \
    --clump-snp-field "ID" \
    --clump-field "P" \
    --allow-extra-chr \
    --allow-no-sex \
    --out "${WORK_DIR}"/"${BASENAME}"

# extract valid SNPs
awk 'NR!=1{print $3}' "${WORK_DIR}"/"${BASENAME}".clumped > "${WORK_DIR}"/"${BASENAME}".valid.snp

# extract P-value
P_COLNUM=$(sed -n '1s/\s/\n/gp' $ASSOC_FILE | grep -nx 'P' | cut -d: -f1)
awk -v p_colnum=$P_COLNUM '{print $3,$p_colnum}' "${ASSOC_FILE}" > "${WORK_DIR}"/"${BASENAME}".pvalue

# list of P-value 
{
    echo "1e-8 0 1e-8"
    echo "1e-7 0 1e-7"
    echo "1e-6 0 1e-6"
    echo "1e-5 0 1e-5"
    echo "1e-4 0 1e-4"
    echo "1e-3 0 1e-3"
    echo "1e-2 0 1e-2"
    echo "1e-1 0 1e-1"
    echo "1 0 1"
} > "$WORK_DIR"/pvalue_ranges

# PRS
BETA_COLNUM=$(sed -n '1s/\s/\n/gp' $ASSOC_FILE | grep -nx 'BETA' | cut -d: -f1)
plink1.9 \
    --bfile "$BFILE" \
    --score "$ASSOC_FILE" 3 6 "${BETA_COLNUM}" header sum \
    --q-score-range "$WORK_DIR"/pvalue_ranges "${WORK_DIR}"/"${BASENAME}".pvalue \
    --extract "${WORK_DIR}"/"${BASENAME}".valid.snp \
    --allow-extra-chr \
    --allow-no-sex \
    --out "${WORK_DIR}"/"${BASENAME}"

# select best fit threshold
if [ -f "$COVAR_FILE" ]; then
    Rscript "$OPT_DIR"/clump_threshold_best_fit.R \
        -i "$BFILE" \
        -c "$COVAR_FILE" \
        -p "$WORK_DIR/$BASENAME" \
        -m "$METHOD" \
        -d "$WORK_DIR" \
        -o "$BASENAME"
else
    Rscript "$OPT_DIR"/clump_threshold_best_fit.R \
        -i "$BFILE" \
        -p "$WORK_DIR/$BASENAME" \
        -m "$METHOD" \
        -d "$WORK_DIR" \
        -o "$BASENAME"
fi