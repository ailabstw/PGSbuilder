#!/bin/bash

# trap error
set -eo pipefail
err_report(){
  echo "Fail at $1: $2"
}
trap 'err_report ${LINENO} "${BASH_COMMAND}"' ERR


### arguments
while getopts 'hi:b:m:d:o:rg:' flag; do
    case $flag in
        h)
            echo "Predict PRS using adjusted beta"
            echo "options:"
            echo "-i, the input bfile prefix"
            echo "-b, the beta file (beta.csv)"
            echo "-m, the method of the PRS (reg or clf)"
            echo "-d, the directory of working and output"
            echo "-o, the basename of output file"
            echo "-r, dedup"
            echo "-g, the GenEpi prediction file"
            ;;
        i) BFILE=$OPTARG;;
        b) MODEL=$OPTARG;;
        m) METHOD=$OPTARG;;
        d) WORK_DIR=$OPTARG;;
        o) BASENAME=$OPTARG;;
        r) DEDUP="true";;
        g) GENEPI_PRED=$OPTARG;;
        *) echo "usage: $0 [-i] [-b] [-m] [-d] [-o] [-r] [-g]"; exit 1;;
    esac
done

REAL_PATH=$(realpath $0)
SRC_DIR=$(dirname ${REAL_PATH})


### dedup
if [ "${DEDUP}" = "true" ];then
    printf "\n\n\n###### Dedup the Bfile ######\n\n\n"
    plink2 \
        --bfile "${BFILE}" \
        --rm-dup force-first \
        --allow-extra-chr \
        --make-bed \
        --out "${WORK_DIR}/${BASENAME}.dedup"
        BFILE="${WORK_DIR}/${BASENAME}.dedup"
fi


### predict by algorithms
# MODEL columns = ['#CHROM','POS','ID','REF','ALT','A1','P','LOG10_P','BETA', ALGO_1, ALGO_2, ...]
read -ra COLS < ${MODEL}
for ((i=9; i<${#COLS[@]}; i++));do
    ALGO=${COLS[i]}
    printf "\n\n\n###### Predict with ${ALGO} ######\n\n\n"
    plink1.9 \
        --bfile "${BFILE}" \
        --score "${MODEL}" 3 6 $((i+1)) header sum \
        --score-no-mean-imputation \
        --allow-extra-chr \
        --allow-no-sex \
        --out "${WORK_DIR}/${BASENAME}.${ALGO}"
done
[[ -f ${GENEPI_PRED} ]] && mv ${GENEPI_PRED} "${WORK_DIR}/${BASENAME}.GenEpi.csv"


### merge predictions
cd ${SRC_DIR} || exit
python3 - << EOF
from utils import *
pred = PRSResults("${BFILE}", "${WORK_DIR}/${BASENAME}", "${METHOD}")
df = pred()
df.to_csv("${WORK_DIR}/prediction.csv", index=False)
EOF
cd ${OUTDIR} || exit


### remove temp files
rm ${WORK_DIR}/${BASENAME}.dedup.* || true