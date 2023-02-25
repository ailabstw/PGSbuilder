#!/bin/bash

# trap error
set -eo pipefail
err_report(){
  echo "Fail at $1: $2"
}
trap 'err_report ${LINENO} "${BASH_COMMAND}"' ERR

# arguments
while getopts 'hi:c:d:o:C:' flag; do
    case $flag in
        h)
            echo "Pipeline for Summary statistics (absolute paths are required)"
            echo "options:"
            echo "-i, the basename of input bfile"
            echo "-c, the filename of covariate"
            echo "-d, the directory of working and output"
            echo "-o, the basename of output bfile"
            echo "-C, the path to config file"
            ;;
        i) BFILE=$OPTARG;;
        c) COVFILE=$OPTARG;;
        d) WORK_DIR=$OPTARG;;
        o) BASENAME=$OPTARG;;
        C) CONFIG=$OPTARG;;
        *) echo "usage: $0 [-i] [-c] [-d] [-o] [-C]"; exit 1;;
    esac
done

if [ ! -f "$BFILE.bed" ] || [ ! -f "$BFILE.bim" ] || [ ! -f "$BFILE.fam" ]; then
    echo "-i missing or designating error"
    exit 1
fi

if [ -z "$WORK_DIR" ];  then
    echo "-d missing, set as $(realpath .)"
    WORK_DIR=$(realpath .)
fi

if [ -z "$BASENAME" ]; then
    echo "-o missing"
    exit 1
fi

if [ -z "$CONFIG" ]; then
    echo "-c missing"
fi

source "${CONFIG}"
REAL_PATH=$(realpath $0)
SRC_DIR=$(dirname ${REAL_PATH})
mkdir -p "$WORK_DIR"
cd "$WORK_DIR" || exit

# pca 
if [ "${PCA_FLAG}" = "true" ];then
    # pruning
    {
    plink2 \
        --bfile "${BFILE}" \
        --indep-pairwise "$PRUNE_WINDOW" "$PRUNE_STEP" "$PRUNE_THRESHOLD" \
        --memory "${MEMORY}" \
        --threads "${THREAD}" \
        --out "${WORK_DIR}/${BASENAME}"
    } 2>&1 | tee "${WORK_DIR}/${BASENAME}.prune.log"

    # pca
    SAMPLE_NUM=$(wc -l ${BFILE}.fam | cut -d" " -f1)
    [ ${SAMPLE_NUM} -gt 5000 ] && PCA_CMD="--pca ${PCA_COUNT} approx" || PCA_CMD="--pca ${PCA_COUNT}"
    {
    plink2 \
        --bfile "${BFILE}" \
        --extract "${WORK_DIR}/${BASENAME}.prune.in" \
        ${PCA_CMD} \
        --memory "${MEMORY}" \
        --threads "${THREAD}" \
        --out "${WORK_DIR}/${BASENAME}" 
    } 2>&1 | tee "${WORK_DIR}/${BASENAME}.pca.log"
    
    PCAFILE="${WORK_DIR}/${BASENAME}.eigenvec"
else
    PCAFILE=""
fi

# covariate
[ ! -f "${COVFILE}" ] && COVFILE=""
python3 "${SRC_DIR}/../split/merge_covariate.py" \
    --fam "${BFILE}.fam" \
    --pca "${PCAFILE}" \
    --cov "${COVFILE}" \
    --out "${WORK_DIR}/${BASENAME}.cov"

# association
{
plink2 \
    --bfile "${BFILE}" \
    --chr 1-22 \
    --glm hide-covar log10 cols=chrom,pos,ref,alt,ax,a1freq,nobs,orbeta,se,ci,tz,p \
    --pfilter 1 \
    --covar "${WORK_DIR}/${BASENAME}.cov" \
    --covar-variance-standardize \
    --out "${WORK_DIR}/${BASENAME}"
} 2>&1 | tee "${WORK_DIR}/${BASENAME}.assoc.log"

# modify summary statistics
SS_FILE=${WORK_DIR}/${BASENAME}.PHENO1.glm.*
python3 "${SRC_DIR}/modify_sumstats.py" -i ${SS_FILE}

# extract significant SNPs
P_COLNUM=$(sed -n '1s/\s/\n/gp' $SS_FILE | grep -nx 'P' | cut -d: -f1)
awk 'NR==1{print $0}' ${SS_FILE} > "${WORK_DIR}/${BASENAME}.significant.txt"
awk -v threshold=$THRESHOLD -v p_colnum=$P_COLNUM '{if ($p_colnum < threshold) print $0}' ${SS_FILE} >> "${WORK_DIR}/${BASENAME}.significant.txt"

# Manhattan plot
THRESHOLD_LOGP=$(awk -v threshold=$THRESHOLD 'BEGIN {print -log(threshold)/log(10)}')
OVER_SIG_COUNT=$(awk -v threshold=1e-100 -v p_colnum=$P_COLNUM 'BEGIN {count=0} {if ($p_colnum < threshold) count+=1} END {print count}' ${SS_FILE})
[ ${OVER_SIG_COUNT} != 0 ] && MAX_Y_LIM_CMD="--max-ylim 100" || MAX_Y_LIM_CMD=""

# Get adjusted
#$PLINK2 --adjust-file ${SS_FILE} \
#    --out "${WORK_DIR}/${BASENAME}"
#grep "lambda" "${WORK_DIR}/${BASENAME}".log | grep -oP "[0-9]\.[0-9]+" > ${WORK_DIR}/${BASENAME}.lambdaGC

manhattan_generator \
    --twopoint ${SS_FILE} \
    --col-chr "#CHROM" \
    --col-name "ID" \
    --col-pos "POS" \
    --col-pvalue "P" \
    --bp \
    --use-pvalues \
    --abline "${THRESHOLD_LOGP}" \
    --significant-threshold "${THRESHOLD_LOGP}" $MAX_Y_LIM_CMD \
    --no-annotation \
    --significant-point-size 2 \
    --point-size 1 \
    --graph-title "${BASENAME}" \
    --chr-text-size 10 \
    --exclude-chr 23,24 \
    --output "${WORK_DIR}/${BASENAME}.manhattan"
