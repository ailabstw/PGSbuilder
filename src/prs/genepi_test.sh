#!/bin/bash

# trap error
set -eo pipefail
err_report(){
  echo "Fail at $1: $2"
}
trap 'err_report ${LINENO} "${BASH_COMMAND}"' ERR

# arguments
while getopts 'hi:t:m:f:p:d:o:' flag; do
    case $flag in
        h)
            echo "GenEpi training"
            echo "options:"
            echo "-i, the input bfile prefix"
            echo "-t, the method ('c' or 'r', 'c' for classification and 'r' for regression), default=c"
            echo "-m, the model file"
            echo "-f, the feature file"
            echo "-p, the predictor file ([SRC_DIR]/GenEpi_predictor.py)"
            echo "-d, the directory of working and output"
            echo "-o, the basename of output file"
            ;;
        i) BFILE=$OPTARG;;
        t) METHOD=$OPTARG;;
        m) MODEL=$OPTARG;;
        f) FEATURE=$OPTARG;;
        p) PREDICTOR=$OPTARG;;
        d) WORK_DIR=$OPTARG;;
        o) BASENAME=$OPTARG;;
        *) echo "usage: $0 [-i] [-t] [-m] [-f] [-p] [-d] [-o]"; exit 1;;
    esac
done

if [ ! -f "$BFILE.bed" ] || [ ! -f "$BFILE.bim" ] || [ ! -f "$BFILE.fam" ]; then
    echo "-i missing or designating error"
    exit 1
fi

if [ -z "$METHOD" ]; then
    METHOD="c"
fi

if [ ! -f "$MODEL" ]; then
    echo "-m no such file"
    exit 1
fi

if [ ! -f "$FEATURE" ]; then
    echo "-f no such file"
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

### remove duplicates
plink2 \
    --bfile "${BFILE}" \
    --rm-dup force-first \
    --chr 1-22,X \
    --allow-extra-chr \
    --make-bed \
    --out "${WORK_DIR}/${BASENAME}.dedup"

# build oxford format
head "${FEATURE}" -n1 | sed -e 's/[,\*]/\n/g' | sed -E 's/_[ATCGatcgNn-]+\.[ATCGatcgNn-]+$//g'| sort | uniq > "${WORK_DIR}/${BASENAME}.snp_list"
if [ ! -f "${WORK_DIR}/${BASENAME}.gen" ] || [ ! -f "${WORK_DIR}/${BASENAME}.sample" ]; then
    plink1.9 \
        --bfile "${WORK_DIR}/${BASENAME}.dedup" \
        --allow-extra-chr \
        --extract "${WORK_DIR}/${BASENAME}.snp_list" \
        --recode oxford \
        --out "${WORK_DIR}/${BASENAME}"
fi

# extract phenotype
if [ ! -f "${WORK_DIR}/${BASENAME}.phenotype.csv" ]; then
    awk 'NR>2{print $5}' "${WORK_DIR}/${BASENAME}.sample" > "${WORK_DIR}/${BASENAME}.phenotype.csv"
fi

# GenEpi
python3 "${PREDICTOR}" \
    -g "${WORK_DIR}/${BASENAME}.gen" \
    -p "${WORK_DIR}/${BASENAME}.phenotype.csv" \
    -m "${MODEL}" \
    -f "${FEATURE}" \
    -o "${WORK_DIR}/${BASENAME}.pred"

if [ ${METHOD} = "c" ]; then
    echo "FID,IID,pred,score" > "${WORK_DIR}/${BASENAME}.pred.csv"
else
    echo "FID,IID,score" > "${WORK_DIR}/${BASENAME}.pred.csv"
fi
awk 'NR==FNR {a[NR]=$1","$2;next} FNR>1 {print a[FNR-1], $0}' OFS=',' "$BFILE.fam" "${WORK_DIR}/${BASENAME}.pred/Prediction.csv" >> "${WORK_DIR}/${BASENAME}.pred.csv"

# remove files
rm ${WORK_DIR}/${BASENAME}.dedup.* || true
EXTS=("gen" "sample" "phenotype.csv" "log" "nosex" "snp_list")
for EXT in "${EXTS[@]}";
do
    TMP_FILE=${WORK_DIR}/${BASENAME}.${EXT}
    rm ${TMP_FILE} || true
done