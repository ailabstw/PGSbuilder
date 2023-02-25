#!/bin/bash

# trap error
set -eo pipefail
err_report(){
  echo "Fail at $1: $2"
}
trap 'err_report ${LINENO} "${BASH_COMMAND}"' ERR

# arguments
while getopts 'hi:m:b:g:ed:o:' flag; do
    case $flag in
        h)
            echo "GenEpi training"
            echo "options:"
            echo "-i, the input bfile prefix"
            echo "-m, the method ('c' or 'r', 'c' for classification and 'r' for regression), default=c"
            echo "-b, the resource directory of genome regions; dir=/volume/prsdata/GenEpi/data/"
            echo "-g, the resource of genome regions (hg19 or hg38), default=hg19"
            echo "-e, apply epistasis"
            echo "-d, the directory of working and output"
            echo "-o, the basename of output file"
            ;;
        i) BFILE=$OPTARG;;
        m) METHOD=$OPTARG;;
        b) DBDIR=$OPTARG;;
        g) GENOME=$OPTARG;;
        e) EPISTASIS="true";;
        d) WORK_DIR=$OPTARG;;
        o) BASENAME=$OPTARG;;
        *) echo "usage: $0 [-i] [-m] [-b] [-g] [-e] [-d] [-o]"; exit 1;;
    esac
done

if [ ! -f "$BFILE.bed" ] || [ ! -f "$BFILE.bim" ] || [ ! -f "$BFILE.fam" ]; then
    echo "-i missing or designating error"
    exit 1
fi

if [ -z "$METHOD" ]; then
    METHOD='c'
fi

if [ ! -d "$DBDIR" ]; then
    echo "-b no such directory"
    exit 1
fi

if [ -z "$GENOME" ]; then
    GENOME='hg19'
fi

if [ ! -d "$WORK_DIR" ]; then
    echo "-d no such directory"
    exit 1
fi

if [ -z "$BASENAME" ]; then
    echo "-o missing"
    exit 1
fi

# build oxford format
if [ ! -f "${WORK_DIR}/${BASENAME}.gen" ] || [ ! -f "${WORK_DIR}/${BASENAME}.sample" ]; then
    plink1.9 \
        --bfile "${BFILE}" \
        --chr 1-22 \
        --recode oxford \
        --out "${WORK_DIR}/${BASENAME}"
fi

# extract phenotype
if [ ! -f "${WORK_DIR}/${BASENAME}.phenotype.csv" ]; then
    awk 'NR>2{print $5}' "${WORK_DIR}/${BASENAME}.sample" > "${WORK_DIR}/${BASENAME}.phenotype.csv"
fi

# replace NA with control
sed -i 's/NA/0/' "${WORK_DIR}/${BASENAME}.phenotype.csv"

# GenEpi
## epistasis
if [ "$EPISTASIS" = "true" ]; then
    echo "Run GenEpi with epistasis"
    GenEpi \
        -g "${WORK_DIR}/${BASENAME}.gen" \
        -p "${WORK_DIR}/${BASENAME}.phenotype.csv" \
        -m "${METHOD}" \
        -s "${DBDIR}/${GENOME}_inter.txt" \
        -o "${WORK_DIR}/${BASENAME}" \
        -b "${GENOME}" \
        -k 5 \
        -t 16
else
    echo "Run GenEpi without epistasis"
    GenEpi \
        -g "${WORK_DIR}/${BASENAME}.gen" \
        -p "${WORK_DIR}/${BASENAME}.phenotype.csv" \
        -m "${METHOD}" \
        -s "${DBDIR}/${GENOME}_inter.txt" \
        -o "${WORK_DIR}/${BASENAME}" \
        -b "${GENOME}" \
        -k 5 \
        -t 16 \
        --noepistasis
fi

# remove unimportant files
EXTS=("gen" "sample" "phenotype.csv" "log" "nosex")
for EXT in "${EXTS[@]}";
do
    TMP_FILE=${WORK_DIR}/${BASENAME}.${EXT}
    rm ${TMP_FILE} || true
done
rm ${WORK_DIR}/${BASENAME}/*.vcf || true
rm -r ${WORK_DIR}/${BASENAME}/snpSubsets || true
mv ${WORK_DIR}/${BASENAME}/${BASENAME}.phenotype/crossGeneResult ${WORK_DIR}/${BASENAME}/
rm -r ${WORK_DIR}/${BASENAME}/${BASENAME}.phenotype || true