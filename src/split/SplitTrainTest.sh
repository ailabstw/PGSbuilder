#!/bin/bash

# trap error
set -eo pipefail
err_report(){
  echo "Fail at $1: $2"
}
trap 'err_report ${LINENO} "${BASH_COMMAND}"' ERR

# arguments
while getopts 'hs:a:b:c:v:m:p:o:i:j:' flag; do
    case $flag in
        h)
            echo "Split training and testing"
            echo "options:"
            echo "-s, the directory of source codes"
            echo "-a, input bed file"
            echo "-b, input bim file"
            echo "-c, input fam file"
            echo "-m, method (clf or reg)"
            echo "-d, drop NA"
            echo "-r, random split"
            echo "-p, the test ratio (default=0.1)"
            echo "-o, output basename"
            echo "-i, the suffix of training"
            echo "-j, the suffix of testing"
            ;;
        s) SRC_DIR=$OPTARG;;
        a) BED=$OPTARG;;
        b) BIM=$OPTARG;;
        c) FAM=$OPTARG;;
        v) COVFILE=$OPTARG;;
        m) METHOD=$OPTARG;;
        p) TEST_RATIO=$OPTARG;;
        o) OUT_BASENAME=$OPTARG;;
        i) TRAIN_SUFFIX=$OPTARG;;
        j) TEST_SUFFIX=$OPTARG;;
        *) echo "usage: $0 [-s] [-a] [-b] [-c] [-v] [-m] [-d] [-r] [-p] [-o] [-i] [-j]"; exit 1;;
    esac
done

if [ ! -f "$BED" ] || [ ! -f "$BIM" ] || [ ! -f "$FAM" ]; then
    echo "bfile missing or designating error"
    exit 1
fi

if [ -z "$METHOD" ]; then
    echo "method (-m) is needed"
    exit 1
fi

if [ -z "$TEST_RATIO" ]; then
    TEST_RATIO=0.1
fi

echo $TRAIN_SUFFIX
echo $OUT_BASENAME
OUT_DIR=$(dirname $OUT_BASENAME)
if [ ! -d "$OUT_DIR" ]; then
    echo "output directory dose not exists"
    exit 1
fi

if [ -z "$TRAIN_SUFFIX" ]; then
    TRAIN_SUFFIX="train"
fi

if [ -z "$TEST_SUFFIX" ]; then
    TEST_SUFFIX="test"
fi

# split list
DROP_CMD="true"
RANDOM_CMD="true"
[ "$DROP_NA" = "true" ] && DROP_CMD="--dropna" || DROP_CMD=""
[ "$RANDOM_SPLIT" = "true" ] && RANDOM_CMD="--random" || RANDOM_CMD=""
python3 "${SRC_DIR}/SplitTrainTest.py" \
    --fam_file "${FAM}" \
    --method "${METHOD}" \
    --testing_ratio "${TEST_RATIO}" $DROP_CMD $RANDOM_CMD \
    --out_basename "${OUT_BASENAME}" \
    --train_suffix "${TRAIN_SUFFIX}" \
    --test_suffix "${TEST_SUFFIX}"

# train
plink1.9 \
        --bed "${BED}" \
        --bim "${BIM}" \
        --fam "${FAM}" \
        --keep "${OUT_BASENAME}.${TRAIN_SUFFIX}.list" \
        --allow-no-sex \
        --allow-extra-chr \
        --keep-allele-order \
        --make-bed \
        --out "${OUT_BASENAME}.${TRAIN_SUFFIX}"

if [ -f "${COVFILE}" ];then
    python3 "${SRC_DIR}/merge_covariate.py" \
        --fam "${OUT_BASENAME}.${TRAIN_SUFFIX}.fam" \
        --cov "${COVFILE}" \
        --out "${OUT_BASENAME}.${TRAIN_SUFFIX}.cov"
fi

# test
plink1.9 \
        --bed "${BED}" \
        --bim "${BIM}" \
        --fam "${FAM}" \
        --keep "${OUT_BASENAME}.${TEST_SUFFIX}.list" \
        --allow-no-sex \
        --allow-extra-chr \
        --keep-allele-order \
        --make-bed \
        --out "${OUT_BASENAME}.${TEST_SUFFIX}"

if [ -f "${COVFILE}" ];then
    python3 "${SRC_DIR}/merge_covariate.py" \
        --fam "${OUT_BASENAME}.${TEST_SUFFIX}.fam" \
        --cov "${COVFILE}" \
        --out "${OUT_BASENAME}.${TEST_SUFFIX}.cov"
fi