#!/bin/bash

# trap error
set -eo pipefail
err_report(){
  echo "Fail at $1: $2"
}
trap 'err_report ${LINENO} "${BASH_COMMAND}"' ERR


### arguments
while getopts 'hi:rg:m:b:f:c:C:d:o:' flag; do
    case $flag in
        h)
            echo "Predict PRS using adjusted beta (GenEpi is not available)"
            echo "options:"
            echo "-i, the input bfile prefix"
            echo "-r, dedup"
            echo "-g, the GenEpi prediction file"
            echo "-m, the method of PRS models; clf or reg"
            echo "-b, the beta file (beta.csv)"
            echo "-f, the reference rank file (rank_ref.csv)"
            echo "-c, the covariate file"
            echo "-C, the directory of covariate model"
            echo "-d, the directory of working and output"
            echo "-o, the basename of output file"
            ;;
        i) BFILE=$OPTARG;;
        r) DEDUP="true";;
        g) GENEPI_PRED=$OPTARG;;
        m) METHOD=$OPTARG;;
        b) MODEL=$OPTARG;;
        f) RANK=$OPTARG;;
        c) COV=$OPTARG;;
        C) COV_DIR=$OPTARG;;
        d) WORK_DIR=$OPTARG;;
        o) BASENAME=$OPTARG;;
        *) echo "usage: $0 [-i] [-r] [-g] [-m] [-b] [-f] [-c] [-C] [-d] [-o]"; exit 1;;
    esac
done

mkdir -p ${WORK_DIR}
REAL_PATH=$(realpath $0)
SRC_DIR=$(dirname ${REAL_PATH})


### predict
[ "${DEDUP}" = "true" ] && DEDUP_CMD="-r" || DEDUP_CMD=""
[ -f "${GENEPI_PRED}" ] && GENEPI_PRED_CMD="-g ${GENEPI_PRED}" || GENEPI_PRED_CMD=""
bash "${SRC_DIR}/predictPRS.sh" \
    -i "${BFILE}" \
    -b "${MODEL}" \
    -m "${METHOD}" \
    -d "${WORK_DIR}" \
    -o "${BASENAME}" ${DEDUP_CMD} ${GENEPI_PRED_CMD}


### analysis
[ -f "${COV}" ] && COV_CMD="--cov ${COV} --cov_ref_dir ${COV_DIR}" || COV_CMD=""
python3 "${SRC_DIR}/analysis.py" \
    --pred_file "${WORK_DIR}/prediction.csv" \
    --method "${METHOD}" \
    --mode "test" \
    --out_dir "${WORK_DIR}" \
    --rank_ref_file "${RANK}" ${COV_CMD}
