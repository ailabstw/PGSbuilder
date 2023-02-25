#!/bin/bash

# trap error
set -eo pipefail
err_report(){
    echo "Fail at $1: $2"
}
trap 'err_report ${LINENO} "${BASH_COMMAND}"' ERR

# arguments
while getopts 'hs:b:t:T:B:c:p:o:C:m:E:' flag; do
    case $flag in
        h)
            echo "Split training and testing"
            echo "options:"
            echo "-s, the directory of source codes"
            echo "-b, bfile prefix"
            echo "-t, test list"
            echo "-T, target list"
            echo "-B, base list"
            echo "-c, covariate file"
            echo "-p, the test ratio (default=0.1)"
            echo "-o, output basename"
            echo "-C, Config file"
            ;;
        s) SRC_DIR=$OPTARG;;
        b) BFILE=$OPTARG;;
        t) TEST_LIST=$OPTARG;;
        T) TARGET_LIST=$OPTARG;;
        B) BASE_LIST=$OPTARG;;
        c) COVFILE=$OPTARG;;
        p) TEST_RATIO=$OPTARG;;
        o) PRS_OUT_DIR=$OPTARG;;
        C) CONFIG_FILE=$OPTARG;;
        m) METHOD=$OPTARG;;
        E) EXT_SUMSTAT_FILE=$OPTARG;;
        *) echo "usage: $0 [-s] [-a] [-b] [-c] [-v] [-m] [-d] [-r] [-p] [-o] [-i] [-j]"; exit 1;;
    esac
done

# redirect log
LOGDIR="${PRS_OUT_DIR}/LOG"
mkdir -p "${LOGDIR}"
LOGFILE=${LOGFILE-"${LOGDIR}/split.log"}
DETAIL_LOG=${DETAIL_LOG-"${LOGDIR}/split.detail.log"}
BASE_FLAG=${BASE_FLAG-"TRUE"}
[ "${BASE_FLAG}" = "TRUE" ] && exec &> >(tee -a "${LOGFILE}" )
[ "${BASE_FLAG}" = "TRUE" ] && exec &> >(tee -a "${DETAIL_LOG}")
BASE_FLAG="FALSE"
export DETAIL_LOG LOGFILE
#echo "SplitPipleline: LOG to ${LOGFILE}, detail log to ${DETAIL_LOG}"

# Check  
TEST_LIST=${TEST_LIST-NOT_PROVIDED}
TARGET_LIST=${TARGET_LIST-NOT_PROVIDED}
BASE_LIST=${BASE_LIST-NOT_PROVIDED}
COVFILE=${COVFILE-NOT_PROVIDED}
EXT_SUMSTAT_FILE=${EXT_SUMSTAT_FILE-NOT_PROVIDED}
TEST_RATIO=${TEST_RATIO-0.1}
BASE_RATIO=${BASE_RATIO-0.25}
[ ! -f "${BFILE}.bed" ] && { echo "SplitPipleline: BFILE ${BFILE} not found!"; TRHOW_AN_ERROR; }
[ ! -f "${BFILE}.fam" ] && { echo "SplitPipleline: BFILE ${BFILE} not found!"; TRHOW_AN_ERROR; }
[ ! -f "${BFILE}.bim" ] && { echo "SplitPipleline: BFILE ${BFILE} not found!"; TRHOW_AN_ERROR; }
[ ! -d "${SRC_DIR}" ] && { echo "SplitPipleline: SRC_DIR ${SRC_DIR} not found!"; TRHOW_AN_ERROR; }
[ ! -d "${PRS_OUT_DIR}" ] && { echo "SplitPipleline: dir ${PRS_OUT_DIR} not found!"; TRHOW_AN_ERROR; }
[ -f "${EXT_SUMSTAT_FILE}" ] && [ -f "${BASE_LIST}" ] && echo "SplitPipleline: BASE_LIST should not be provided when EXT_GWAS is present. BASE_LIST will have no effect"
[ -f "${EXT_SUMSTAT_FILE}" ] && [ -f "${TARGET_LIST}" ] && echo "SplitPipleline: TARGET_LIST should not be provided when EXT_GWAS is present. TARGET_LIST will have no effect"

DROP_NA='true'
RANDOM_SPLIT='true'
[ "$DROP_NA" = "true" ] && DROP_CMD="-d" || DROP_CMD=""
[ "$RANDOM_SPLIT" = "true" ] && RANDOM_CMD="-r" || RANDOM_CMD=""
IN_BASENAME=$(basename "${BFILE}")
mkdir -p "${PRS_OUT_DIR}/split"
mkdir -p "${PRS_OUT_DIR}/qc"


split_file_by_list(){
    ARG_BFILE="$1"
    ARG_LIST="$2"
    ARG_REMOVE="$3"
    ARG_OUT="$4"
    ARG_COV="$5"

    keep_cmd="I will be error"
    [ "${ARG_REMOVE}" =  "KEEP" ] && keep_cmd="--keep"
    [ "${ARG_REMOVE}" =  "REMOVE" ] && keep_cmd="--remove"
    plink2 --bfile "${ARG_BFILE}" "${keep_cmd}" "${ARG_LIST}" \
        --make-bed --out "${ARG_OUT}" >> ${DETAIL_LOG} 2>&1 || \
        { echo "SPLIT: split on ${ARG_BFILE} failed"; TRHOW_AN_ERROR; }
    
    if [ -f "${ARG_COV}" ];then
        python3 "${SRC_DIR}/split/merge_covariate.py" \
            --fam "${ARG_OUT}.fam" \
            --cov "${ARG_COV}" \
            --out "${ARG_OUT}.cov" >> ${DETAIL_LOG} 2>&1 || \
        { echo "SPLIT: merge_covariate on ${ARG_BFILE} failed"; TRHOW_AN_ERROR; }
    fi

    [[ $(wc -l <"${ARG_OUT}.fam") -ge 10 ]] || { echo "SplitPipleline: Less than 10 lines found in ${ARG_OUT}.fam, there may be some trouble your ind list"; TRHOW_AN_ERROR; }

}



#echo "SplitPipleline: Start prs split pipeline"
#echo "SplitPipleline: RAW DATA -> ${BFILE}"
#echo "SplitPipleline: EXT GWAS -> ${EXT_SUMSTAT_FILE}"
#echo "SplitPipleline: COV FILE -> ${COVFILE}"



# Split train test 
#echo "SplitPipleline: Split test/train from raw"
if [ -f "${TEST_LIST}" ]; then 
    echo "SPLIT: Split test/train by list"
    split_file_by_list "${BFILE}" "${TEST_LIST}" "KEEP" \
        "${PRS_OUT_DIR}/split/${IN_BASENAME}.test" "${COVFILE}" >> ${DETAIL_LOG} 2>&1 || \
        { echo "SPLIT: split test failed"; TRHOW_AN_ERROR; }

    split_file_by_list "${BFILE}" "${TEST_LIST}" "REMOVE" \
        "${PRS_OUT_DIR}/split/${IN_BASENAME}.train" "${COVFILE}"  >> ${DETAIL_LOG} 2>&1 || \
        { echo "SPLIT: split train failed"; TRHOW_AN_ERROR; }

else
    echo "SPLIT: Split test/train by random"
    bash ${SRC_DIR}/split/SplitTrainTest.sh \
        -s "${SRC_DIR}/split" \
        -a "${BFILE}.bed" \
        -b "${BFILE}.bim" \
        -c "${BFILE}.fam" \
        -v "${COVFILE}" \
        -m "${METHOD}" \
        -p ${TEST_RATIO} \
        -o "${PRS_OUT_DIR}/split/${IN_BASENAME}" \
        -i "train" \
        -j "test"  >> ${DETAIL_LOG} 2>&1 || \
        { echo "SPLIT: split train/test failed"; TRHOW_AN_ERROR; }
fi


# QC on train base
echo "SPLIT: Start QC on train"

bash ${SRC_DIR}/qc/gwas_qc.plink2.sh \
    -i "${PRS_OUT_DIR}/split/${IN_BASENAME}.train" \
    -o "${PRS_OUT_DIR}/qc/" \
    -C "${CONFIG_FILE}" >> ${DETAIL_LOG} 2>&1 || \
    { echo "SPLIT: QC on train failed"; TRHOW_AN_ERROR; }



# Split base target 
#echo "SplitPipleline: Split base/target from train after QC"
if [ -f "${EXT_SUMSTAT_FILE}" ]; then
    echo "SPLIT: Use external GWAS, base/target not split"
    if [ -f "${COVFILE}" ];then
        python3 "${SRC_DIR}/split/merge_covariate.py" \
            --fam "${PRS_OUT_DIR}/qc/${IN_BASENAME}.train.QC.fam" \
            --cov "${COVFILE}" \
            --out "${PRS_OUT_DIR}/qc/${IN_BASENAME}.train.QC.cov" >> ${DETAIL_LOG} 2>&1 || \
        { echo "SPLIT: QC on train failed"; TRHOW_AN_ERROR; }
    fi

elif [ -f "${TARGET_LIST}" ] && [ -f "${BASE_LIST}" ]; then 
    echo "SPLIT: Split base/target by base/target list"
    split_file_by_list "${PRS_OUT_DIR}/split/${IN_BASENAME}.train" "${BASE_LIST}" "KEEP" \
        "${PRS_OUT_DIR}/split/${IN_BASENAME}.train.QC.base" "${COVFILE}"  >> ${DETAIL_LOG} 2>&1 || \
        { echo "SPLIT: split base failed"; TRHOW_AN_ERROR; }

    split_file_by_list "${PRS_OUT_DIR}/split/${IN_BASENAME}.train" "${TARGET_LIST}" "KEEP" \
        "${PRS_OUT_DIR}/split/${IN_BASENAME}.train.QC.target" "${COVFILE}"  >> ${DETAIL_LOG} 2>&1 || \
        { echo "SPLIT: split target failed"; TRHOW_AN_ERROR; }

elif [ -f "${TARGET_LIST}" ] && [ ! -f "${BASE_LIST}" ]; then 
    echo "SPLIT: Split base/target by target list"
    split_file_by_list "${PRS_OUT_DIR}/split/${IN_BASENAME}.train" "${TARGET_LIST}" "REMOVE" \
        "${PRS_OUT_DIR}/split/${IN_BASENAME}.train.QC.base" "${COVFILE}"  >> ${DETAIL_LOG} 2>&1 || \
        { echo "SPLIT: split base failed"; TRHOW_AN_ERROR; }

    split_file_by_list "${PRS_OUT_DIR}/split/${IN_BASENAME}.train" "${TARGET_LIST}" "KEEP" \
        "${PRS_OUT_DIR}/split/${IN_BASENAME}.train.QC.target" "${COVFILE}"  >> ${DETAIL_LOG} 2>&1 || \
        { echo "SPLIT: split target failed"; TRHOW_AN_ERROR; }

elif [ ! -f "${TARGET_LIST}" ] && [ -f "${BASE_LIST}" ]; then 
    echo "SPLIT: Split base/target by base list"
    split_file_by_list "${PRS_OUT_DIR}/split/${IN_BASENAME}.train" "${BASE_LIST}" "KEEP" \
        "${PRS_OUT_DIR}/split/${IN_BASENAME}.train.QC.base" "${COVFILE}"  >> ${DETAIL_LOG} 2>&1 || \
        { echo "SPLIT: split base failed"; TRHOW_AN_ERROR; }

    split_file_by_list "${PRS_OUT_DIR}/split/${IN_BASENAME}.train" "${BASE_LIST}" "REMOVE" \
        "${PRS_OUT_DIR}/split/${IN_BASENAME}.train.QC.target" "${COVFILE}"  >> ${DETAIL_LOG} 2>&1 || \
        { echo "SPLIT: split target failed"; TRHOW_AN_ERROR; }

else
    echo "SPLIT: Split base/target by random"
    bash ${SRC_DIR}/split/SplitTrainTest.sh \
        -s "${SRC_DIR}/split" \
        -a "${PRS_OUT_DIR}/qc/${IN_BASENAME}.train.QC.bed" \
        -b "${PRS_OUT_DIR}/qc/${IN_BASENAME}.train.QC.bim" \
        -c "${PRS_OUT_DIR}/qc/${IN_BASENAME}.train.QC.fam" \
        -v "${COVFILE}" \
        -m "${METHOD}"  \
        -p ${BASE_RATIO} \
        -o "${PRS_OUT_DIR}/split/${IN_BASENAME}.train.QC" \
        -i "base" \
        -j "target"  >> ${DETAIL_LOG} 2>&1 || \
        { echo "SPLIT: split base/target failed"; TRHOW_AN_ERROR; }
fi



#echo "SplitPipleline: Get sample list"

SAMPLE_FILE="${PRS_OUT_DIR}/${IN_BASENAME}.sample.csv"
echo "FID,IID,dataset" > "${SAMPLE_FILE}"
awk '{ OFS=","; print $1,$2,"test" }' "${PRS_OUT_DIR}/split/${IN_BASENAME}.test.fam" >> "${SAMPLE_FILE}"
if [ -f "${EXT_SUMSTAT_FILE}" ]; then
    awk '{ OFS=","; print $1,$2,"target" }' "${PRS_OUT_DIR}/qc/${IN_BASENAME}.train.QC.fam" >> "${SAMPLE_FILE}"
else
    awk '{ OFS=","; print $1,$2,"base" }' "${PRS_OUT_DIR}/split/${IN_BASENAME}.train.QC.base.fam" >> "${SAMPLE_FILE}"
    awk '{ OFS=","; print $1,$2,"target" }' "${PRS_OUT_DIR}/split/${IN_BASENAME}.train.QC.target.fam" >> "${SAMPLE_FILE}"
fi
python3 - << EOF
import pandas as pd
df1 = pd.read_csv("${BFILE}.fam", sep='\s+', names=['FID','IID','FA','MA','SEX','PHE'])
df2 = pd.read_csv("${SAMPLE_FILE}")
merge_df = df1[['FID','IID']].merge(df2, on=['FID', 'IID'], how='left')
merge_df['dataset'] = merge_df['dataset'].fillna('removed')
merge_df.to_csv("${SAMPLE_FILE}", index=False)
EOF


echo "SPLIT: Split pipeline successfully complete"
#echo "SplitPipleline: TEST DATA -> ${TEST_DATA}"
#echo "SplitPipleline: BASE DATA -> ${BASE_DATA}"
#echo "SplitPipleline: TARGET DATA -> ${TARGET_DATA}"
export TEST_DATA BASE_DATA TARGET_DATA

