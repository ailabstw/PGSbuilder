#!/bin/bash

# trap error
set -eo pipefail
err_report(){
  echo "Fail at $1: $2"
}
trap 'err_report ${LINENO} "${BASH_COMMAND}"' ERR


### arguments
while getopts 'hb:t:v:C:c:o:a:m:d:' flag; do
    case $flag in
        h)
            echo "options:"
            echo "-b, the base bfile prefix (optional; only for GenEpi)"
            echo "-t, the target bfile prefix"
            echo "-v, the isolated validation bfile prefix"
            echo "-C, the config file"
            echo "-c, the covariate file of target data"
            echo "-o, the covariate file of test data"
            echo "-a, the summary statistics file of base data"
            echo "-m, the method of the PRS (reg or clf)"
            echo "-d, the output directory"
            ;;
        b) BASE=$OPTARG;;
        t) TARGET=$OPTARG;;
        v) TEST=$OPTARG;;
        C) CONFIG=$OPTARG;;
        c) TARGET_COV=$OPTARG;;
        o) TEST_COV=$OPTARG;;
        a) SS=$OPTARG;;
        m) METHOD=$OPTARG;;
        d) OUTDIR=$OPTARG;;
        *) echo "usage: $0 [-b] [-t] [-v] [-C] [-c] [-o] [-a] [-m] [-d]"; exit 1;;
    esac
done

if [ ! -f "$BASE.bed" ] || [ ! -f "$BASE.bim" ] || [ ! -f "$BASE.fam" ]; then
    echo "PRS: no base data"
    RUN_BASE="false"
else
    RUN_BASE="true"
fi

if [ ! -f "$TARGET.bed" ] || [ ! -f "$TARGET.bim" ] || [ ! -f "$TARGET.fam" ]; then
    echo "PRS: -t missing or designating error"
    exit 1
fi

if [ ! -f "$TEST.bed" ] || [ ! -f "$TEST.bim" ] || [ ! -f "$TEST.fam" ]; then
    RUN_TEST="false"
else
    RUN_TEST="true"
fi

if [ ! -f $SS ]; then
    echo "PRS: -a missing or designating error"
    exit 1
fi

if [ -z "$METHOD" ]; then
    echo "PRS: -m missing; assign reg or clf"
    exit 1
fi

if [ -z "$OUTDIR" ]; then
    echo "PRS: -d missing"
    exit 1
fi

if [ -z "$TOOLS" ]; then
    TOOLS="CandT,PRSice2,Lassosum,LDpred2,PRScs,GenEpi"
fi

# Basename
if [ "$RUN_BASE" = "true" ]; then
    BASE_BASENAME=$(basename "$BASE")
fi
TARGET_BASENAME=$(basename "$TARGET")
if [ "$RUN_TEST" = "true" ]; then
    TEST_BASENAME=$(basename "$TEST")
fi

source "${CONFIG}"
REAL_PATH=$(realpath $0)
SRC_DIR=$(dirname ${REAL_PATH})
mkdir -p "${OUTDIR}"
mkdir -p "${OUTDIR}/prediction"
mkdir -p "${OUTDIR}/analysis"
cd ${OUTDIR} || exit


##### ADD by YILUN 20220615 for try and catch
LOGDIR="${OUTDIR}/LOG"
mkdir -p "${LOGDIR}"
LOGFILE="${LOGDIR}/prs.log"
DETAIL_LOG=${DETAIL_LOG-"${LOGDIR}/prs.detail.log"}
exec &> >(tee "${LOGFILE}")
exec &> >(tee -a "${DETAIL_LOG}")

# https://stackoverflow.com/questions/22009364/is-there-a-try-catch-command-in-bash
function my_try(){
  # Check if out bash is set -e or +e and save the condition
  [[ $- = *e* ]]; SAVED_OPT_E=$?
  # turn off set -e 
  set +e
}

function my_catch(){
  EXIT_CODE=$?
  ALGO_NAME=$1
  if [ "${EXIT_CODE}" = 0 ];then 
    echo "PRS: ${ALGO_NAME} successfully complete"
  else 
    echo "PRS: ${ALGO_NAME} failed"
    rm -rf "${OUTDIR}/${ALGO_NAME}" || true
  fi
  # return back to origin set e condidctionm
  [ ${SAVED_OPT_E} = 0 ] && set -e
}


### PRS models
#echo "==========================================================="
#printf "Building PRS Models\n"
#echo "==========================================================="

# Clumping and Thresholding
if [[ ${TOOLS} =~ "CandT" ]]; then
my_try
(   
    set -e
    printf "\n###### Clumping and Thresholding ######\n"
    mkdir -p "$OUTDIR/CandT"

    # target
    bash "${SRC_DIR}/clump_threshold_train.sh" \
        -i "$TARGET" \
        -a $SS \
        -m "$METHOD" \
        -p "$SRC_DIR" \
        -d "$OUTDIR/CandT" \
        -o "$TARGET_BASENAME"

) 2>&1  | tee ${LOGDIR}/CandT.log >> "${DETAIL_LOG}"
my_catch "CandT"
fi

# PRSice2
if [[ ${TOOLS} =~ "PRSice2" ]]; then
my_try
(   
    set -e
    printf "\n###### PRSice2 ######\n"
    mkdir -p "$OUTDIR/PRSice2"

    # target
    bash "${SRC_DIR}/PRSice2_train.sh" \
        -r "$TARGET" \
        -a $SS \
        -m "$METHOD" \
        -d "$OUTDIR/PRSice2" \
        -o "$TARGET_BASENAME"
) 2>&1  | tee ${LOGDIR}/PRSice2.log >> "${DETAIL_LOG}"
my_catch "PRSice2"
fi


# Lassosum
if [[ ${TOOLS} =~ "Lassosum" ]]; then
my_try
(   
    set -e
    printf "\n###### Lassosum ######\n"
    mkdir -p "$OUTDIR/Lassosum"

    # target
    Rscript "${SRC_DIR}/lassosum_train.R" \
        -i "$TARGET" \
        -a $SS \
        -l "$POPULATION_PRS" \
        -g "$GENOME" \
        -d "$OUTDIR/Lassosum" \
        -o "$TARGET_BASENAME"

    # remove temp files
    rm "./Rplots.pdf" || true
) 2>&1  | tee ${LOGDIR}/Lassosum.log >> "${DETAIL_LOG}"
my_catch "Lassosum"
fi


# LDpred2
if [[ ${TOOLS} =~ "LDpred2" ]]; then
my_try
(   
    set -e
    printf "\n###### LDpred2 ######\n"
    mkdir -p "$OUTDIR/LDpred2"

    # mode
    P_COLNUM=$(sed -n '1s/\s/\n/gp' $SS | grep -nx 'P' | cut -d: -f1)
    SIG_COUNT=$(awk -v n="$P_COLNUM" '$n<1e-8{c++}END{print c}' $SS)
    [ "$SIG_COUNT" -gt 10 ] && MODE=4 || MODE=1

    # target
    Rscript "${SRC_DIR}/ldpred2_train.R" \
        -i "$TARGET" \
        -a $SS \
        -d "$OUTDIR/LDpred2" \
        -o "$TARGET_BASENAME" \
        -l "$LIFTOVER_REF_DIR" \
        -b "$LDPRED_REF_DIR" \
        -g "$GENOME" \
        -t "$MODE"

    # remove temp files
    rm "${OUTDIR}/LDpred2/${TARGET_BASENAME}.bk" || true
    rm "${OUTDIR}/LDpred2/${TARGET_BASENAME}.rds" || true

) 2>&1  | tee ${LOGDIR}/LDpred2.log >> "${DETAIL_LOG}"
my_catch "LDpred2"
fi


# PRScs: must be rsID
if [[ ${TOOLS} =~ "PRScs" ]]; then
my_try
(   
    set -e
    printf "\n###### PRScs ######\n"
    mkdir -p "$OUTDIR/PRScs"


    ####### for test
    if [ "$POPULATION_PRS" == "ASN" ];then
        POPULATION_LOWER="eas"
    else
        POPULATION_LOWER=$(echo "$POPULATION_PRS" | awk '{print tolower($0)}')
    fi

    # target
    bash "${SRC_DIR}/PRScs_train.sh" \
        -i "$TARGET" \
        -a $SS \
        -r "$PRSCS_REF_DIR/ldblk_1kg_$POPULATION_LOWER" \
        -p "$PRSCS_SRC" \
        -d "$OUTDIR/PRScs" \
        -o "$TARGET_BASENAME"
) 2>&1  | tee ${LOGDIR}/PRScs.log >> "${DETAIL_LOG}"
my_catch "PRScs"
fi


# GenEpi
if [[ ${TOOLS} =~ "GenEpi" ]] && [ "$RUN_BASE" = "true" ]; then
my_try
(   
    set -e
    printf "\n###### GenEpi ######\n"
    if [ "${METHOD}" == "clf" ]; then
        METHOD_CODE="c"
        METHOD_BASENAME="Classifier"
    else
        METHOD_CODE="r"
        METHOD_BASENAME="Regressor"
    fi

    ## target
    mkdir -p "${OUTDIR}/GenEpi"
    bash "${SRC_DIR}/genepi_train.sh" \
        -i "${BASE}" \
        -m "${METHOD_CODE}" \
        -b "${GENEPI_REF_DIR}" \
        -g "${GENOME}" \
        -d "${OUTDIR}/GenEpi" \
        -o "${BASE_BASENAME}" \
        -e
    
    bash "${SRC_DIR}/genepi_test.sh" \
        -i "${TARGET}" \
        -t "${METHOD_CODE}" \
        -m "${OUTDIR}/GenEpi/${BASE_BASENAME}/crossGeneResult/${METHOD_BASENAME}.pkl" \
        -f "${OUTDIR}/GenEpi/${BASE_BASENAME}/crossGeneResult/Feature.csv" \
        -p "${SRC_DIR}/GenEpi_predictor.py" \
        -d "${OUTDIR}/GenEpi" \
        -o "${TARGET_BASENAME}"
    
    ## testing
    if [ "$RUN_TEST" = "true" ]; then
        bash "${SRC_DIR}/genepi_test.sh" \
            -i "${TEST}" \
            -t "${METHOD_CODE}" \
            -m "${OUTDIR}/GenEpi/${BASE_BASENAME}/crossGeneResult/${METHOD_BASENAME}.pkl" \
            -f "${OUTDIR}/GenEpi/${BASE_BASENAME}/crossGeneResult/Feature.csv" \
            -p "${SRC_DIR}/GenEpi_predictor.py" \
            -d "${OUTDIR}/GenEpi" \
            -o "${TEST_BASENAME}"
    fi
) 2>&1  | tee ${LOGDIR}/GenEpi.log >> "${DETAIL_LOG}"
my_catch "GenEpi"
fi

{
# check and merge beta
cd ${SRC_DIR} || exit
SS_STR=$(ls ${SS})
python3 ${SRC_DIR}/CollectBeta.py \
    -t "${TARGET}" \
    -s "${SS_STR}" \
    -o "${OUTDIR}" \
    -a "${TOOLS}"

# Get freq 
plink1.9 \
    --bfile "${TARGET}" \
    --freq \
    --allow-extra-chr \
    --allow-no-sex \
    --out "${TARGET}"
cd ${OUTDIR} || exit

# remove failed algo 
for ALGO_NAME in $(cat ${LOGDIR}/fail_algo_list);do
    [ -d "${OUTDIR}/${ALGO_NAME}" ] && rm -rf "${OUTDIR}/${ALGO_NAME}" || true
    echo "Remove ${OUTDIR}/${ALGO_NAME} since $ALGO_NAME failed"
done
}  >> ${DETAIL_LOG} 2>&1 || \
        { echo "PRS: Collect beta failed"; TRHOW_AN_ERROR; }

echo "PRS: Build PRS model complete"


{
### prediction
echo "==========================================================="
printf "Predicting Target and Test Sets \n"
echo "==========================================================="

printf "###### Predicting Target ######\n"
mkdir -p ${OUTDIR}/prediction/target
mkdir -p ${OUTDIR}/analysis/target

bash ${SRC_DIR}/predictPRS.sh \
    -i "${TARGET}" \
    -b "${OUTDIR}/beta.tsv" \
    -m "${METHOD}" \
    -d "${OUTDIR}/prediction/target" \
    -o "${TARGET_BASENAME}" \
    -r \
    -g "${OUTDIR}/GenEpi/${TARGET_BASENAME}.pred.csv"

mv "${OUTDIR}/prediction/target/prediction.csv" "${OUTDIR}/analysis/target/prediction.csv"

if [ "$RUN_BASE" = "true" ]; then
    printf "###### Predicting Base ######\n"
    mkdir -p ${OUTDIR}/prediction/base
    mkdir -p ${OUTDIR}/analysis/base

    bash ${SRC_DIR}/predictPRS.sh \
        -i "${BASE}" \
        -b "${OUTDIR}/beta.tsv" \
        -m "${METHOD}" \
        -d "${OUTDIR}/prediction/base" \
        -o "${BASE_BASENAME}" \
        -r \
        -g "${OUTDIR}/GenEpi/${BASE_BASENAME}.pred.csv"

    mv "${OUTDIR}/prediction/base/prediction.csv" "${OUTDIR}/analysis/base/prediction.csv"
fi


if [ "$RUN_TEST" = "true" ]; then
    printf "###### Predicting Test ######\n"
    mkdir -p ${OUTDIR}/prediction/test
    mkdir -p ${OUTDIR}/analysis/test

    bash ${SRC_DIR}/predictPRS.sh \
        -i "${TEST}" \
        -b "${OUTDIR}/beta.tsv" \
        -m "${METHOD}" \
        -d "${OUTDIR}/prediction/test" \
        -o "${TEST_BASENAME}" \
        -r \
        -g "${OUTDIR}/GenEpi/${TEST_BASENAME}.pred.csv"

    mv "${OUTDIR}/prediction/test/prediction.csv" "${OUTDIR}/analysis/test/prediction.csv"
fi

printf "###### Predicting Target and Test Sets Complete ######\n"



### analysis: cohort reference, covariates, performance
echo "==========================================================="
printf "Analyzing PRS Results \n"
echo "==========================================================="

printf "###### Analyzing Target ######\n"
[ -f "${TARGET_COV}" ] && TARGET_COV_CMD="--cov ${TARGET_COV}" || TARGET_COV_CMD=""
python3 ${SRC_DIR}/analysis.py \
    --pred_file "${OUTDIR}/analysis/target/prediction.csv" \
    --method "${METHOD}" \
    --mode "target" \
    --out_dir "${OUTDIR}/analysis/target" ${TARGET_COV_CMD} \
    --run_performance

mv "${OUTDIR}/analysis/target/rank_ref.csv" "${OUTDIR}/rank_ref.csv"
mv "${OUTDIR}/analysis/target/hist_ref.csv" "${OUTDIR}/hist_ref.csv"


if [ "$RUN_TEST" = "true" ]; then
    printf "###### Analyzing Test ######\n"
    [ -f "${TEST_COV}" ] && TEST_COV_CMD="--cov ${TEST_COV}" || TEST_COV_CMD=""
    python3 ${SRC_DIR}/analysis.py \
        --pred_file "${OUTDIR}/analysis/test/prediction.csv" \
        --method "${METHOD}" \
        --mode "test" \
        --out_dir "${OUTDIR}/analysis/test" \
        --rank_ref_file "${OUTDIR}/rank_ref.csv" \
        --cov_ref_dir "${OUTDIR}/analysis/target/cov" ${TEST_COV_CMD} \
        --run_performance
fi


if [ "$RUN_BASE" = "true" ]; then
    printf "###### Analyzing Base ######\n"
    BASE_COV="${BASE}.cov"
    [ -f "${BASE_COV}" ] && BASE_COV_CMD="--cov ${BASE_COV}" || BASE_COV_CMD=""
    python3 ${SRC_DIR}/analysis.py \
        --pred_file "${OUTDIR}/analysis/base/prediction.csv" \
        --method "${METHOD}" \
        --mode "test" \
        --out_dir "${OUTDIR}/analysis/base" \
        --rank_ref_file "${OUTDIR}/rank_ref.csv" \
        --cov_ref_dir "${OUTDIR}/analysis/target/cov" ${BASE_COV_CMD} \
        --run_performance
fi



#printf "###### Analyzing PRS Performance Complete ######\n"
}  >> ${DETAIL_LOG} 2>&1 || \
        { echo "PRS: Prediction and evaluation failed"; TRHOW_AN_ERROR; }





echo "PRS: Prediction and evaluation complete"
