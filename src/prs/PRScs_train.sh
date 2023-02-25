#!/bin/bash

## pip3 install scipy h5py

# trap error
set -eo pipefail
err_report(){
  echo "Fail at $1: $2"
}
trap 'err_report ${LINENO} "${BASH_COMMAND}"' ERR

# arguments
while getopts 'hi:a:r:p:d:o:c:' flag; do
    case $flag in
        h)
            echo "PRScs training (the last character of directory path should not be '/')"
            echo "options:"
            echo "-i, the input bfile prefix"
            echo "-a, the filename of summary statistics"
            echo "-r, the directory of LD reference"
            echo "-p, the src of PRScs (PRScs.py)"
            echo "-d, the directory of working and output"
            echo "-o, the basename of output file"
            echo "-c, chromosome, default: '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22'"
            ;;
        i) BFILE=$OPTARG;;
        a) ASSOC_FILE=$OPTARG;;
        r) LD_REF_DIR=$OPTARG;;
        p) SRC=$OPTARG;;
        d) WORK_DIR=$OPTARG;;
        o) BASENAME=$OPTARG;;
        c) CHR=$OPTARG;;
        *) echo "usage: $0 [-i] [-a] [-r] [-p] [-d] [-o] [-c]"; exit 1;;
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

if [ ! -d "$LD_REF_DIR" ]; then
    echo "-r missing or designating error"
    exit 1
fi

if [ ! -f "$SRC" ]; then
    echo "-p missing or designating error"
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

# summary statistics
if [ ! -f "$WORK_DIR/summary_statistics_no_na.txt" ]; then
    tail -n +2 "$ASSOC_FILE" > "$WORK_DIR/summary_statistics.txt"
    echo -e "SNP\tA1\tA2\tOR\tP" > "$WORK_DIR/summary_statistics_no_na.txt"
    A1_COLNUM=$(sed -n '1s/\s/\n/gp' $ASSOC_FILE | grep -nx 'A1' | cut -d: -f1)
    A2_COLNUM=$(sed -n '1s/\s/\n/gp' $ASSOC_FILE | grep -nx 'A2' | cut -d: -f1)
    OR_COLNUM=$(sed -n '1s/\s/\n/gp' $ASSOC_FILE | grep -nx 'OR' | cut -d: -f1)
    P_COLNUM=$(sed -n '1s/\s/\n/gp' $ASSOC_FILE | grep -nx 'P' | cut -d: -f1)
    awk -v a1_colnum=$A1_COLNUM -v a2_colnum=$A2_COLNUM -v or_colnum=$OR_COLNUM -v p_colnum=$P_COLNUM '{if ($or_colnum != "NA") print $3"\t"$a1_colnum"\t"$a2_colnum"\t"$or_colnum"\t"$p_colnum}' "$WORK_DIR/summary_statistics.txt" >> "$WORK_DIR/summary_statistics_no_na.txt"
fi

# get sample size
OBS_CT_COLNUM=$(sed -n '1s/\s/\n/gp' $ASSOC_FILE | grep -nx 'OBS_CT' | cut -d: -f1)
SAMPLE_SIZE=$(awk -v max=0 -v obs_ct_colnum=$OBS_CT_COLNUM '{if($obs_ct_colnum>max){want=$obs_ct_colnum;max=$obs_ct_colnum}}END{print want}' "$WORK_DIR/summary_statistics.txt")

# keep autosomal chromosomes
plink2 \
    --bfile "${BFILE}" \
    --chr 1-22 \
    --allow-extra-chr \
    --make-bed \
    --out "${WORK_DIR}/${BASENAME}.auto"

# run PRScs in parallel
if [ -z "$CHR" ]; then
    CHR="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22"
fi
IFS=',' read -ra CHR_ARR <<< $CHR

for i in ${CHR_ARR[@]}
do
    python3 "$SRC" \
        --ref_dir="$LD_REF_DIR" \
        --bim_prefix="${WORK_DIR}/${BASENAME}.auto" \
        --sst_file="$WORK_DIR/summary_statistics_no_na.txt" \
        --n_gwas="$SAMPLE_SIZE" \
        --chrom="$i" \
        --phi=1e-2 \
        --n_iter=1000 \
        --n_burnin=500 \
        --out_dir="${WORK_DIR}/${BASENAME}"
done

# merge all
MERGE=true
for i in ${CHR_ARR[@]}
do
    if [ ! -f ${WORK_DIR}/${BASENAME}_*_chr${i}.txt ]; then
        MERGE=false
    fi
done

if $MERGE; then
    rm ${WORK_DIR}/effect_size.txt || true
    for i in ${CHR_ARR[@]}
    do
        cat ${WORK_DIR}/${BASENAME}_*_chr${i}.txt >> ${WORK_DIR}/effect_size.txt
    done
    rm ${WORK_DIR}/${BASENAME}.auto.* || true
fi
