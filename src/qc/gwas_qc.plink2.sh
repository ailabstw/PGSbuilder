#!/bin/bash

# trap error
set -eo pipefail
err_report(){
  echo "Fail at $1: $2"
}
trap 'err_report ${LINENO} "${BASH_COMMAND}"' ERR

# arguments
while getopts 'hi:o:C:' flag; do
    case $flag in
        h)
            echo "Pipeline for GWAS QC (absolute paths are required)"
            echo "options:"
            echo "-i, the basename of input bfile"
            echo "-o, the directory of output"
            echo "-C, the path to config file"
            ;;
        i) IN_FILENAME=$OPTARG;;
        o) OUT_DIR=$OPTARG;;
        C) CONFIG=$OPTARG;;
        *) echo "usage: $0 [-i] [-o] [-C]"; exit 1;;
    esac
done

if [ -z "$CONFIG" ]; then
    echo "-C missing"
    exit 1
fi

if [ ! -f "$IN_FILENAME.bed" ] || [ ! -f "$IN_FILENAME.bim" ] || [ ! -f "$IN_FILENAME.fam" ]; then
    echo "-i missing or designating error"
    exit 1
fi

if [ -z "$OUT_DIR" ]; then
    echo "-o missing"
    exit 1
fi


# preparation
source "${CONFIG}"
REAL_PATH=$(realpath $0)
SRC_DIR=$(dirname ${REAL_PATH})
IN_BASENAME=$(basename "${IN_FILENAME}")
HM_BASENAME=$(basename "${HM_FILENAME}")
mkdir -p "$OUT_DIR"
cd "$OUT_DIR" || exit
mkdir -p QualityControl
mkdir -p PopulationStratification


##########################################################
########################### QC ###########################
##########################################################
cd QualityControl || exit

echo "==========================================================="
printf "Perfoming QC Step 1 -- Dedup & Valid Chrom\n"
echo "==========================================================="
plink2 \
    --bfile "${IN_FILENAME}" \
    --rm-dup force-first \
    --chr 1-22,X \
    --allow-extra-chr \
    --make-bed \
    --write-snplist \
    --out "${IN_BASENAME}.dedup"
mv "${IN_BASENAME}.dedup.snplist" "${IN_BASENAME}.QC.1.dedup.snplist" 


echo "==========================================================="
printf "Perfoming QC Step 2 -- SNP Filter\n"
echo "==========================================================="
# allele frequency, missing rate
{
plink2 \
    --bfile "${IN_BASENAME}.dedup" \
    --allow-extra-chr \
    --freq \
    --missing variant-only \
    --out "${IN_BASENAME}.QC"
} 2>&1 | tee "${IN_BASENAME}.QC.2.snp.1.log"

# H-W, flip scan
[ "${FLIPSCAN}" = "true" ] && FLIPSCAN_CMD="--flip-scan" || FLIPSCAN_CMD=""
{
plink1.9 \
    --bfile "${IN_BASENAME}.dedup" \
    --allow-extra-chr \
    --allow-no-sex \
    --hardy ${FLIPSCAN_CMD} \
    --out "${IN_BASENAME}.QC"
} 2>&1 | tee "${IN_BASENAME}.QC.2.snp.2.log"

# SNP report
python3 "${SRC_DIR}/Report/SNP_QC_report.py" \
    -i "${IN_BASENAME}.QC" \
    -C "${CONFIG}" \
    -o "${IN_BASENAME}.QC"


echo "==========================================================="
printf "Perfoming QC Step 3 -- Sample Filter\n"
echo "==========================================================="
# Missing, Heterozygosity
{
plink2 \
    --bfile "${IN_BASENAME}.dedup" \
    --allow-extra-chr \
    --extract "${IN_BASENAME}.QC.snp_qc.snp_list" \
    --missing sample-only \
    --read-freq "${IN_BASENAME}.QC.afreq" \
    --het \
    --out "${IN_BASENAME}.QC"
} 2>&1 | tee "${IN_BASENAME}.QC.3.ind.1.log"

# Sex
if [ "${CHECK_SEX}" = "true" ]; then
    {
    plink1.9 \
        --bfile "${IN_BASENAME}.dedup" \
        --allow-extra-chr \
        --extract "${IN_BASENAME}.QC.snp_qc.snp_list" \
        --check-sex \
        --out "${IN_BASENAME}.QC" || echo "Check sex fail, possibly there is no chrX"
    } 2>&1 | tee "${IN_BASENAME}.QC.3.ind.2.log"
fi

# Individual report
python3 "${SRC_DIR}/Report/IND_QC_report.py" \
    -i "${IN_BASENAME}.QC" \
    -f "${IN_FILENAME}.fam" \
    -C "${CONFIG}" \
    -o "${IN_BASENAME}.QC"


echo "==========================================================="
printf "Perfoming QC Step 4 -- Kinship Filter\n"
echo "==========================================================="
{
if [ $KINGCUTOFF = "0" ]; then
    ln -s "${IN_BASENAME}.QC.ind_qc.ind_list" "${IN_BASENAME}.QC.king.cutoff.in.id"
    touch ${IN_BASENAME}.QC.king.cutoff.out.id 
    echo "KINGCUTOFF is set to 0, skip prune and kingship filter"

else
    plink2 \
        --bfile "${IN_BASENAME}.dedup" \
        --allow-extra-chr \
        --extract "${IN_BASENAME}.QC.snp_qc.snp_list" \
        --keep "${IN_BASENAME}.QC.ind_qc.ind_list" \
        --indep-pairwise "${PRUNE_WINDOW}" "${PRUNE_STEP}" "${PRUNE_THRESHOLD}" \
        --read-freq "${IN_BASENAME}.QC.afreq" \
        --memory "${MEMORY}" \
        --threads "${THREAD}" \
        --out "${IN_BASENAME}.QC" 2>&1 | tee "${IN_BASENAME}.QC.4.prune.log"

    plink2 \
        --bfile "${IN_BASENAME}.dedup" \
        --read-freq "${IN_BASENAME}.QC.afreq" \
        --allow-extra-chr \
        --extract "${IN_BASENAME}.QC.prune.in" \
        --keep "${IN_BASENAME}.QC.ind_qc.ind_list" \
        --king-cutoff "${KINGCUTOFF}" \
        --memory "${MEMORY}" \
        --threads "${THREAD}" \
        --out "${IN_BASENAME}.QC"
fi 
} 2>&1 | tee "${IN_BASENAME}.QC.4.kinship.log"


echo "==========================================================="
printf "Perfoming QC Step 5 -- Making post-QC bfile\n"
echo "==========================================================="
{
plink2 \
    --bfile "${IN_BASENAME}.dedup" \
    --allow-extra-chr \
    --make-bed \
    --extract "${IN_BASENAME}.QC.snp_qc.snp_list" \
    --keep "${IN_BASENAME}.QC.king.cutoff.in.id" \
    --keep-allele-order \
    --read-freq "${IN_BASENAME}.QC.afreq" \
    --out "${IN_BASENAME}.QC"
} 2>&1 | tee "${IN_BASENAME}.QC.5.bfile.log"


##########################################################
#################Population Stratification################
##########################################################
cd ../PopulationStratification || exit

ln -s "${HM_FILENAME}.bed" "${HM_BASENAME}.bed" || echo "${HM_BASENAME}.bed exits"
ln -s "${HM_FILENAME}.fam" "${HM_BASENAME}.fam" || echo "${HM_BASENAME}.fam exits"
awk 'NR==FNR {A[$1 $4]=$2; next;}; {if(($1 $4) in A) K=A[$1 $4]; else K=$2}; { print $1 OFS K OFS $3 OFS $4 OFS $5 OFS $6}' OFS='\t' "../QualityControl/${IN_BASENAME}.QC.bim" "${HM_FILENAME}.bim" > "${HM_BASENAME}.bim"

echo "==========================================================="
printf "Population Stratification Step 1 -- QC for Hapmap\n"
echo "==========================================================="
{
plink2 \
    --bfile "${HM_BASENAME}" \
    --rm-dup force-first \
    --make-bed \
    --allow-extra-chr \
    --chr 1-22 \
    --out "${HM_BASENAME}.dedup"
} 2>&1 | tee "${HM_BASENAME}.PS.1.dedup.log"

[ "${STRICT_HWE}" = "true" ] && HWE_QC_HAPMAP=${HWE_QC_CONTROL} || HWE_QC_HAPMAP=${HWE_QC}
{
plink2 \
    --bfile "${HM_BASENAME}.dedup" \
    --maf "${MAF_QC}" \
    --hwe ${HWE_QC_HAPMAP} \
    --geno "${GENO_QC}" \
    --mind "${MIND_QC}" \
    --write-snplist \
    --make-just-fam \
    --out "${HM_BASENAME}.PS"
} 2>&1 | tee "${HM_BASENAME}.PS.1.qc.log"


echo "==========================================================="
printf "Population Stratification Step 2 -- Merging Input and Hapmap\n"
echo "==========================================================="
awk '{print $2}' "../QualityControl/${IN_BASENAME}.QC.bim" | grep -F -x -f "${HM_BASENAME}.PS.snplist" > "common_rsid.txt"

{
plink2 \
    --bfile "../QualityControl/${IN_BASENAME}.QC" \
    --extract "common_rsid.txt" \
    --make-just-bim \
    --allow-no-sex \
    --out "${IN_BASENAME}.PS"
} 2>&1 | tee "${IN_BASENAME}.PS.2.1.comm.log"

{
plink2 \
    --bfile "${HM_BASENAME}.dedup" \
    --keep "${HM_BASENAME}.PS.fam" \
    --extract "common_rsid.txt" \
    --update-map "${IN_BASENAME}.PS.bim" 4 2 \
    --make-just-bim \
    --out "${HM_BASENAME}.PS"
} 2>&1 | tee "${HM_BASENAME}.PS.2.1.comm.log"

# Resolve strand issues
## get difference between input and heatmap
awk '{print $2,$5,$6}' "${IN_BASENAME}.PS.bim" > "${IN_BASENAME}.snp.txt"
awk '{print $2,$5,$6}' "${HM_BASENAME}.PS.bim" > "${HM_BASENAME}.snp.txt"
sort "${HM_BASENAME}.snp.txt" "${IN_BASENAME}.snp.txt" | uniq -u | awk '{print $1}' | sort -u > "flip.snp.txt"

## flip SNPs for resolving difference
{
plink1.9 \
    --bfile "../QualityControl/${IN_BASENAME}.QC" \
    --extract "common_rsid.txt" \
    --flip "flip.snp.txt" \
    --make-just-bim \
    --allow-no-sex \
    --out "${IN_BASENAME}.PS"
} 2>&1 | tee "${IN_BASENAME}.PS.2.2.flip.log"

## get unresolved SNPs
awk '{print $2,$5,$6}' "${IN_BASENAME}.PS.bim" > "${IN_BASENAME}.snp.txt"
sort "${HM_BASENAME}.snp.txt" "${IN_BASENAME}.snp.txt" | uniq -u | awk '{print $1}' | sort -u > "exclusion.snp.txt"

## exclude unresolved SNPs
{
plink1.9 \
    --bfile "../QualityControl/${IN_BASENAME}.QC" \
    --extract "common_rsid.txt" \
    --exclude "exclusion.snp.txt" \
    --flip "flip.snp.txt" \
    --allow-no-sex \
    --keep-allele-order \
    --make-bed \
    --out "${IN_BASENAME}.PS"
} 2>&1 | tee "${IN_BASENAME}.PS.2.3.unsol.log"

{
plink2 \
    --bfile "${HM_BASENAME}.dedup" \
    --keep "${HM_BASENAME}.PS.fam" \
    --extract "common_rsid.txt" \
    --exclude "exclusion.snp.txt" \
    --update-map "${IN_BASENAME}.PS.bim" 4 2 \
    --keep-allele-order \
    --make-bed \
    --out "${HM_BASENAME}.PS"
} 2>&1 | tee "${HM_BASENAME}.PS.2.3.unsol.log"

# Merge input and heatmap
{
plink1.9 \
    --bfile "${IN_BASENAME}.PS" \
    --bmerge "${HM_BASENAME}.PS" \
    --allow-no-sex \
    --keep-allele-order \
    --make-bed \
    --out "merge.PS"
} 2>&1 | tee "merge.PS.2.4.merge.log"


echo "==========================================================="
printf "Population Stratification Step 3 -- Performing PCA\n"
echo "==========================================================="
{
plink2 \
    --bfile "${IN_BASENAME}.PS" \
    --indep-pairwise "${PRUNE_WINDOW}" "${PRUNE_STEP}" "${PRUNE_THRESHOLD}" \
    --memory "${MEMORY}" \
    --threads "${THREAD}" \
    --out "${IN_BASENAME}.PS"
} 2>&1 | tee "${IN_BASENAME}.PS.3.1.prune.log"

PCA_COUNT=3
SAMPLE_NUM=$(wc -l merge.PS.fam | cut -d" " -f1)
[ ${SAMPLE_NUM} -gt 5000 ] && PCA_CMD="--pca ${PCA_COUNT} approx" || PCA_CMD="--pca ${PCA_COUNT}"
{
plink2 \
    --bfile "merge.PS" \
    --extract "${IN_BASENAME}.PS.prune.in" \
    ${PCA_CMD} \
    --memory "${MEMORY}" \
    --threads "${THREAD}" \
    --out "merge.PS"
} 2>&1 | tee "merge.PS.3.2.pca.log"

Rscript "${SRC_DIR}/PCA.R" \
    -e "merge.PS.eigenvec" \
    -v "merge.PS.eigenval" \
    -f "merge.PS.fam" \
    -i "${HM_INFO}" \
    -p "${POPULATION}" \
    -n "${PCA_COUNT}" \
    -m "SD" \
    -s "${POP_SD}" \
    -o "merge.PS" 2>&1 | tee "merge.PS.3.3.pca.log"


echo "==========================================================="
printf "Population Stratification Step 4 -- Removing Failed Samples\n"
echo "==========================================================="
plink2 \
    --bfile "../QualityControl/${IN_BASENAME}.QC" \
    --keep "merge.PS.pca_qc.ind_list" \
    --keep-allele-order \
    --make-bed \
    --out "../${IN_BASENAME}.QC"


echo "==========================================================="
printf "QC Report\n"
echo "==========================================================="
cd ../
python3 "${SRC_DIR}/QC_report.py" \
    --bfile "${IN_FILENAME}" \
    --dedup_record "QualityControl/${IN_BASENAME}.QC.1.dedup.snplist" \
    --snp_record "QualityControl/${IN_BASENAME}.QC.snp_qc.json" \
    --ind_record "QualityControl/${IN_BASENAME}.QC.ind_qc.json" \
    --kinship_record "QualityControl/${IN_BASENAME}.QC.king.cutoff.out.id" \
    --population_record "PopulationStratification/merge.PS.pca_qc.ind_list" \
    --fig_maf "QualityControl/${IN_BASENAME}.QC.maf.hist.png" \
    --fig_vmiss "QualityControl/${IN_BASENAME}.QC.geno.hist.png" \
    --fig_hwe "QualityControl/${IN_BASENAME}.QC" \
    --fig_het_vs_smiss "QualityControl/${IN_BASENAME}.QC.hetxmiss.png" \
    --fig_population_prefix "PopulationStratification/merge.PS" \
    --maf "${MAF_QC}" \
    --vmiss "${GENO_QC}" \
    --hwe_all "${HWE_QC}" \
    --hwe_strict "${STRICT_HWE}" \
    --hwe_control "${HWE_QC_CONTROL}" \
    --hwe_case "${HWE_QC_CASE}" \
    --flipscan "${FLIPSCAN}" \
    --heter "${HETER_SD}" \
    --smiss "${MIND_QC}" \
    --sex "${CHECK_SEX}" \
    --kinship "${KINGCUTOFF}" \
    --prune_window "${PRUNE_WINDOW}" \
    --prune_step "${PRUNE_STEP}" \
    --prune_threshold "${PRUNE_THRESHOLD}" \
    --population "${POPULATION}" \
    --population_sd "${POP_SD}" \
    --work_dir ./

