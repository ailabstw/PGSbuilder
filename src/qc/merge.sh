#!/bin/bash

# trap error
set -eo pipefail
err_report(){
  echo "Fail at $1: $2"
}
trap 'err_report ${LINENO} "${BASH_COMMAND}"' ERR

# arguments
while getopts ':h:c:l:o:f:' flag; do
    case $flag in
        h)
            echo "Pipeline for merging CASE CONTROL (absolute paths are required)"
            echo "options:"
            echo "-c, the basename of case bfile"
            echo "-l, the basename of control bfile"
            echo "-o, the directory of output"
            echo "-f, the output filename"
            ;;
        c) IN_CASE=$OPTARG;;
        l) IN_CONTROL=$OPTARG;;
        o) OUT_DIR=$OPTARG;;
        f) OUT_FILE=$OPTARG;;
        *) echo "usage: $0 [-c] [-l] [-o] [-f]"; exit 1;;
    esac
done

CASE_BED=$IN_CASE.bed
CASE_BIM=$IN_CASE.bim
CASE_FAM=$IN_CASE.fam

CONTROL_BED=$IN_CONTROL.bed
CONTROL_BIM=$IN_CONTROL.bim
CONTROL_FAM=$IN_CONTROL.fam

if [ ! -f "$CASE_BED" ] || [ ! -f "$CASE_BIM" ] || [ ! -f "$CASE_FAM" ]; then
    echo "-case missing or designating error"
    exit 1
fi

if [ ! -f "$CONTROL_BED" ] || [ ! -f "$CONTROL_BIM" ] || [ ! -f "$CONTROL_FAM" ]; then
    echo "-control missing or designating error"
    exit 1
fi

if [ -z "$OUT_DIR" ]; then
    echo "-o missing"
    exit 1
fi

# input dir
if [[ $IN_CASE =~ ([^[:space:]]+)/ ]]; then IN_CASE_DIR=${BASH_REMATCH[1]}; fi
if [[ $IN_CONTROL =~ ([^[:space:]]+)/ ]]; then IN_CONTROL_DIR=${BASH_REMATCH[1]}; fi
# Set Working Directory
cd $OUT_DIR
mkdir -p WorkingDir
cd WorkingDir
echo "==========================================================="
printf "STEP 1: Preparing case and control data\n"
echo "==========================================================="
printf "Finding duplicate SNPs ID\n"
printf ".\n.\n.\n"
cat $IN_CASE.bim | awk '{print $2}' > all_snps.snplist
cat all_snps.snplist | sort | uniq -d > duplicated_snps.snplist
plink1.9 --bfile $IN_CASE \
         --exclude duplicated_snps.snplist \
		 --keep-allele-order \
         --make-bed \
         --out unique_CASE \
         --allow-no-sex

cat $IN_CONTROL.bim | awk '{print $2}' > all_snps.snplist
cat all_snps.snplist | sort | uniq -d > duplicated_snps.snplist
plink1.9 --bfile $IN_CONTROL \
         --exclude duplicated_snps.snplist \
		 --keep-allele-order \
         --make-bed \
         --out unique_CONTROL \
         --allow-no-sex

printf "Finding Positional and Allelic Duplicates\n"
printf ".\n.\n.\n"
plink1.9 --bfile unique_CASE \
         --list-duplicate-vars ids-only suppress-first \
         --out Dups2Remove_CASE \
         --allow-no-sex
plink1.9 --bfile unique_CONTROL \
         --list-duplicate-vars ids-only suppress-first \
         --out Dups2Remove_CONTROL \
         --allow-no-sex

printf "Removing Positional and Allelic Duplicates if they exist\n"
if [[ $(wc -l <Dups2Remove_CASE.dupvar) > 0 ]]
then 
    printf ".\n.\n.\n"
    plink1.9 --bfile unique_CASE \
             --exclude Dups2Remove_CASE.dupvar \
		     --keep-allele-order \
             --make-bed \
             --out CASE_Dups2Remove \
             --allow-no-sex
    printf "Done \n"
    mv CASE_Dups2Remove.bed unique_CASE.bed
    mv CASE_Dups2Remove.bim unique_CASE.bim
    mv CASE_Dups2Remove.fam unique_CASE.fam
fi

if [[ $(wc -l <Dups2Remove_CONTROL.dupvar) > 0 ]]
then 
    printf ".\n.\n.\n"
    plink1.9 --bfile unique_CONTROL \
             --exclude Dups2Remove_CONTROL.dupvar \
		     --keep-allele-order \
             --make-bed \
             --out CONTROL_Dups2Remove \
             --allow-no-sex
    printf "Done \n"
    mv CONTROL_Dups2Remove.bed unique_CONTROL.bed
    mv CONTROL_Dups2Remove.bim unique_CONTROL.bim
    mv CONTROL_Dups2Remove.fam unique_CONTROL.fam
fi

echo "==========================================================="
printf "STEP 2: Bim file modification (SNPs ID)\n"
echo "==========================================================="
# for case
mv unique_CASE.bim CASE_Step2.bim
mv unique_CASE.bed CASE_Step2.bed
mv unique_CASE.fam CASE_Step2.fam

# for control 
mv unique_CONTROL.bim CONTROL_Step2.bim
mv unique_CONTROL.bed CONTROL_Step2.bed
mv unique_CONTROL.fam CONTROL_Step2.fam

echo "==========================================================="
printf "STEP 3: extract common snps\n"
echo "==========================================================="
grep -f CASE_Step2.bim CONTROL_Step2.bim | cut -f2 > common_snps.txt

# for case
plink1.9 --bfile CASE_Step2 \
         --extract common_snps.txt \
		 --keep-allele-order \
         --make-bed \
         --out CASE_Step3 \
         --allow-no-sex

# for control
plink1.9 --bfile CONTROL_Step2 \
         --extract common_snps.txt \
		 --keep-allele-order \
         --make-bed \
         --out CONTROL_Step3 \
         --allow-no-sex

# merge bfile (keep executing when "exclude 3+ alleles" error occurs)
plink1.9 --bfile CASE_Step3 --bmerge CONTROL_Step3 --make-bed --out CASE_CONTROL_Step5 --allow-no-sex --keep-allele-order || true

echo "==========================================================="
printf "STEP 4: exclude 3+ alleles\n"
echo "==========================================================="
if [ -s CASE_CONTROL_Step5-merge.missnp ]
then
    # for case
    plink1.9 --bfile CASE_Step3 \
             --exclude CASE_CONTROL_Step5-merge.missnp \
		     --keep-allele-order \
             --make-bed \
             --out CASE_Step4 \
             --allow-no-sex

    # for control
    plink1.9 --bfile CONTROL_Step3 \
             --exclude CASE_CONTROL_Step5-merge.missnp \
		     --keep-allele-order \
             --make-bed \
             --out CONTROL_Step4 \
             --allow-no-sex

    # merge bfile
    plink1.9 --bfile CASE_Step4 \
             --bmerge CONTROL_Step4 \
		     --keep-allele-order \
             --make-bed \
             --out CASE_CONTROL_Step5 \
             --allow-no-sex
fi

echo "==========================================================="
printf "STEP 5: exclude INDEL (insertion deletion)\n"
echo "==========================================================="
# extract SNPs with INDEL
cat CASE_CONTROL_Step5.bim | awk '{if (length($5)>1 || length($6)>1) print $2}' > INDEL.txt

plink1.9 --bfile CASE_CONTROL_Step5 \
         --exclude INDEL.txt \
		 --keep-allele-order \
         --make-bed \
         --out $OUT_DIR/$OUT_FILE \
         --allow-no-sex


