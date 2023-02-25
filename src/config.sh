MEMORY=20000 #MB
THREAD=8

##### Path
HAPMAP_REF=
LIFTOVER_REF_DIR=
LDPRED_REF_DIR=
PRSCS_SRC=
PRSCS_REF_DIR=
GENEPI_REF_DIR=

##### Split option
TEST_RATIO=0.1
TARGET_RATIO=0.25

##### QC option
GENOME="hg38" # hg19 or hg38
SEP_FLAG="false" # whether filter QC for case and control seperately

# Variant filter
FLIPSCAN="true" # only case-control study
MAF_QC=0.05
GENO_QC=0.02
HWE_QC=0 # 1e-6 for quantitative phenotypes
STRICT_HWE="true" # only case-control study
HWE_MISSING="false"
HWE_QC_CONTROL=1e-6 # only case-control study
HWE_QC_CASE=1e-10 # only case-control study

# Individual filter
MIND_QC=0.02
HETER_SD=3
KINGCUTOFF=0.0442
CHECK_SEX="false"

# LD prune
PRUNE_WINDOW=50
PRUNE_STEP=5
PRUNE_THRESHOLD=0.2

# For population selection: CEU, ASW, MKK, MEX, CHD, CHB, JPT, LWK, TSI, GIH, YRI
POPULATION="CHD-CHB-JPT" # use "-" to separate the population
POP_SD=3 # the threshold of standard deviation 

# Hapmap ref
if [ "$GENOME" = "hg19" ]; then
    HM_FILENAME="${HAPMAP_REF}/hapmap3_r1_b37_fwd_consensus.hapmap3"
else
    HM_FILENAME="${HAPMAP_REF}/hapmap3_hg38.hapmap3"
fi
HM_INFO="${HAPMAP_REF}/relationships_w_pops_121708.hapmap3"


###### GWAS option
# pvalue threshold
THRESHOLD=1e-8

# covariate
PCA_FLAG="true"
PCA_COUNT=10


###### PRS option
POPULATION_PRS="ASN" # population: ASN, EUR, AFR
TOOLS="CandT,Lassosum,LDpred2" # CandT,PRSice2,Lassosum,LDpred2,PRScs,GenEpi