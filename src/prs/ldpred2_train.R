#!/usr/bin/Rscript

### prepare library
## install.packages(c("R.utils", "optparse", "magrittr", "data.table", "remotes", "dplyr"), dependencies=TRUE)
## library(remotes)
## remotes::install_github("https://github.com/privefl/bigsnpr.git")
library(R.utils)
library(optparse)
library(data.table)
library(magrittr)
library(dplyr)
library(parallel)
library(bigsnpr)


### snp_modifyBuild: for liftOver
make_executable <- function(exe) {
  Sys.chmod(exe, mode = (file.info(exe)$mode | "111"))
}

snp_modifyBuild <- function(info_snp, liftOver, db_dir, from = "hg18", to = "hg19") {

  if (!all(c("chr", "pos") %in% names(info_snp)))
    stop2("Please use proper names for variables in 'info_snp'. Expected %s.",
          "'chr' and 'pos'")

  # Make sure liftOver is executable
  liftOver <- normalizePath(liftOver)
  make_executable(liftOver)

  # Need BED UCSC file for liftOver
  BED <- tempfile(fileext = ".BED")
  info_BED <- with(info_snp, data.frame(
    paste0("chr", chr), pos0 = pos - 1L, pos, id = rows_along(info_snp)))
  bigreadr::fwrite2(info_BED, BED, col.names = FALSE, sep = " ")

  # Need chain file
  chain_name <- paste0(from, "To", tools::toTitleCase(to), ".over.chain.gz")
  if (file.exists(file.path(db_dir, chain_name))){
      chain <- file.path(db_dir, chain_name)
  } else {
      url <- paste0("ftp://hgdownload.cse.ucsc.edu/goldenPath/", from, "/liftOver/",
                    from, "To", tools::toTitleCase(to), ".over.chain.gz")
      chain <- tempfile(fileext = ".over.chain.gz")
      utils::download.file(url, destfile = chain)
  }

  # Run liftOver (usage: liftOver oldFile map.chain newFile unMapped)
  lifted <- tempfile(fileext = ".BED")
  unmapped <- tempfile(fileext = ".txt")
  system(paste(liftOver, BED, chain, lifted, unmapped))

  # readLines(lifted, n = 5)
  new_pos <- bigreadr::fread2(lifted)

  # readLines(unmapped, n = 6)
  bad <- grep("^#", readLines(unmapped), value = TRUE, invert = TRUE)
  print(paste0(length(bad), " variants have not been mapped."))

  info_snp$pos <- NA
  info_snp$pos[new_pos$V4] <- new_pos$V3
  info_snp
}


### arguments
option_list <- list(
    make_option(c('-i', '--input'), help='the input bfile prefix'),
    make_option(c('-a', '--assoc'), help='the summary statisitics file'),
    make_option(c('-d', '--out_dir'), help='the working directory'),
    make_option(c('-o', '--output'), help='the output basename'),
    make_option(c('-l', '--liftover_dir'), help='the directory of the liftOver database', default=''),
    make_option(c('-b', '--db_dir'), help='the directory of 1000 Genome map'),
    make_option(c('-m', '--hapmap'), help='use hapmap3 as SNP filter', action='store_true', default=F),
    make_option(c('-g', '--genome'), help='the reference genome: hg18, hg19, or hg38 [default %default]', default='hg19'),
    make_option(c('-t', '--method'), help='infinitesimal(1), grid-sparse(2), grid-no-sparse(3), or auto(4) [default %default]', default=3, type='integer')
)

arg <- parse_args(OptionParser(option_list=option_list)) # load arguments
ncores <- if (detectCores() > 8) 8 else detectCores() # ncores = min(8, available)


####################################
# Preprocess data
####################################
### data
# load summary statistics
ss <- fread(arg$assoc)
origin <- c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'A1', 'A2', 'A1_FREQ', 'OBS_CT', 'OR', 'BETA', 'BETA_SE', 'STAT', 'P', 'LOG10_P')
rename <- c('chr', 'pos', 'rsid', 'ref', 'alt', 'a1', 'a0', 'a1freq', 'n_eff', 'OR', 'beta', 'beta_se', 'stat', 'p', 'logp')
colnames(ss) <- dplyr::recode(colnames(ss), !!!setNames(rename, origin))
ss <- ss[!ss$p == 'null',]
ss$chr <- as.integer(ss$chr)

# load HapMap reference
hapmap <- readRDS(paste0(arg$db_dir, '/hapmap3.rds'))
if (arg$genome == 'hg19'){
    hapmap <- hapmap[c('chr', 'pos', 'rsid', 'a0', 'a1')]
} else {
    hapmap <- hapmap[c('chr', paste0('pos_', arg$genome), 'rsid', 'a0', 'a1')]
    names(hapmap) <- c('chr', 'pos', 'rsid', 'a0', 'a1')
}
hapmap$chr <- as.integer(hapmap$chr)

# load genotype (bfile)
if (!(file.exists(paste0(arg$out_dir, '/', arg$output, '.rds')))) {
    snp_readBed(paste0(arg$input, '.bed'), backingfile=paste0(arg$out_dir, '/', arg$output)) # build .rds file
}
obj.bigSNP <- snp_attach(paste0(arg$out_dir, '/', arg$output, '.rds'))
G <- obj.bigSNP$genotypes
y <- obj.bigSNP$fam$affection

# extract the SNP information from input
map <- obj.bigSNP$map[-3]
names(map) <- c("chr", "rsid", "pos", "a1", "a0")
map$chr <- as.integer(map$chr)


### SNP-matching dataset
# match SNPs between summary statistics, genotype, (and hapmap)
if (arg$hapmap){
    match_snp <- snp_match(ss, hapmap, join_by_pos=TRUE, match.min.prop=0.25)
    match_snp <- match_snp[c('chr', 'pos', 'rsid', 'a0', 'a1', 'n_eff', 'OR', 'beta_se', 'p', 'beta')]
    match_snp <- snp_match(match_snp, map, join_by_pos=TRUE, match.min.prop=0.25)
} else {
    match_snp <- snp_match(ss, map, join_by_pos=TRUE, match.min.prop=0.25)
}

# liftOver to hg19
if (arg$genome == 'hg38'){
    info_snp <- snp_modifyBuild(match_snp, 'liftOver', arg$liftover_dir, from='hg38', to='hg19')
} else if (arg$genome == 'hg18'){
    info_snp <- snp_modifyBuild(match_snp, 'liftOver', arg$liftover_dir, from='hg18', to='hg19')
} else {
    info_snp <- match_snp
}

# sort chr and pos of info_snp
ord <- with(info_snp, order(chr, pos))
info_snp_sort <- info_snp[ord,]


####################################
# Calculate the LD matrix
####################################
# open a temporary file
tmp <- tempfile(tmpdir = paste0(arg$out_dir, "/tmp-data"))
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

# get the CM information from 1000 Genome
POS <- snp_asGeneticPos(info_snp_sort$chr, info_snp_sort$pos, dir = arg$db_dir)

### calculate LD
corr <- NULL
ld <- NULL
cc = 1
for (chr in unique(info_snp_sort$chr)) {
    # extract target chromosome
    ind.chr <- which(info_snp_sort$chr == chr) # index of merged dataframe
    ind.chr2 <- info_snp_sort$`_NUM_ID_`[ind.chr] # index of the original genotype

    # calculate the LD
    corr0 <- snp_cor(
            G,
            ind.col = ind.chr2,
            ncores = ncores,
            infos.pos = POS[ind.chr],
            size = 3 / 1000
        )
    if (cc == 1) {
        df_beta <- info_snp_sort[ind.chr, c("beta", "beta_se", "n_eff", "_NUM_ID_")]
        ld <- Matrix::colSums(corr0^2)
        corr <- as_SFBM(corr0, tmp, compact=TRUE)
    } else {
        df_beta <- rbind(df_beta, info_snp_sort[ind.chr, c("beta", "beta_se", "n_eff", "_NUM_ID_")])
        ld <- c(ld, Matrix::colSums(corr0^2))
        corr$add_columns(corr0, nrow(corr))
    }
    cc = cc + 1
}


####################################
# Perform LD score regression and get heritability
####################################
ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2, 
                               sample_size = n_eff, blocks = NULL))
h2_est <- ldsc[["h2"]]
if (h2_est <= 0){
    h2_est <- 0.01
}


####################################
# Compute LDpred2 models
####################################
# fill missing genotype with rounded mean to 2 decimal places
G2 <- snp_fastImputeSimple(G, method='mean2')

### infinitesimal mode
if (arg$method == 1){
    print("execute infinitesimal mode...")
    best_beta <- snp_ldpred2_inf(corr, df_beta, h2_est)

### grid mode
} else if ((arg$method == 2) || (arg$method == 3)){
    # params
    if (arg$method == 2){
        sparse <- TRUE
        print("execute grid-sparse mode...")
    } else {
        sparse <- FALSE
        print("execute grid-no-sparse mode...")
    }
    h2_seq <- round(h2_est * c(0.3, 0.7, 1, 1.4), 4)
    p_seq <- signif(seq_log(1e-5, 1, length.out=21), 2)
    params <- expand.grid(p=p_seq, h2=h2_seq, sparse=sparse)

    # grid
    beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores=ncores)
    params$sparsity <- colMeans(beta_grid==0)

    # pred
    bigparallelr::set_blas_ncores(ncores)
    pred_grid <- big_prodMat(G2, beta_grid, ind.col=df_beta$`_NUM_ID_`)
    params$score <- apply(pred_grid, 2, cor, y=y)

    # best beta
    best_beta <- params %>%
        mutate(id = row_number()) %>%
        arrange(desc(score)) %>%
        slice(1) %>%
        pull(id) %>%
        beta_grid[, .]

### auto mode
} else {
    print("execute auto mode...")
    multi_auto <- snp_ldpred2_auto(
        corr,
        df_beta,
        h2_init = h2_est,
        vec_p_init = seq_log(1e-4, 0.5, length.out=30),
        allow_jump_sign = FALSE,
        shrink_corr = 0.95,
        ncores = ncores
    )
    beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)
    pred_auto <- big_prodMat(G2, beta_auto, ind.col=df_beta$`_NUM_ID_`)
    pred_scaled <- apply(pred_auto, 2, sd)
    best_beta <- rowMeans(beta_auto[,abs(pred_scaled - median(pred_scaled)) < 3 * mad(pred_scaled)])
}


####################################
# Compute Lassosum2 models
####################################
## lassosum2
#beta_lassosum2 <- snp_lassosum2(
#    corr,
#    df_beta,
#    ncores=ncores
#)
#params <- attr(beta_lassosum2, "grid_param")

## pred
#bigparallelr::set_blas_ncores(ncores)
#pred_lassosum2 <- big_prodMat(G2, beta_lassosum2, ind.col=df_beta$`_NUM_ID_`)
#params$score <- apply(pred_lassosum2, 2, cor, y=y)

## best beta
#best_lassosum2 <- params %>%
#  mutate(id = row_number()) %>%
#  arrange(desc(score)) %>%
#  slice(1) %>%
#  pull(id) %>%
#  beta_lassosum2[, .]


####################################
# save models
####################################
# save
saveRDS(best_beta, file=paste0(arg$out_dir, '/', arg$output, '.beta.rds'))
saveRDS(df_beta, file=paste0(arg$out_dir, '/', arg$output, '.snp.rds'))

# remove tmp files
unlink(paste0(arg$out_dir, '/tmp-data'), recursive=T)
