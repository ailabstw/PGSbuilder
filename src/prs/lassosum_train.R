#!/usr/bin/Rscript

## apt install build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev
## install.packages(c("devtools","RcppArmadillo", "data.table", "Matrix"), dependencies=TRUE)
## library(devtools)
## install_github("tshmak/lassosum")
#!/usr/bin/Rscript

### prepare library
## install.packages(c("optparse", "data.table", "parallel", "R.utils"), dependencies=TRUE)
library(optparse)
library(data.table)
library(lassosum)
library(parallel)

### arguments
option_list <- list(
    make_option(c('-i', '--input'), help='the input bfile prefix'),
    make_option(c('-a', '--assoc'), help='the summary statisitics file'),
    make_option(c('-l', '--ld'), default='EUR', help='population for LD: EUR, ASN, or AFR [default %default]'),
    make_option(c('-g', '--genome', default='hg19', help='the reference genome: hg19, hg38 [default %default]')),
    make_option(c('-d', '--dir'), help='the output directory'),
    make_option(c('-o', '--output'), help='the output basename')
)
arg <- parse_args(OptionParser(option_list=option_list)) # load arguments
threads <- if (detectCores() > 8) 8 else detectCores() # threads = min(8, available)
cl <- makeCluster(threads) # load threads
ld <- paste0(arg$ld, '.', arg$genome) # LD source


### summary statistics
# load summary statistics
ss <- fread(arg$assoc)
ss <- ss[!P == 0] # remove P-value = 0, which causes problem in the transformation

# Transform the P-values into correlation
cor <- p2cor(p=ss$P, n=ss$OBS_CT, sign=log(ss$OR), min.n=min(ss$OBS_CT))

# sample size
sample.num <- as.integer(system2('wc', args=c('-l', paste0(arg$input, '.fam'), " | awk '{print $1}'"), stdout=TRUE))
if (sample.num > 20000) {
    sample <- 5000
} else {
    sample <- NULL
}


### PRS
# run the lassosum pipeline
out <- lassosum.pipeline(
    cor = cor,
    chr = ss$"#CHROM",
    pos = ss$POS,
    A1 = ss$A1,
    ref.bfile = arg$input,
    sample = sample,
    LDblocks = ld, 
    cluster = cl
)

# train
res <- validate(out)


### save
# model
model <- subset(out, s=res$best.s, lambda=res$best.lambda)
saveRDS(model, file=paste0(arg$dir, '/', arg$output, '.model.rds'))

# beta
beta <- replicate(length(model$test.extract), 0)
beta[model$test.extract] <- model$beta[[1]] # model$beta[[model$s]]
write(beta, file=paste0(arg$dir, '/', arg$output, '.beta'), ncolumns=1)