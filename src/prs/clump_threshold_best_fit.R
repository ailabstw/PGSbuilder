#!/usr/bin/Rscript

# libraries
## install.packages(c("data.table","rms"), dependencies=TRUE)
library(optparse)
library(data.table)
library(rms)

# arguments
option_list <- list(
    make_option(c('-i', '--input'), help='the input bfile prefix'),
    make_option(c('-c', '--covar'), default=NULL, help='the covariate file'),
    make_option(c('-p', '--pfile'), help='the basename of prediction files'),
    make_option(c('-m', '--method'), help='classificaion(clf) or regression(reg) [default %default]', default='clf'),
    make_option(c('-d', '--dir'), help='the output directory'),
    make_option(c('-o', '--output'), help='the output basename')
)
arg <- parse_args(OptionParser(option_list=option_list)) # load arguments
p.threshold <- c('1', '1e-1', '1e-2', '1e-3', '1e-4', '1e-5', '1e-6', '1e-7', '1e-8')

# Read in the phenotype file 
phenotype <- read.table(paste0(arg$input, '.fam'), header=F)
colnames(phenotype) <- c('FID', 'IID', 'FATHER', 'MOTHER', 'SEX', 'PHENOTYPE')
phenotype <- subset(phenotype, select=c('FID', 'IID', 'PHENOTYPE'))
phenotype$PHENOTYPE[phenotype$PHENOTYPE==1] <- 0
phenotype$PHENOTYPE[phenotype$PHENOTYPE==2] <- 1

# add covariates and calculate the null R-square
if (!is.null(arg$covar)){
    cov <- read.table(arg$covar, header=T) # load covariates
    merge.pheno <- merge(phenotype, cov, by=c('FID', 'IID'), all=F)

    # We can then calculate the null model (model with PRS) using a linear/logistic regression
    if (arg$method == 'clf'){
        # for classification
        null.model <- lrm(PHENOTYPE~., data=merge.pheno[!colnames(merge.pheno)%in%c('FID','IID')])
        null.r2 <- null.model$stats['R2']
    } else {
        # for regression
        null.model <- lm(PHENOTYPE~., data=merge.pheno[!colnames(merge.pheno)%in%c('FID','IID')])
        null.r2 <- summary(null.model)$r.squared
    }

} else {
    merge.pheno <- phenotype
    null.r2 <- 0
}

# PRS
prs.result <- NULL
for (i in p.threshold){
    # check file existence
    filename <- paste0(arg$pfile, '.', i, '.profile')
    if (file.exists(filename) == FALSE){
        next
    }

    # load PRS of each threshold
    prs <- read.table(filename, header=T)
    merge.prs <- merge(merge.pheno, prs[c("FID","IID", "SCORESUM")], by=c('FID', 'IID'), all=F)
    # min-max normalization
    merge.prs$SCORESUM <- (merge.prs$SCORESUM - min(merge.prs$SCORESUM)) / (max(merge.prs$SCORESUM) - min(merge.prs$SCORESUM))
    if (arg$method == 'clf'){
        # for classification
        model <- lrm(PHENOTYPE~., data=merge.prs[!colnames(merge.prs)%in%c('FID','IID')], maxit=100)
        model.r2 <- model$stats['R2']
    } else {
        # for regression
        model <- lm(PHENOTYPE~., data=merge.prs[!colnames(merge.prs)%in%c('FID','IID')], maxit=100)
        model.r2 <- summary(model)$r.squared
    }
    prs.r2 <- model.r2 - null.r2
    prs.result <- rbind(prs.result, data.frame(Threshold=i, R2=prs.r2))
}
write.table(prs.result, file=paste0(arg$dir, '/', arg$output, '.r2.summary'), sep='\t', row.names=F, quote=F)

# best result
best.threshold <- prs.result[which.max(prs.result$R2),]$Threshold
print(paste0("Best fit threshold: ", best.threshold))
write(paste0(best.threshold, " 0 ", best.threshold), file=paste0(arg$dir, '/best_pvalue_range'))
