library(data.table)
pheno <- 'lymphocyte.count'

f1 <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_subset/ukbb.',paste0(pheno,'.ALL'),'.vGWAS.txt')
f2 <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_subset/ukbb.',paste0(pheno,'.rint.ALL'),'.vGWAS.txt')
f3 <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/imputed/results/ukbb.',paste0(pheno,'.ALL'),'.results.txt')

var.raw <- fread(f1,data.table = F,stringsAsFactors = F)
var.rint <- fread(f2,data.table = F,stringsAsFactors = F)
mean.raw <- fread(f3,data.table = F,stringsAsFactors = F); colnames(mean.raw)[1] <- 'rs'

f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.',pheno,'.merged_subset.bim')
bim <- fread(f,data.table = F,stringsAsFactors = F)
var.raw.sub <- subset(var.raw,rs %in% bim[,2])
var.rint.sub <- subset(var.rint,rs %in% bim[,2])
mean.raw.sub <- subset(mean.raw,rs %in% bim[,2])

mean.raw.sub <- mean.raw.sub[,c('rs','CHR','BP','A1','A2','MAF','BETA','P')]
colnames(mean.raw.sub)[(ncol(mean.raw.sub)-1):ncol(mean.raw.sub)] <- paste0(colnames(mean.raw.sub)[(ncol(mean.raw.sub)-1):ncol(mean.raw.sub)],'.MEAN')
var.raw.sub <- var.raw.sub[,c('rs','BETA','P')]; colnames(var.raw.sub)[2:3] <- paste0(colnames(var.raw.sub)[2:3],'.VAR.RAW')
var.rint.sub <- var.rint.sub[,c('rs','BETA','P')]; colnames(var.rint.sub)[2:3] <- paste0(colnames(var.rint.sub)[2:3],'.VAR.RINT')

df.mg <- merge(merge(mean.raw.sub,var.raw.sub,by='rs'),var.rint.sub,by='rs')

#################################################################

f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.',pheno,'.merged_subset.LD.ld')
ld <- fread(f,data.table = F,stringsAsFactors = F)
ld <- subset(ld,R2 > 0.1)
df.mg.ld <- subset(df.mg, rs %in% c(ld$SNP_A,ld$SNP_B))
df.mg.ld <- df.mg.ld[order(df.mg.ld$CHR),]

df.mg2 <- df.mg[,c('rs','CHR','BP','P.MEAN','P.VAR.RAW','P.VAR.RINT')]
