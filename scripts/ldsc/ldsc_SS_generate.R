library(data.table)

DIR <- '/athena/elementolab/scratch/anm2868/vQTL/'

pheno='bmi'

mean <- fread(paste0(DIR,'ukb_vqtl/output/imputed/results/ukbb.bmi.ALL.results.txt'),data.table = F,stringsAsFactors = F)
var.raw <- fread(paste0(DIR,'ukb_vqtl/output/vGWAS_subset/ukbb.bmi.ALL.vGWAS.txt'),data.table = F,stringsAsFactors = F)
var.rint <- fread(paste0(DIR,'ukb_vqtl/output/vGWAS_subset/ukbb.bmi.rint.ALL.vGWAS.txt'),data.table = F,stringsAsFactors = F)
HLMM <- fread(paste0(DIR,'ukb_vqtl/output/GWAS/HLMM_results/ukbb.bmi.rint.ALL.HLMM.dispersion_nochr6.txt'),data.table = F,stringsAsFactors = F)

colnames(mean)[1] <- 'rs'
mean <- mean[,c('rs','CHR','BP','A1','A2','MAF','NMISS','BETA','SE','T','P')]
colnames(mean)[(ncol(mean)-3):ncol(mean)] <- paste0(colnames(mean)[(ncol(mean)-3):ncol(mean)],'.MEAN')
var.raw <- var.raw[,c('rs','BETA','SE','T','P')]; colnames(var.raw)[2:5] <- paste0(colnames(var.raw)[2:5],'.VAR.RAW')
var.rint <- var.rint[,c('rs','BETA','SE','T','P')]; colnames(var.rint)[2:5] <- paste0(colnames(var.rint)[2:5],'.VAR.RINT')
HLMM.sub <- HLMM[,c('SNP','dispersion','dispersion_se','dispersion_t','dispersion_pval')]

df.mg <- merge(mean,var.raw,by='rs')
df.mg <- merge(df.mg,var.rint,by='rs')
df.mg <- merge(df.mg,HLMM.sub,by.x='rs',by.y='SNP')

mean.sumstat <- df.mg[,c('rs','CHR','BP','NMISS','A1','A2','MAF','BETA.MEAN','SE.MEAN','T.MEAN','P.MEAN')]
var.sumstat <- df.mg[,c('rs','CHR','BP','NMISS','A1','A2','MAF','BETA.VAR.RAW','SE.VAR.RAW','T.VAR.RAW','P.VAR.RAW')]
var.rint.sumstat <- df.mg[,c('rs','CHR','BP','NMISS','A1','A2','MAF','BETA.VAR.RINT','SE.VAR.RINT','T.VAR.RINT','P.VAR.RINT')]
disp.sumstat <- df.mg[,c('rs','CHR','BP','NMISS','A1','A2','MAF','dispersion','dispersion_se','dispersion_t','dispersion_pval')]

colnames(mean.sumstat) <- c('SNP','CHR','BP','N','A1','A2','MAF','BETA','SE','Z','P')
colnames(var.sumstat) <- c('SNP','CHR','BP','N','A1','A2','MAF','BETA','SE','Z','P')
colnames(var.rint.sumstat) <- c('SNP','CHR','BP','N','A1','A2','MAF','BETA','SE','Z','P')
colnames(disp.sumstat) <- c('SNP','CHR','BP','N','A1','A2','MAF','BETA','SE','Z','P')

fwrite(mean.sumstat,'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/ldsc/mean.sumSS.txt',quote = F,na = 'NA',sep = '\t',row.names = F,col.names = T)
fwrite(var.sumstat,'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/ldsc/var.sumSS.txt',quote = F,na = 'NA',sep = '\t',row.names = F,col.names = T)
fwrite(var.rint.sumstat,'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/ldsc/var_rint.sumSS.txt',quote = F,na = 'NA',sep = '\t',row.names = F,col.names = T)
fwrite(disp.sumstat,'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/ldsc/disp.sumSS.txt',quote = F,na = 'NA',sep = '\t',row.names = F,col.names = T)

