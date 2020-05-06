# library(data.table)
# args = commandArgs(trailingOnly=TRUE)
# phenotype=args[1]
# phenotype='monocyte.count.rint.ALL'
# phenotype='lymphocyte.count.rint.ALL'

library(data.table)
phenotype='bmi.rint.ALL'

f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/HLMM_results/ukbb.',phenotype,'.HLMM.txt')
results <- fread(f.out,data.table = F,stringsAsFactors = F)

f.mean <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/imputed/results/ukbb.lymphocyte.count.rint.ALL.results.txt'
df.mean <- fread(f.mean,data.table = F,stringsAsFactors = F)
chr6 <- subset(df.mean,CHR==6)$SNP
MAF_0.05 <- subset(df.mean,MAF > 0.05)$SNP

results.nochr6 <- subset(results,frequency > 0.05 & !(SNP%in%chr6))
results.maf <- subset(results,SNP %in% MAF_0.05)

# sum(!(results$SNP%in%df.mean$SNP))
# sum(!(df.mean$SNP%in%results$SNP))
# length(MAF_0.05);nrow(results.maf)

##
# source('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/bin/estimate_dispersion_effects.R')

require(MASS)
r_var_mean=rlm(var~0+add,data=results.nochr6)

mean_noise=mean(results.nochr6$add_se^2,na.rm=T)
noise_adjustment=1+mean_noise/(var(results.nochr6$add,na.rm=T)-mean_noise)

r_av=r_var_mean$coefficients[1]*noise_adjustment

results.maf$dispersion=results.maf$var-r_av*results.maf$add
results.maf$dispersion_se=sqrt(results.maf$var_se^2+(r_av^2)*results.maf$add_se^2)
results.maf$dispersion_t=results.maf$dispersion/results.maf$dispersion_se
results.maf$dispersion_pval=-log10(pchisq(results.maf$dispersion_t^2,1,lower.tail=F))

###

f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/HLMM_results/ukbb.',phenotype,'.HLMM.dispersion_nochr6.txt')
fwrite(results,f.out,row.names = F,col.names = T,quote = F,sep = '\t',na = 'NA')

# /athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/HLMM_results/ukbb.lymphocyte.count.rint.ALL.HLMM.dispersion_nochr6.txt
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################


# cor(results[,c('add','var','dispersion')],use='p')
# cor(results[,c('add_pval','var_pval','av_pval','dispersion_pval')])
# x <- subset(results,av_pval > -log10(5e-8) & dispersion_pval > -log10(1e-3))
# dim(subset(results,av_pval > -log10(5e-8) & dispersion_pval > -log10(1e-3)))
# dim(subset(results,var_pval < -log10(5e-8) & dispersion_pval > -log10(1e-5)))
# dim(subset(results,add_pval < -log10(5e-8) & dispersion_pval > -log10(1e-5)))
