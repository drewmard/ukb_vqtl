library(data.table)
args = commandArgs(trailingOnly=TRUE)
phenotype=args[1]
phenotype='lymphocyte.count.rint.ALL'

library(data.table)
phenotype='bmi.rint.ALL'

f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/HLMM_results/ukbb.',phenotype,'.HLMM.txt')
results <- fread(f.out,data.table = F,stringsAsFactors = F)
# results <- fread(f.out,data.table = F,stringsAsFactors = F,fill=TRUE)

f.mean <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/imputed/results/ukbb.lymphocyte.count.rint.ALL.results.txt'
df.mean <- fread(f.mean,data.table = F,stringsAsFactors = F)
chr6 <- subset(df.mean,CHR==6)$SNP

results <- subset(results,frequency > 0.05 & !(SNP%in%chr6))
source('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/bin/estimate_dispersion_effects.R')
cor(results[,c('add','var','dispersion')])
cor(results[,c('add_pval','var_pval','av_pval','dispersion_pval')])

results <- results[order(results$dispersion_pval,decreasing = T),]; head(results)
results <- results[order(results$av_pval,decreasing = T),]; head(results)
df.sub <- subset(results,av_pval < 1e-5 & dispersion_pval > 1e-3)

f.var <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_subset/','ukbb.',phenotype,'.vGWAS.txt')
df.var <- fread(f.var,data.table = F,stringsAsFactors = F)

results.mg <- merge(results,df.var[,c('rs','BETA','P')],by.x='SNP',by.y='rs')

results <- subset(results,frequency > 0.05 & )
cor(results.mg[,c('add','var','dispersion','BETA')])
cor(results.mg[,c('add_pval','var_pval','dispersion_pval','av_pval','P')])
df.sub <- df.sub[order(df.sub$var_pval,decreasing = T),]; head(df.sub)
df.sub <- df.sub[order(df.sub$dispersion_pval,decreasing = T),]; head(df.sub)


f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_subset/','ukbb.',phenotype,'.vGWAS.txt')
print(paste0('Writing: ',f.out))
fwrite(df.res.save,f.out,sep='\t',quote=F,col.names = T,row.names = F,na="NA")

