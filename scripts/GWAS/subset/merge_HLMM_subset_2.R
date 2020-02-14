library(data.table)
args = commandArgs(trailingOnly=TRUE)
phenotype=args[1]
# phenotype='lymphocyte.count.rint.ALL'

f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/HLMM_results/ukbb.',phenotype,'.HLMM.txt')
results <- fread(f.out,data.table = F,stringsAsFactors = F)

f.var <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_subset/','ukbb.',phenotype,'.vGWAS.txt')
df.var <- fread(f.var,data.table = F,stringsAsFactors = F)

results.mg <- merge(results,df.var[,c('rs','BETA','P')],by.x='SNP',by.y='rs')

df.sub <- subset(results.mg,frequency > 0.05)
cor(results.mg[,c('add','var','dispersion','BETA')])
cor(results.mg[,c('add_pval','var_pval','dispersion_pval','av_pval','P')])
df.sub <- df.sub[order(df.sub$var_pval,decreasing = T),]; head(df.sub)
df.sub <- df.sub[order(df.sub$dispersion_pval,decreasing = T),]; head(df.sub)


f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_subset/','ukbb.',phenotype,'.vGWAS.txt')
print(paste0('Writing: ',f.out))
fwrite(df.res.save,f.out,sep='\t',quote=F,col.names = T,row.names = F,na="NA")

