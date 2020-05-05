library(data.table)
df <- fread('/home/kulmsc/athena/ukbiobank/qc/ukb_sample_qc.txt',data.table = F,stringsAsFactors = F)
df2 <- df[,c(4,5,6,19,20,23,24,3,26:45)] # 19=exc heterozyg, 20 = aneuploidy, 23 = exc rel,24 = brits,3= genotyping array
colnames(df2)[1:8] <- c(
  'batch',
  'plate',
  'well',
  'het.missing.outliers',
  'putative.sex.chromosome.aneuploidy',
  'excess.relatives',
  'in.white.British.ancestry.subset',
  'genotyping.array')
colnames(df2)[9:ncol(df2)] <- paste0('PC',1:20)
i.remove <- which(df2$excess.relatives==1 |
                    df2$putative.sex.chromosome.aneuploidy==1 |
                    df2$het.missing.outliers==1)
df2$QC_In <- 1; df2$QC_In[i.remove] <- 0

fwrite(df2,'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/ukb_sample_qc.txt',na='NA',row.names = F,col.names = T,quote = F,sep = '\t')

