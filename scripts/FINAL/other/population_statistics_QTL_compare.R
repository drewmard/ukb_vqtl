library(data.table)
f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/andrew_copies/maf/ukbb.','ALL','.impute.frqx')
df.frqx <- fread(f,data.table = F,stringsAsFactors = F)

f <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/sig_results/bmi.sig.txt'
df <- fread(f,data.table = F,stringsAsFactors = F,header = T)

df.sub <- subset(df.frqx,SNP %in% df$rs)
df.sub <- merge(df,df.frqx[,c('SNP','hom_a1','MISS')],by.x='rs',by.y='SNP')

t.test(
  subset(df.sub,Mean.QTL==1)$hom_a1,
  subset(df.sub,Raw.vQTL==1)$hom_a1
)
t.test(
  subset(df.sub,Mean.QTL==1)$hom_a1,
  subset(df.sub,Rint.vQTL==1)$hom_a1
)
t.test(
  subset(df.sub,Mean.QTL==1)$hom_a1,
  subset(df.sub,dQTL==1)$hom_a1
)

#

t.test(
  subset(df.sub,Mean.QTL==1)$MISS,
  subset(df.sub,Raw.vQTL==1)$MISS
)
t.test(
  subset(df.sub,Mean.QTL==1)$MISS,
  subset(df.sub,Rint.vQTL==1)$MISS
)
t.test(
  subset(df.sub,Mean.QTL==1)$MISS,
  subset(df.sub,dQTL==1)$MISS
)

#

t.test(
  subset(df.sub,Mean.QTL==1)$MAF,
  subset(df.sub,Raw.vQTL==1)$MAF
)
t.test(
  subset(df.sub,Mean.QTL==1)$MAF,
  subset(df.sub,Rint.vQTL==1)$MAF
)
t.test(
  subset(df.sub,Mean.QTL==1)$MAF,
  subset(df.sub,dQTL==1)$MAF
)


# 
# f <- '/athena/elementolab/scratch/anm2868/open_targets/query/bmi.mean.txt'
# df <- fread(f,data.table = F,stringsAsFactors = F,header = F)
# 
# mean.raw <- subset(df.frqx,SNP%in%df[,1])
# mean.raw <- subset(mean.raw,!(SNP%in%var.raw[,1]))
# t.test(mean.raw[,2],var.raw[,2])
# t.test(mean.raw[,3],var.raw[,3])
# t.test(mean.raw[,4],var.raw[,4])
