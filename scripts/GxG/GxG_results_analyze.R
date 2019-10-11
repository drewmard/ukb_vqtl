library(data.table)
phenoName <- 'lymphocyte.count.rint.ALL'
# phenoName <- 'monocyte.count.rint.ALL'
df <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG/ukbb.',phenoName,'.ALL.sub.GxG.epi.qt'),data.table=F,stringsAsFactors = F)

df.sub <- df
# df.sub <- subset(df,SNP2=='rs146125856' & CHR1!=6)
df.sub <- subset(df,CHR1!=CHR2)
# df.sub <- subset(df,CHR1!=CHR2 & SNP2!='rs146125856')
# df.sub <- subset(df,CHR1!=CHR2 & CHR1!=6 & CHR2!=6)
df.sub[order(df.sub$P)[1:5],] 

