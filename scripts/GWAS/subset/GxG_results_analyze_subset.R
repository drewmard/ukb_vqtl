library(data.table)
# phenoName <- 'lymphocyte.count.rint.ALL'
phenoName <- 'monocyte.count.rint.ALL'
# phenoName <- 'neutrophil.count.rint.ALL'
# phenoName <- 'wbc.leukocyte.count.rint.ALL'
# phenoName <- 'rbc.erythrocyte.count.rint.ALL'

s <- '80'
df.80 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_GxG/ukbb.',phenoName,'.ALL.sub.GxG.',s,'.epi.qt'),data.table=F,stringsAsFactors = F)
colnames(df.80)[5:7] <- paste0(colnames(df.80)[5:7],'.',s)

s <- '20'
df.20 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_GxG/ukbb.',phenoName,'.ALL.sub.GxG.',s,'.epi.qt'),data.table=F,stringsAsFactors = F)
colnames(df.20)[5:7] <- paste0(colnames(df.20)[5:7],'.',s)

df <- merge(df.80,df.20,by=c('CHR1','SNP1','CHR2','SNP2'))
df.sub <- df
df.sub[order(df.sub$P.80)[1:5],] 

df.sub <- subset(df,CHR1!=CHR2)
df.sub[order(df.sub$P.80)[1:5],] 

df.sub <- subset(df,P.80 < 1e-3 & P.20 < 0.05 & CHR1 != CHR2 & sign(BETA_INT.80) == sign(BETA_INT.20))
df.sub[order(df.sub$P.80)[1:5],] 

df.sub <- subset(df,P.80 < 1e-3 & P.20 < 0.05 & CHR1 == CHR2 & sign(BETA_INT.80) == sign(BETA_INT.20))
df.sub[order(df.sub$P.80)[1:5],] 


df.sub <- df
# df.sub <- subset(df,SNP2=='rs146125856' & CHR1!=6)
# df.sub <- subset(df,CHR1!=CHR2)
# df.sub <- subset(df,CHR1!=CHR2 & SNP2!='rs146125856')
# df.sub <- subset(df,CHR1!=CHR2 & CHR1!=6 & CHR2!=6)
df.sub[order(df.sub$P)[1:5],] 

df.sub <- subset(df,CHR1!=CHR2)
df.sub[order(df.sub$P)[1:5],] 

df.sub <- subset(df,CHR1==15 | CHR2==15)
