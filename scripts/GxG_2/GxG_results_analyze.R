library(data.table)
phenoName <- 'lymphocyte.count'
# phenoName <- 'monocyte.count.rint.ALL'
# phenoName <- 'neutrophil.count.rint.ALL'
# phenoName <- 'wbc.leukocyte.count.rint.ALL'
# phenoName <- 'rbc.erythrocyte.count.rint.ALL'

s <- '80'
f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.',phenoName,'.ALL.sub.GxG.',s,'.epi.qt')
# f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_GxG/ukbb.',phenoName,'.ALL.sub.GxG.epi.qt')
df.80 <- fread(f,data.table=F,stringsAsFactors = F)
colnames(df.80)[5:7] <- paste0(colnames(df.80)[5:7],'.',s)
df.80[order(df.80$P.80)[1:5],] 
# unique(subset(df.80,P.80 < 0.05/nrow(df.80)/5)[,'CHR1'])
# monocyte: 4
# lymphocyte: 28
# neutrophil: 1

subset(df.80,P.80 < 5e-8)

s <- '20'
df.20 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.',phenoName,'.ALL.sub.GxG.',s,'.epi.qt'),data.table=F,stringsAsFactors = F)
colnames(df.20)[5:7] <- paste0(colnames(df.20)[5:7],'.',s)

df <- merge(df.80,df.20,by=c('CHR1','SNP1','CHR2','SNP2'))
df.sub <- df
df.sub[order(df.sub$P.80)[1:5],] 

f.freq <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.lymphocyte.count.ALL.sub.GxG.frq'
df.freq <- fread(f.freq,data.table = F,stringsAsFactors = F)
df2 <- merge(df,df.freq[,c('SNP','MAF')],by.x='SNP1',by.y='SNP')
df2 <- merge(df2,df.freq[,c('SNP','MAF')],by.x='SNP2',by.y='SNP')
colnames(df2)[which(colnames(df2) %in% c('MAF.x','MAF.y'))] <- c('MAF1','MAF2')

df.sub <- subset(df2,CHR1!=CHR2 & MAF1 > 0.3 & MAF2 > 0.1 & P.80 < 0.05)
mean(sign(df.sub$BETA_INT.80) == sign(df.sub$BETA_INT.20))

df.sub <- subset(df2,CHR1!=CHR2 & MAF1 > 0.2 & MAF2 > 0.2 & P.80 < 0.05 & (sign(BETA_INT.80) == sign(BETA_INT.20)))
df.sub[order(df.sub$P.80)[1:5],] 


df <- df.80
# df.sub <- subset(df,CHR1!=6)
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
