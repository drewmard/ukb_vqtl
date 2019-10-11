library(data.table)

SNP1 <- 'rs2242659'
SNP2 <- 'rs146125856'
phenoName <- 'lymphocyte.count.rint.ALL'

SNP1 <- 'rs887468'
SNP2 <- 'rs2516491'
phenoName <- 'lymphocyte.count.rint.ALL'



############
library(data.table)
library(stringr)
phenoName <- 'lymphocyte.count.rint.ALL'
x <- strsplit(phenoName,'\\.')[[1]]; phenoName2 <- paste0(c(x[-length(x)],'na'),collapse='.')

s <- '80'
f.res <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG/ukbb.lymphocyte.count.rint.ALL.ALL.sub.GxG.epi.qt')
results <- fread(f.res,data.table = F,stringsAsFactors = F)
res.sub <- subset(results,P < 0.05/nrow(results))
# res.sub <- res.sub[which(!duplicated(res.sub$SNP1)),]

res.sub <- data.table(res.sub)
res.sub <- res.sub[, .SD[which.min(P)],by=SNP1]
res.sub <- res.sub[, .SD[which.min(P)],by=SNP2]

s <- '20'
f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/results/full_data_gxe.',s,'.txt')
df <- fread(f,data.table = F,stringsAsFactors = F)

score.gxg <- 0
for (i in 1:nrow(res.sub)) {
# for (i in 1:5) {
  BETA.GXG <- res.sub$BETA_INT[i]
  SNP1 <- str_replace_all(res.sub$SNP1[i],'-','.')
  SNP2 <- str_replace_all(res.sub$SNP2[i],'-','.')
  score.gxg <- score.gxg + BETA.GXG*df[,SNP1]*df[,SNP2]
}
cor(score.gxg,df[,phenoName2],use='p')




BETA.SNP1.GXG <- BETA.SNP1
BETA.SNP2.GXG <- BETA.SNP2
# BETA.GXG <- -0.0773763
BETA.GXG <- 0.0658181

score.gxg <- BETA.SNP1*df[,SNP1] + 
  BETA.SNP2*df[,SNP2] +
  BETA.GXG*df[,SNP1]*df[,SNP2]





BETA.SNP1 <- subset(results,SNP==SNP1)$BETA.x
BETA.SNP2 <- subset(results,SNP==SNP2)$BETA.x
score <- BETA.SNP1*df[,SNP1] + 
  BETA.SNP2*df[,SNP2]
cor(score,df[,phenoName2],use='p')

BETA.SNP1.GXG <- BETA.SNP1
BETA.SNP2.GXG <- BETA.SNP2
# BETA.GXG <- -0.0773763
BETA.GXG <- 0.0658181

score.gxg <- BETA.SNP1*df[,SNP1] + 
  BETA.SNP2*df[,SNP2] +
  BETA.GXG*df[,SNP1]*df[,SNP2]
  
cor(score.gxg,df[,phenoName2],use='p')


