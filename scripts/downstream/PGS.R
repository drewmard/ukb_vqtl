library(data.table)

SNP <- 'rs887468'
ENV <- 'age'

s <- '80'
# s <- '20'
f.res <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/results/ukbb.gxe.',s,'.txt')
results <- fread(f.res,data.table = F,stringsAsFactors = F)
results[order(results$P_GxE),][1,]
subset(results,vQTL==SNP & E==ENV)

s <- '20'
f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/results/full_data_gxe.',s,'.txt')
df <- fread(f,data.table = F,stringsAsFactors = F)


BETA.SNP <- -3.98e-2
BETA.ENV <- 4.417e-2
score <- BETA.ENV*df[,ENV] + 
  BETA.SNP*df[,SNP]
cor(score,df$lymphocyte.count.na.rint,use='p')

# BETA.SNP.gxe <- -0.144
# BETA.ENV.gxe <- 0.00673
# BETA.GxE.gxe <- 0.00180
BETA.SNP.gxe <- -1.371e-1
BETA.ENV.gxe <- 4.315e-2
BETA.GxE.gxe <- 1.715e-3

score.gxe <- BETA.GxE.gxe*df[,SNP]*df[,ENV] + 
  BETA.ENV.gxe*df[,ENV] + 
  BETA.SNP.gxe*df[,SNP]
cor(score.gxe,df$lymphocyte.count.na.rint,use='p')


