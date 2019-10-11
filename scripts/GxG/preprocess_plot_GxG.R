library(data.table)

s <- '80'
f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/results/full_data_gxe.',s,'.txt')
df <- fread(f,data.table = F,stringsAsFactors = F)

SNP1 <- 'rs887468'
SNP2 <- 'rs2516491'
phenoName <- 'lymphocyte.count.rint.ALL'

table(df[,c(SNP1,SNP2)])
cor(df[,c(SNP1,SNP2)])
# df.sub <- (df[,c(SNP1,SNP2,'resid1')])
# 
# x <-paste0('resid1~',SNP1,'*',SNP2)
df.sub <- (df[,c(SNP1,SNP2,'lymphocyte.count.rint.na')])

x <-paste0('lymphocyte.count.rint.na','~',SNP1,'*',SNP2)

x <- formula(x)
mod <- lm(x,data=df.sub)
summary(mod)$coef

x <-paste0('lymphocyte.count.rint.na','~',SNP1)
x <- formula(x)
mod <- lm(x,data=df.sub)
summary(mod)$coef

x <-paste0('lymphocyte.count.rint.na','~',SNP2)
x <- formula(x)
mod <- lm(x,data=df.sub)
summary(mod)$coef

x <-paste0('lymphocyte.count.rint.na','~',SNP1)
x <- formula(x)
mod <- lm(x,data=subset(df.sub,df.sub[,SNP2]==0))
summary(mod)$coef
mod <- lm(x,data=subset(df.sub,df.sub[,SNP2]==1))
summary(mod)$coef
mod <- lm(x,data=subset(df.sub,df.sub[,SNP2]==2))
summary(mod)$coef

x <-paste0('lymphocyte.count.rint.na','~',SNP2)
x <- formula(x)
mod <- lm(x,data=subset(df.sub,df.sub[,SNP1]==0))
summary(mod)$coef
mod <- lm(x,data=subset(df.sub,df.sub[,SNP1]==1))
summary(mod)$coef
mod <- lm(x,data=subset(df.sub,df.sub[,SNP1]==2))
summary(mod)$coef

fwrite(df.sub,'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG/plot/GxG_data.txt',row.names = F,na = 'NA',quote = F,col.names = T,sep = '\t')






