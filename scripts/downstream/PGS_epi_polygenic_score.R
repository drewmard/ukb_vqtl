library(data.table)
library(stringr)

phenoName2 <- 'lymphocyte.count.na'

for (i in 1:22) {
  print(i)
  # f <- paste0('/home/kulmsc/athena/forAndrew/clumping/ukbb.',i,'.profile')
  f <- paste0('/home/kulmsc/athena/forAndrew/noCorrection/ukbb.',i,'.profile')
  df <- fread(f,data.table = F,stringsAsFactors = F)
  if (i == 1) {
    score.save <- df
  } else {
    score.save$SCORE <-score.save$SCORE+df$SCORE
  }
}

s <- '20'
# s <- '80'
f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/results/full_data_gxe.',s,'.txt')
df <- fread(f,data.table = F,stringsAsFactors = F)

cor(df[,c('rs2242659','rs146125856')],use='p')

s <- '80'
f.res <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG/ukbb.lymphocyte.count.rint.ALL.ALL.sub.GxG.epi.qt')
results <- fread(f.res,data.table = F,stringsAsFactors = F)
res.sub <- subset(results,P < 0.05/nrow(results))
res.sub <- res.sub[order(res.sub$P),]
head(subset(res.sub,CHR1!=CHR2))
subset(res.sub,SNP1=='rs2242659' & SNP2=='rs146125856')
# res.sub <- subset(res.sub,!(CHR1==CHR2 & (SNP1=='rs2242659') | (SNP2=='rs2242659')))

SNP_list <- c()
for (i in 1:nrow(res.sub)) {
# for (i in 1:184) {
  SNP1 <- res.sub$SNP1[i]
  SNP2 <- res.sub$SNP2[i]
  if (i == 1) {
    res.sub.save <- res.sub[i,]
    SNP_list <- c(SNP_list,SNP1,SNP2)
  } else if (!(SNP1 %in% SNP_list) & !(SNP2 %in% SNP_list)) {
    res.sub.save <- rbind(res.sub.save,res.sub[i,])
    SNP_list <- c(SNP_list,SNP1,SNP2)
  }
}
s <- '20'
f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/results/full_data_gxe.',s,'.txt')
df <- fread(f,data.table = F,stringsAsFactors = F)

score.gxg <- 0
for (i in 1:nrow(res.sub.save)) {
  # for (i in 1:5) {
  BETA.GXG <- res.sub.save$BETA_INT[i]
  SNP1 <- str_replace_all(res.sub.save$SNP1[i],'-','.')
  SNP2 <- str_replace_all(res.sub.save$SNP2[i],'-','.')
  score.gxg <- score.gxg + BETA.GXG*df[,SNP1]*df[,SNP2]
}
cor(score.gxg,df[,phenoName2],use='p')

df$PGS.epi <- score.gxg

df.save <- merge(df,score.save[,c('IID','SCORE')],by='IID')

phenoName <- 'lymphocyte.count.rint.ALL'
x <- strsplit(phenoName,'\\.')[[1]]; phenoName2 <- paste0(c(x[-length(x)],'na'),collapse='.')

cor(df.save$PGS.epi,df.save[,phenoName2],use='p')
cor(df.save$SCORE,df.save[,phenoName2],use='p')

mod.formula <- paste0(phenoName2,'~','PGS.epi+SCORE')
mod.formula <- formula(mod.formula)
mod.epi <- lm(mod.formula,data=df.save)
# summary(mod.epi)
mod.formula <- paste0(phenoName2,'~','SCORE')
mod.formula <- formula(mod.formula)
mod <- lm(mod.formula,data=df.save)
# summary(mod)
summary(mod.epi)$adj.r.squared/summary(mod)$adj.r.squared


mod$adj.r.squared



mod.formula <- paste0(phenoName2,'~','PGS.epi')
mod.formula <- formula(mod.formula)
mod <- lm(mod.formula,data=df.save)
summary(mod)$coef


mod.formula <- paste0(phenoName2,'~','PGS.epi+SCORE+age+as.factor(menopause2) ')
mod.formula <- formula(mod.formula)
mod <- lm(mod.formula,data=df.save)
summary(mod)



cor(df.save[,c('SCORE','PGS.epi')],use='p')



