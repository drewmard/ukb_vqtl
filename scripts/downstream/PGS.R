library(data.table)

SNP <- 'rs887468'
ENV <- 'age'

s <- '80'
# s <- '20'
phenoName <- 'lymphocyte.count.rint.ALL'
f.res <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/results/ukbb.gxe.',phenoName,'.',s,'.txt')
# f.res <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/results/ukbb.gxe.',s,'.txt')
results <- fread(f.res,data.table = F,stringsAsFactors = F)
results[order(results$P_GxE),][1,]

s <- '20'
f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/results/full_data_gxe.',s,'.txt')
df <- fread(f,data.table = F,stringsAsFactors = F)

mod.formula.1 <- paste0('lymphocyte.count.rint.na',' ~ genotyping.array+
       PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
                        PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+
                        menopause2+
                        bmi2.dummy+bmi2+
                        Smoking+Smoking.dummy+
                        alcohol.freq2+alcohol.freq2.dummy')
mod.formula.1.ext <- paste0(mod.formula.1,'+age+',SNP)
mod.formula.1.ext <- formula(mod.formula.1.ext)
mod <- lm(mod.formula.1.ext,data=fam3)
mod.sum <- summary(mod)$coef

BETA.ENV <- mod.sum[ENV,1]
BETA.SNP <- mod.sum[SNP,1]
score <- BETA.ENV*df[,ENV] + 
  BETA.SNP*df[,SNP]
cor(score,df$lymphocyte.count.rint.na,use='p')

# BETA.SNP.gxe <- -0.144
# BETA.ENV.gxe <- 0.00673
# BETA.GxE.gxe <- 0.00180
x <- subset(results,vQTL==SNP & E==ENV)

BETA.SNP.gxe <- x[1,'BETA_vQTL']
BETA.ENV.gxe <- x[1,'BETA_E']
BETA.GxE.gxe <- x[1,'BETA_GxE']

score.gxe <- BETA.GxE.gxe*df[,SNP]*df[,ENV] + 
  BETA.ENV.gxe*df[,ENV] + 
  BETA.SNP.gxe*df[,SNP]
cor(score.gxe,df$lymphocyte.count.rint.na,use='p')


