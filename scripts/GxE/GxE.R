library(data.table)
library(sandwich)
library(lmtest)

source('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/scripts/GxE/gen_datafiles.R')

# initialize
# phenoName <- 'lymphocyte.count.rint'
phenoName <- 'lymphocyte.count.rint.ALL'
# phenoName <- 'lymphocyte.count.ALL'
phenoName2 <- 'lymphocyte.count'

PHENOTYPE_NAMES <- phenoName

# PHENOTYPE_NAMES <- c('lymphocyte.count','monocyte.count','neutrophil.count','neutrophil.percentage','wbc.leukocyte.count')
# ENVIR_NAMES <- c("PA","SB","Smoking.E","sleep.duration","getting.up.morning","nap.during.day",
#                  "time.spent.outdoors.summer","time.spent.outdoors.winter","time.spent.outdoors",
#                  "tobacco.smoke.exposure","alcohol.freq.E","childhood.sunburn.occasions",
#                  "age.started.smoking","hormone.replacement.therapy",
#                  'age')#,
                 # 'sex')
ENVIR_NAMES <- c("PA","SB","Smoking.E","sleep.duration","alcohol.freq.E",
                 'age')


#,'bmi2','menopause','time.since.period2')
# ENVIR_NAMES <- c('PA','SB','Smoking.E','alcohol.freq.E','childhood.sunburn.occasions','age','sex') #,'bmi2','menopause','time.since.period2')

# loading:
source('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/scripts/GxE/gen_datafiles.R');
gen_datafiles(phenoName,phenoName2)
# PHENOTYPE_NAMES <- paste0(rep(PHENOTYPE_NAMES,each=3),c('','.log','.rint'),'.na')
# phenoName2 <- paste0(phenoName,'.na')

for (s in c('80','20')) {
# for (s in c('80')) {
  
  s <- '80'
  
  f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/results/full_data_gxe.',s,'.txt')
  df2 <- fread(f,data.table = F,stringsAsFactors = F)
  df2$menopause2 <- as.factor(df2$menopause2)
  
  index <- colnames(df2)[2:(which(colnames(df2)=='FID')-1)]

  fam3 <- df2
  fam3$bmi2 <- df2$bmi2.with_outliers
  fam3$time.since.period2 <- df2$time.since.period2.with_outliers
  
  source('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/scripts/GxE/GxE_acrossPhenotypes_2.R')
  source('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/scripts/GxE/GxE_acrossEnvirFactor.R')
  
  print('GxE...')
  START=TRUE
  lapply(1:2,GxE_acrossSNPs)
  
  
  # for (k in 112:135) {
  for (k in 1:length(index)) {
    vQTL <- index[k]
    print(paste0('SNP ',k,'/',length(index),': ',vQTL)) # for debugging
    y2 <- (lapply(which(PHENOTYPE_NAMES%in%
                          # c('monocyte.count.na','monocyte.count.rint.na')
                        c(paste0(phenoName))
                        ),GxE_acrossPhenotypes))
    y2.df <- as.data.frame(do.call(rbind,y2))
    for (l in 4:length(y2.df)) {
      y2.df[,l] <- as.numeric(levels(y2.df[,l]))[y2.df[,l]]
    }
    if (START) {
      results <- y2.df
      START<-FALSE
    } else {
      results <- rbind(results,y2.df)
    }
  }
  
  print(results[order(results$P_GxE),][1:10,])
  # subset(results[order(results$P_GxE),],E!='age')[1:3,]
  
  # SAVE # # # #
  f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/results/ukbb.gxe.',phenoName,'.',s,'.txt')
  fwrite(results,f.out,sep='\t',na='NA',row.names = F,col.names = T,quote = F)
  
  
}

########################################################################
########################################################################
########################################################################
########################################################################


f.in <- paste0('/athena/elementolab/scratch/anm2868/vQTL/UKB/results/ukbb.bmi.vqtl_gxe.txt')
results <- fread(f.in,data.table = F,stringsAsFactors = F)
vQTL <- 'rs1421085'
phenotype <- 'bmi'
envir.factor <- 'PA'
x <- subset(results,Phenotype=='bmi')
x.sub <- subset(x,P_GxE < 0.05/nrow(x)); x.sub

# library(plyr); df.aggre <- aggregate(df2[,paste0(phenotype,'.na')],by=list(df2[,vQTL],round_any(df2[,envir.factor],10)),mean,na.rm=T) # for age
df.aggre <- aggregate(df2[,paste0(phenotype,'.na')],by=list(df2[,vQTL],df2[,envir.factor]),mean,na.rm=T)
# aggregate(df2[,paste0(phenotype,'.na')],by=list(df2[,vQTL],df2[,envir.factor]),function(x) {sum(!is.na(x))})
colnames(df.aggre) <- c(vQTL,envir.factor,phenotype)
f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/UKB/subset/gxe/ukbb.',vQTL,'.',phenotype,'.aggre.txt')
fwrite(df.aggre,f.out,sep='\t',na='NA',row.names = F,col.names = T,quote = F)
df.save <- df2[,c(vQTL,envir.factor,
  paste0(phenotype,'.na')
)]
f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/UKB/subset/gxe/ukbb.',vQTL,'.',phenotype,'.txt')
fwrite(df.save,f.out,sep='\t',na='NA',row.names = F,col.names = T,quote = F)

# replicate in hypertension
sig <- paste(x.sub$vQTL,x.sub$E,sep=' - ')
y <- subset(results,Phenotype=='hypertension')
y.sub <- y[which(paste(y$vQTL,y$E,sep=' - ') %in% sig),]
hypertension.sig <- subset(y.sub,P_GxE < 0.05/nrow(y.sub)); hypertension.sig
# replicate in other phenotypes
subset(results,paste(vQTL,E,sep=' - ') %in% paste(hypertension.sig$vQTL,hypertension.sig$E,sep=' - ')) 

# extract beta
mod.formula <- formula(paste(paste0(phenotype,'.na'),' ~ sex+age+age2+age*sex+age2*sex+genotyping.array+
                             PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
                             PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+',
                             vQTL))
mod1 <- lm(mod.formula,
           data=df2,na.action=na.exclude)


mod.formula <- formula(paste(paste0(phenotype,'.na'),' ~ sex+age+age2+age*sex+age2*sex+genotyping.array+
                             PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
                             PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+',
                             envir.factor))
mod2 <- lm(mod.formula,
           data=df2,na.action=na.exclude)

mod.formula <- formula(paste(paste0(phenotype,'.na'),' ~ sex+age+age2+age*sex+age2*sex+genotyping.array+
                             PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
                             PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+',
                             vQTL,'+',envir.factor))
mod3 <- lm(mod.formula,
           data=df2,na.action=na.exclude)

f.in <- paste0('/athena/elementolab/scratch/anm2868/vQTL/UKB/results/ukbb.bmi.vqtl_gxe.txt')
results <- fread(f.in,data.table = F,stringsAsFactors = F); x1 <- vQTL; x2 <- envir.factor; mod4 <- subset(results,vQTL==x1 & E==x2 & Phenotype=='bmi')
mod.df <- data.frame(matrix(NA,nrow=3,ncol=4)); colnames(mod.df) <- c('model','vqtl','e','vqtl*e')
mod.df[1,] <- c('alone',as.numeric(mod1$coef[vQTL]),as.numeric(mod2$coef[envir.factor]),NA)
mod.df[2,] <- c('together',as.numeric(mod3$coef[vQTL]),as.numeric(mod3$coef[envir.factor]),NA)
mod.df[3,] <- c('interaction',as.numeric(mod4$BETA_vQTL),as.numeric(mod4$BETA_E),as.numeric(mod4$BETA_GxE))
mod.df






x <- subset(results,Phenotype=='bmi')
x <- x[,c('vQTL','E','P_GxE')]
x$P_GxE <- -log10(x$P_GxE)
fwrite(x,'/athena/elementolab/scratch/anm2868/vQTL/UKB/GxE.bmi.txt',sep='\t',na='NA',row.names = F,col.names = T,quote = F)

x <- subset(results,Phenotype!='bmi')
x.sub <- subset(x,P_GxE < 0.05/nrow(x))
x.sub