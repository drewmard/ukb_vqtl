library(data.table)
library(sandwich)
library(lmtest)

source('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/scripts/GxE/gen_datafiles.R')

# initialize
phenoName <- 'monocyte.count.rint'
gen_datafiles(phenoName)
# phenoName2 <- 'monocyte.count.na.rint'

for (s in c('80','20')) {
  
  f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/results/full_data_gxe.',s,'.txt')
  df2 <- fread(f,data.table = F,stringsAsFactors = F)

  PHENOTYPE_NAMES <- paste0(rep(PHENOTYPE_NAMES,each=3),c('','.log','.rint'),'.na')

  GxE_acrossPhenotypes <- function(i) {
    phenoName <- PHENOTYPE_NAMES[i]
    print(paste0('Phenotype ',i,'/',length(PHENOTYPE_NAMES),': ',phenoName )) # for debugging
    
    GxE <- function(j) {
      ENVIR_FACTOR <- ENVIR_NAMES[j]
      print(paste0('Environmental factor ',j,'/',length(ENVIR_NAMES),': ',ENVIR_FACTOR)) # for debugging
      
      # Create LM model w/ GxE interaction:
      mod1 <- tryCatch({
        mod.formula <- paste(paste0(phenoName),' ~ age+age2+genotyping.array+sex+
                                   PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
                             PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+',
                             vQTL,'*',ENVIR_FACTOR)
        
        # time.since.period2+time.since.period2.dummy+menopause2+ # women only
        # bmi2.dummy+bmi2+bmi2*age+',
        if (ENVIR_FACTOR!='Smoking.E') {
          mod.formula <- paste0(mod.formula,'+Smoking+Smoking.dummy')
        } 
        if (ENVIR_FACTOR!='alcohol.freq.E') { 
          mod.formula <- paste0(mod.formula,'+alcohol.freq2+alcohol.freq2.dummy')
        }
        mod.formula <- formula(mod.formula)
        
        mod1 <- lm(mod.formula,
                   data=df2,na.action=na.exclude)
        mod1
      },error=function(e) {
        mod.formula <- paste(paste0(phenoName),' ~ age+age2+sex+
                                   PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
                             PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+',
                             vQTL,'*',ENVIR_FACTOR)
        
        # time.since.period2+time.since.period2.dummy+menopause2+ # women only
        # bmi2.dummy+bmi2+bmi2*age+',
        if (ENVIR_FACTOR!='Smoking.E') {
          mod.formula <- paste0(mod.formula,'+Smoking+Smoking.dummy')
        } 
        if (ENVIR_FACTOR!='alcohol.freq.E') { 
          mod.formula <- paste0(mod.formula,'+alcohol.freq2+alcohol.freq2.dummy')
        }
        mod.formula <- formula(mod.formula)
        
        mod1 <- lm(mod.formula,
                   data=df2,na.action=na.exclude)
        mod1
      })

      # mod1.coef <- summary(mod1)$coef
      mod1.coef <- coeftest(mod1, vcov = vcovHC(mod1))
      
      # }
      res1 <- mod1.coef[paste0(vQTL),]
      res2 <- mod1.coef[paste0(ENVIR_FACTOR),]
      res3 <- tryCatch({mod1.coef[paste0(vQTL,':',ENVIR_FACTOR),]},
                       error=function(e) { mod1.coef[paste0(ENVIR_FACTOR,':',vQTL),]})
      
      return(c(phenoName,vQTL,ENVIR_FACTOR,res1,res2,res3))
    }
    
    y <- lapply(1:length(ENVIR_NAMES),GxE)
    y.df <- as.data.frame(do.call(rbind,y))
    colnames(y.df) <- c('Phenotype','vQTL','E','BETA_vQTL','SE_vQTL','T_vQTL','P_vQTL',
                        'BETA_E','SE_E','T_E','P_E',
                        'BETA_GxE','SE_GxE','T_GxE','P_GxE')
    return(y.df)
  }
  
  print('GxE...')
  # df2 <- df.80
  # for (k in 1:nrow(index)) {
    # vQTL=index[k,2]
  for (k in 1:length(index)) {
    vQTL <- index[k]
    print(paste0('SNP ',k,'/',length(index),': ',vQTL)) # for debugging
    # y2 <- (lapply(1:length(PHENOTYPE_NAMES),GxE_acrossPhenotypes))
    # y2 <- (lapply(which(PHENOTYPE_NAMES==phenoName2),GxE_acrossPhenotypes))
    y2 <- (lapply(which(PHENOTYPE_NAMES=='monocyte.count.na'),GxE_acrossPhenotypes))
    y2.df <- as.data.frame(do.call(rbind,y2))
    for (l in 4:length(y2.df)) {
      y2.df[,l] <- as.numeric(levels(y2.df[,l]))[y2.df[,l]]
    }
    if (k==1) {
      results <- y2.df
    } else {
      results <- rbind(results,y2.df)
    }
  }
  
  # SAVE # # # #
  f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/results/ukbb.gxe.',s,'.txt')
  fwrite(results,f.out,sep='\t',na='NA',row.names = F,col.names = T,quote = F)
  
  print(results[order(results$P_GxE),][1:10,])
  
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