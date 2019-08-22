library(data.table)
library(sandwich)
library(lmtest)

# initialize
user_direc <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl'
phenoName <- 'lymphocyte.count'
f.geno <- paste0(user_direc,'/output/GWAS/subset/',phenoName,'/ukbb.ALL_vQTL.raw')
df.geno <- fread(f.geno,data.table = F,stringsAsFactors = F)

# read in data
fam2.80 <- fread('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/full_data.80.txt',data.table = F,stringsAsFactors = F)
fam2.20 <- fread('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/full_data.20.txt',data.table = F,stringsAsFactors = F)

# var hits
# f <- paste0(user_direc,'/output/GWAS/results2/',phenoName,'/ukbb.',phenoName,'.results.var.clumped.cut.txt')
# index <- fread(f,data.table = F,stringsAsFactors = F)

# mean hits
# f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/UKB/results/ukbb.',phenoName,'.results.mean.clumped.cut.txt')
# index <- fread(f,data.table = F,stringsAsFactors = F)

# for (s in c('80','20')) {
for (s in c('20')) {
    
  print("reading in phenotypes & environmental data, then merging...")
  
  if (s=='80') {
    fam2 <- fam2.80 
  } else if (s=='20') {
    fam2 <- fam2.20
  }
  
  fam2$menopause2 <- as.factor(fam2$menopause2)
  
  source('/home/anm2868/scripts/Useful_scripts/rntransform.R')
  PHENOTYPE_NAMES <- c('lymphocyte.count','monocyte.count','neutrophil.count','neutrophil.percentage','wbc.leukocyte.count')
  for (k in 1:length(PHENOTYPE_NAMES)) {
    phenoName <- PHENOTYPE_NAMES[k]
    
    # remove outliers
    i.outlier <- which(abs(scale(fam2[,paste0(phenoName,'.na')])) > 5)
    fam2[,paste0(phenoName,'.na')][i.outlier] <- NA
    
    # apply transformations
    suffix<-'.rint'; fam2[,paste0(phenoName,'.na',suffix)] <- rntransform(fam2[,paste0(phenoName,'.na')])
    suffix<-'.log'; fam2[,paste0(phenoName,'.na',suffix)] <- log10(fam2[,paste0(phenoName,'.na')]+1)
  }
  #########################################
  
  df.pheno <- fam2
  
  # P data:
  # need to switch phenotype file directory:
  df <- merge(df.geno,df.pheno,by='IID')
  
  # E data:
  df.envir <- fread('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/envir_data.txt',data.table = F,stringsAsFactors = F)
  ENVIR_NAMES <- c(colnames(df.envir)[-1],'age','sex') #,'bmi2','menopause','time.since.period2')
  # df.envir <- df.envir[,-which(colnames(df.envir) %in% c('sex','age'))]
  df2 <- merge(df,df.envir,by.x='IID',by.y='eid')
  
  # PHENOTYPE_NAMES <- paste0(rep(PHENOTYPE_NAMES,each=3),'.na.',c('','log','rint'))
  # PHENOTYPE_NAMES <- paste0(rep(PHENOTYPE_NAMES,each=1),'.na.',c('rint'))
  # PHENOTYPE_NAMES <- paste0(PHENOTYPE_NAMES[1:2],'.na') # if only want to run specific phenotypes
  PHENOTYPE_NAMES <- 'lymphocyte.count.na.rint' # if only want to run one phenotype
  
  GxE_acrossPhenotypes <- function(i) {
    phenoName <- PHENOTYPE_NAMES[i]
    print(paste0('Phenotype ',i,'/',length(PHENOTYPE_NAMES),': ',phenoName )) # for debugging
    
    GxE <- function(j) {
      ENVIR_FACTOR <- ENVIR_NAMES[j]
      print(paste0('Environmental factor ',j,'/',length(ENVIR_NAMES),': ',ENVIR_FACTOR)) # for debugging
      
      # Create LM model w/ GxE interaction:
      # might need to change column names: for example, ENVIR_FACTOR will be smoking in GxE testing and need to impute medians to keep full data
      mod.formula <- paste(paste0(phenoName),' ~ age+age2+genotyping.array+sex+
                                   PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
                                   PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+
                                   time.since.period2+time.since.period2.dummy+
                                   menopause2+bmi2.dummy+bmi2+bmi2*age+',
                           vQTL,'*',ENVIR_FACTOR)
      if (ENVIR_FACTOR!='Smoking.E') {
        mod.formula <- paste0(mod.formula,'+Smoking+Smoking.dummy')
      } else if (ENVIR_FACTOR!='alcohol.freq.E') { 
        mod.formula <- paste0(mod.formula,'+alcohol.freq2+alcohol.freq2.dummy')
      }
      mod.formula <- formula(mod.formula)
      
      # need to switch glm
      # if (i != 1) {
      #   mod1 <- glm(mod.formula,
      #               data=df2,na.action=na.exclude,family=binomial(link="logit"))
      #   mod1.coef <- summary(mod1)$coef
      # } else {
      
      mod1 <- tryCatch({
        lm(mod.formula,
           data=df2,na.action=na.exclude)
      },error=function(e) {
        print('wtf????')
        mod.formula <- formula(paste(paste0(phenoName,'.na'),' ~ sex+age+age2+age*sex+age2*sex+
                                 PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
                                 PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+',
                                     vQTL,'*',ENVIR_FACTOR))
        mod1 <- lm(mod.formula,
                   data=df2,na.action=na.exclude)
        return(mod1)
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
  for (k in 1:nrow(index)) {
    vQTL=index[k,2]
    print(paste0('SNP ',k,'/',nrow(index),': ',vQTL)) # for debugging
    y2 <- (lapply(1:length(PHENOTYPE_NAMES),GxE_acrossPhenotypes))
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
}
results[order(results$P_GxE),][1:10,]

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