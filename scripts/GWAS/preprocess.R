library(data.table)

# read in data
fam2.80 <- fread('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/full_data.80.txt',data.table = F,stringsAsFactors = F)
fam2.20 <- fread('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/full_data.20.txt',data.table = F,stringsAsFactors = F)

# Create phenotype file
phenotypeDataFile.80 <- fam2.80[,c('FID','IID')]
phenotypeDataFile.20 <- fam2.20[,c('FID','IID')]

#initialize
# PHENOTYPE_NAMES <- c('lymphocyte.count','monocyte.count','neutrophil.count','neutrophil.percentage','wbc.leukocyte.count')
PHENOTYPE_NAMES <- c('lymphocyte.count','monocyte.count','neutrophil.count','neutrophil.percentage','wbc.leukocyte.count',
                     'rbc.erythrocyte.count','platelet.count','eosinophil.count','basophil.count')


############

for (s in c('80','20')) {
  
  #########################################
  # identify data set #####################
  if (s=='80') {
    fam2 <- fam2.80 
    phenotypeDataFile <- phenotypeDataFile.80
  } else if (s=='20') {
    fam2 <- fam2.20
    phenotypeDataFile <- phenotypeDataFile.20
  }
  #########################################
  fam2$menopause2[which(fam2$sex==1)] <- 4
  fam2$menopause2 <- as.factor(fam2$menopause2)
  # table(fam2[,c('sex','menopause2')])
  #########################################
  
  
  # phenotype preprocessing ###############
  source('/home/anm2868/scripts/Useful_scripts/rntransform.R')
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
  
  #########################################
  # fit model to no covariate outliers but measure residuals in all individuals
  fam3 <- fam2
  fam3$bmi2 <- fam3$bmi2.with_outliers
  # fam3$time.since.period2 <- fam3$time.since.period2.with_outliers
  #########################################  

  # run modeling
  suffix='';k=1
  for (suffix in c('','.log','.rint')) {
    
    for (k in 1:length(PHENOTYPE_NAMES)) {
      
      phenoName <- PHENOTYPE_NAMES[k]
      
      # menopause sex related variable
      i <- which(fam2$sex==1); j <- which(fam2$sex==0); 
      mod.formula.1 <- formula(paste(paste0(phenoName,'.na',suffix),' ~ age+age2+genotyping.array+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
               PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+
               Smoking+Smoking.dummy+alcohol.freq2+alcohol.freq2.dummy+bmi2.dummy+bmi2'))
      mod.formula.2 <- formula(paste(paste0(phenoName,'.na',suffix),' ~ age+age2+genotyping.array+
                PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
                PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+
                 menopause2+
                 Smoking+Smoking.dummy+alcohol.freq2+alcohol.freq2.dummy+bmi2.dummy+bmi2'))

      # fitting models, not including individuals that are covariate outliers
      mod1 <- lm(mod.formula.1,
                 data=fam2[i,],na.action=na.exclude)
      mod2 <- lm(mod.formula.2,
                 data=fam2[j,],na.action=na.exclude)
      mod3 <- lm(mod.formula.2,
                 data=fam2,na.action=na.exclude)
      
      # but calculate residuals in all individuals
      mod1.pred <- predict.lm(mod1,newdata = fam3[i,])
      mod2.pred <- predict.lm(mod2,newdata = fam3[j,])
      mod3.pred <- predict.lm(mod3,newdata = fam3)
      resid1 <- fam3[i,paste0(phenoName,'.na',suffix)] - mod1.pred
      resid2 <- fam3[j,paste0(phenoName,'.na',suffix)] - mod2.pred
      resid3 <- fam3[,paste0(phenoName,'.na',suffix)] - mod3.pred
      resid1 <- scale(resid1)
      resid2 <- scale(resid2)
      resid3 <- scale(resid3)
      
      
      # Supplement to phenotype file, and switch NA to -9 for plink:
      x <- paste0(phenoName,suffix,'.M')
      phenotypeDataFile[,x] <- NA
      phenotypeDataFile[,x][i] <- resid1
      phenotypeDataFile[,x][which(is.na(phenotypeDataFile[,x]))] <- -9
      x <- paste0(phenoName,suffix,'.F')
      phenotypeDataFile[,x] <- NA
      phenotypeDataFile[,x][j] <- resid2
      phenotypeDataFile[,x][which(is.na(phenotypeDataFile[,x]))] <- -9
      x <- paste0(phenoName,suffix,'.ALL')
      phenotypeDataFile[,x] <- resid3
      phenotypeDataFile[,x][which(is.na(phenotypeDataFile[,x]))] <- -9
      
    }
  }
  fwrite(phenotypeDataFile,paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/phenotypes_processed.',s,'.txt'),sep='\t',quote = F,col.names = T,row.names = F)
}

# library(data.table)
# s <- '80'
# phenotypeDataFile <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/phenotypes_processed.',s,'.txt'),data.table = F,stringsAsFactors = F)
fwrite(as.data.frame(colnames(phenotypeDataFile)[-c(1:2)]),'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/phenotype_names.txt',col.names = F,row.names = F,sep = '\t',quote = F)



