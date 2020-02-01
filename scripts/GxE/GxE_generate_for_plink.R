library(data.table)
source('/home/anm2868/scripts/Useful_scripts/rntransform.R')

library(data.table)

# read in data
fam2.80 <- fread('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/full_data.80.txt',data.table = F,stringsAsFactors = F)
fam2.20 <- fread('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/full_data.20.txt',data.table = F,stringsAsFactors = F)

# Create phenotype file
phenotypeDataFile.80 <- fam2.80[,c('FID','IID')]
phenotypeDataFile.20 <- fam2.20[,c('FID','IID')]

#initialize
# PHENOTYPE_NAMES <- c('lymphocyte.count','monocyte.count','neutrophil.count','neutrophil.percentage','wbc.leukocyte.count',
#                      'rbc.erythrocyte.count','platelet.count','eosinophil.count','basophil.count')
PHENOTYPE_NAMES <- c('lymphocyte.count','monocyte.count','neutrophil.count','eosinophil.count','basophil.count')


############

for (s in c('80','20')) {
  
  #########################################
  # identify data set #####################
  if (s=='80') {
    fam2 <- fam2.80 
    fam2 <- fam2[-which(duplicated(fam2$IID)),]
    phenotypeDataFile <- phenotypeDataFile.80
    pheno_df <- fam2[,c('FID','IID')]
    
  } else if (s=='20') {
    fam2 <- fam2.20
    fam2 <- fam2[-which(duplicated(fam2$IID)),]
    phenotypeDataFile <- phenotypeDataFile.20
    pheno_df <- fam2[,c('FID','IID')]
    
  }
  #########################################
  fam2$menopause2[which(fam2$sex==1)] <- 4
  fam2$menopause2 <- as.factor(fam2$menopause2)

  # E data:
  df.envir <- fread('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/envir_data.txt',data.table = F,stringsAsFactors = F)
  ENVIR_NAMES <- c(colnames(df.envir)[-1],'age','sex') #,'bmi2','menopause','time.since.period2')
  # df.envir <- df.envir[,-which(colnames(df.envir) %in% c('sex','age'))]
  fam2 <- merge(fam2,df.envir,by.x='IID',by.y='eid')
  
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
    # suffix<-'.log'; fam2[,paste0(phenoName,'.na',suffix)] <- log10(fam2[,paste0(phenoName,'.na')]+1)
    
  }
  #########################################
  
  #########################################
  # fit model to no covariate outliers but measure residuals in all individuals
  fam3 <- fam2
  fam3$bmi2 <- fam3$bmi2.with_outliers
  # fam3$time.since.period2 <- fam3$time.since.period2.with_outliers
  #########################################  
  
  for (phenoName in PHENOTYPE_NAMES) {
    for (iter in 1:4) {
      mod.formula <- paste0(paste0(phenoName,'.na','.rint'),' ~ genotyping.array+
                              PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
                              PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+
                              menopause2+
                              bmi2.dummy+bmi2')
      if (iter!=2) { #(ENVIR_FACTOR!='Smoking.E') {
        mod.formula <- paste0(mod.formula,'+Smoking+Smoking.dummy')
      } 
      
      if (iter!=3) { #if (ENVIR_FACTOR!='alcohol.freq.E') { 
        mod.formula <- paste0(mod.formula,'+alcohol.freq2+alcohol.freq2.dummy')
      }
      if (iter!=4) {
        mod.formula <- paste0(mod.formula,'+age+age2')
      }
      
      mod.formula <- formula(mod.formula)
      mod3 <- lm(mod.formula,
                 data=fam2,na.action=na.exclude)
      mod3.pred <- predict.lm(mod3,newdata = fam3)
      resid3 <- fam3[,paste0(phenoName,'.na',suffix)] - mod3.pred
      resid3 <- scale(resid3)
      
      pheno_df[,paste0(phenoName,'.rint','.ALL','.resid.',iter)] <- resid3
      pheno_df[,paste0(phenoName,'.rint','.ALL','.resid.',iter)][which(is.na(pheno_df[,paste0(phenoName,'.rint','.ALL','.resid.',iter)]))] <- -9
      
      # if (iter==2) {
      #   fam2
      # }
      # 
    }
  }
  f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/EnvData/pheno.',s,'.txt')
  fwrite(pheno_df,file=f,sep='\t',quote = F,col.names = T,row.names = F)
  
  dir.create('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/EnvData',showWarnings = F)
  for (i in 1:length(ENVIR_NAMES)) {
    fwrite(
      fam2[,c('FID','IID',ENVIR_NAMES[i])],
      paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/EnvData/',ENVIR_NAMES[i],'.',s,'.txt'),
      sep='\t',quote = F,col.names = T,row.names = F,na='NA'
    )
  }
  
}




  
  