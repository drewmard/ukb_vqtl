library(data.table)

gen_datafiles <- function(phenoName) {
  
  # initialize
  # phenoName <- 'monocyte.count'
  phenoName2 <- paste0(phenoName,'.na')
  user_direc <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl'
  f.geno <- paste0(user_direc,'/output/GWAS/subset/',phenoName,'/ukbb.ALL_vQTL.raw')
  df.geno <- fread(f.geno,data.table = F,stringsAsFactors = F,check.names = T)
  index <- colnames(df.geno)[-1]
  
  # read in data
  fam2.80 <- fread('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/full_data.80.txt',data.table = F,stringsAsFactors = F)
  fam2.20 <- fread('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/full_data.20.txt',data.table = F,stringsAsFactors = F)
  
  for (s in c('80','20')) {
    
    print("reading in phenotypes & environmental data, then merging...")
    
    if (s=='80') {
      fam2 <- fam2.80 
    } else if (s=='20') {
      fam2 <- fam2.20
    }
    
    fam2$menopause2[which(fam2$sex==1)] <- 4
    fam2$menopause2 <- as.factor(fam2$menopause2)

    source('/home/anm2868/scripts/Useful_scripts/rntransform.R')
    PHENOTYPE_NAMES <- c('lymphocyte.count','monocyte.count','neutrophil.count','neutrophil.percentage','wbc.leukocyte.count')
    for (k in 1:length(PHENOTYPE_NAMES)) {
      phenoName <- PHENOTYPE_NAMES[k]
      
      # remove outliers
      i.outlier <- which(abs(scale(fam2[,paste0(phenoName,'.na')])) > 5)
      fam2[,paste0(phenoName,'.na')][i.outlier] <- NA
      
      # apply transformations
      suffix<-'.rint'; fam2[,paste0(phenoName,suffix,'.na')] <- rntransform(fam2[,paste0(phenoName,'.na')])
      suffix<-'.log'; fam2[,paste0(phenoName,suffix,'.na')] <- log10(fam2[,paste0(phenoName,'.na')]+1)
    }
    #########################################
    
    df.pheno <- fam2
    
    # P data:
    df <- merge(df.geno,df.pheno,by='IID')
    
    # E data:
    df.envir <- fread('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/envir_data.txt',data.table = F,stringsAsFactors = F)
    ENVIR_NAMES <- c(colnames(df.envir)[-1],'age','sex') #,'bmi2','menopause','time.since.period2')
    # df.envir <- df.envir[,-which(colnames(df.envir) %in% c('sex','age'))]
    df2 <- merge(df,df.envir,by.x='IID',by.y='eid')
    
    ########################################
    # build residual files
    fam3 <- df2
    fam3$bmi2 <- df2$bmi2.with_outliers
    fam3$time.since.period2 <- df2$time.since.period2.with_outliers
    
    for (iter in 1:4) {
      i <- which(df2$sex==1); j <- which(df2$sex==0); 
      mod.formula.1 <- (paste(paste0(phenoName),' ~ genotyping.array+
                              PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
                              PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+
                              menopause2+
                              bmi2.dummy+bmi2'))

      if (iter!=3) { #(ENVIR_FACTOR!='Smoking.E') {
        mod.formula.1 <- paste0(mod.formula.1,'+Smoking+Smoking.dummy')
      } 
      if (iter!=2) { #if (ENVIR_FACTOR!='alcohol.freq.E') { 
        mod.formula.1 <- paste0(mod.formula.1,'+alcohol.freq2+alcohol.freq2.dummy')
      }
      if (iter!=4) {
        mod.formula.1 <- paste0(mod.formula.1,'+age+age2')
      }
      mod.formula.1 <- formula(mod.formula.1)

      # fitting models, not including individuals that are outliers
      mod1 <- lm(mod.formula.1,
                 data=df2,na.action=na.exclude)

      # but calculate residuals in all individuals
      mod1.pred <- predict.lm(mod1,newdata = fam3)
      resid1 <- fam3[,paste0(phenoName)] - mod1.pred
      resid1 <- scale(resid1)

      # Supplement to phenotype file
      df2[,paste0('resid',iter)] <- resid1

    }

    f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/results/full_data_gxe.',s,'.txt')
    fwrite(df2,f.out,sep='\t',na='NA',row.names = F,col.names = T,quote = F)
    
  }
}
