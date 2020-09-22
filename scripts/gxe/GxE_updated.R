library(data.table)
library(BEDMatrix)
source('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/bin/rntransform.R')
use_rint <- TRUE

for (s in c('80','20')) {
  
  environmental_factors <- c('SB',
                             'PA',
                             'Smoking.E',
                             'age',
                             'sex',
                             'Alcohol_intake_frequency'
                             )
  f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/full_data.',s,'.txt')
  fam <- fread(f,data.table = F,stringsAsFactors = F)
  
  pheno <- 'bmi'
  fam[,paste0(pheno,'.na')] <- fam[,pheno]
  fam[,paste0(pheno,'.na')][which(fam$In==0)] <- NA
  fam[,paste0(pheno,'.na')][which(fam$QC_In==0)] <- NA
  if (use_rint) {
    fam[,paste0(pheno,'.na.RINT')] <- rntransform(fam[,paste0(pheno,'.na')])
    phenoName <- paste0(pheno,'.na.RINT')
  } else {
    phenoName <- paste0(pheno,'.na')
  }
  
  # covariates data
  covariate_dataset <- fam[,c('IID','age','sex','age2','genotyping.array',
                              paste0('PC',1:20))]
  
  # environmental factors
  f <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/envir_data.txt'
  environmental_dataset <- fread(f,data.table = F,stringsAsFactors =F )
  colnames(environmental_dataset)[1] <- 'IID'
  environmental_factors_to_keep <- colnames(environmental_dataset)[-1]
  
  # environmental_factors_to_keep <- colnames(environmental_dataset)[colnames(environmental_dataset) %in% environmental_factors]
  # environmental_dataset <- environmental_dataset[,c("IID",environmental_factors_to_keep)]
  # environmental_factors_to_keep <- colnames(environmental_dataset)[colnames(environmental_dataset) %in% environmental_factors]
  # environmental_dataset <- environmental_dataset[,c('IID',environmental_factors_to_keep)]
  
  # rearrange and finalize
  colnames(environmental_dataset) <- gsub("-| |/",'_',c('IID',environmental_factors_to_keep))
  # environmental_factors <- gsub("-| |/",'_',c(environmental_factors))
  
  # merge covariate & environmental
  covariate_environmental_dataset <- merge(covariate_dataset,environmental_dataset,by='IID')
  
  # merge phenotypic data
  phenotype_dataset <- fam[,c('IID',phenoName)] # read in all phenotypes at once? 
  full_dataset <- merge(covariate_environmental_dataset,phenotype_dataset,by='IID')

  # genetic data
  f.geno <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.',pheno,'.merged_subset2')
  geno <- BEDMatrix(f.geno)

  geno_names <- unlist(lapply(strsplit(rownames(geno),'_'),function(x) {return(x[2])}))
  
  # merge stuff
  ind <- which(geno_names %in% full_dataset$IID)
  full_dataset <- full_dataset[match(geno_names[ind],full_dataset$IID),]
  fwrite(full_dataset,paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/full_data.',pheno,'.GxE.',s,'.txt'),quote = F,sep = '\t',na = 'NA',col.names = T,row.names = F)
  
  environmental_factors <- c('age','Alcohol_intake_frequency',
                             'PA','SB','sex','Smoking.E')

  library(parallel)
  GxE <- function(i) {
    print(i)
    full_dataset$SNP <- geno[ind,i]
    mod.formula <- formula(paste(phenoName,' ~ age+age2+genotyping.array+sex+age*sex+age2*sex+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
               PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+',
                                 envir_name,'*SNP'))
    mod <- lm(mod.formula,data=full_dataset)
    res <- summary(mod)$coef[paste0(envir_name,':SNP'),c("Estimate","Pr(>|t|)")]
    return(res)
  }
  
  Fit_Model <- function(start=1,p=5000) {
    df.save <- mclapply(start:p,GxE,mc.cores=8)
    df.save <- do.call(rbind,df.save)
    df.save <- as.data.frame(df.save)
    df.save$SNP <- colnames(geno)[start:p]
    colnames(df.save)[1:2] <- c('Estimate','Pr(>|t|)')
    df.save <- df.save[,c(3,1:2)]
    return(df.save)
  }
  
  for (k in 2:length(environmental_factors)) {
    envir_name <- environmental_factors[k]
    df.results <- Fit_Model(1,length(colnames(geno)))
    # df.results <- Fit_Model(60,60)
    df.results$E <- envir_name
    if (k==1) {
      df.results.save <- df.results
    } else {
      df.results.save <- rbind(df.results.save,df.results)
    }
  }
  if (use_rint) {
    fwrite(df.results.save,paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.',s,'.ext.more_snp.RINT.txt'),quote = F,sep = '\t',na = 'NA',row.names = F,col.names = T)
  } else {
    fwrite(df.results.save,paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.',s,'.ext.more_snp.txt'),quote = F,sep = '\t',na = 'NA',row.names = F,col.names = T)
  }
}




