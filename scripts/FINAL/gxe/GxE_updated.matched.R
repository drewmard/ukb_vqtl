library(data.table)
library(BEDMatrix)
source('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/bin/rntransform.R')

for (s in c('80','20')) {
  
  environmental_factors <- c('SB',
                             'PA',
                             'Smoking.E',
                             'age',
                             'sex',
                             paste0('DIET_PC',1:10),
                             c('Cooked vegetable intake',
                               'Salad / raw vegetable intake',
                               'Fresh fruit intake',
                               'Dried fruit intake',
                               'Bread intake',
                               'Cereal intake',
                               'Tea intake',
                               'Coffee intake',
                               'Water intake',
                               'Oily fish intake',
                               'Non-oily fish intake',
                               'Processed meat intake',
                               'Poultry intake',
                               'Beef intake',
                               'Lamb/mutton intake',
                               'Pork intake',
                               'Cheese intake',
                               'Salt added to food',
                               'Variation in diet',
                               'Alcohol intake frequency',
                               c('produce_intake','veggie_intake','fruit_intake','tea_coffee_intake','fish_intake','meat_intake','red_meat_intake','fish_meat_intake')
                             ))
  # s <- '20'
  f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/full_data.',s,'.txt')
  fam <- fread(f,data.table = F,stringsAsFactors = F)
  
  pheno <- 'bmi'
  # pheno <- 'lymphocyte.count'
  fam[,paste0(pheno,'.na')] <- fam[,pheno]
  fam[,paste0(pheno,'.na')][which(fam$In==0)] <- NA
  fam[,paste0(pheno,'.na')][which(fam$QC_In==0)] <- NA
  # i.outlier <- which(abs(scale(fam[,paste0(pheno,'.na')])) > 6)
  # fam[,paste0(pheno,'.na')][i.outlier] <- NA
  phenoName <- paste0(pheno,'.na')
  
  # phenoName <- 'bmi.log.ALL'
  # fam[,phenoName] <- log(fam[,paste0(pheno,'.na')])
  
  # covariates data
  covariate_dataset <- fam[,c('IID','age','sex','age2','genotyping.array',
                              paste0('PC',1:20))]
  
  # environmental factors
  # need to fix: remove individuals w/ missing data??????? yes.
  f <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/envir_data.txt'
  # f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/results/full_data_gxe.',s,'.txt') #### new file needed
  environmental_dataset <- fread(f,data.table = F,stringsAsFactors =F )
  colnames(environmental_dataset)[1] <- 'IID'
  environmental_factors_to_keep <- colnames(environmental_dataset)[colnames(environmental_dataset) %in% environmental_factors]
  environmental_dataset <- environmental_dataset[,c("IID",environmental_factors_to_keep)]
  
  # diet PCA data
  pca <- fread('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/diet_pca.txt',data.table = F,stringsAsFactors = F)
  colnames(pca)[1] <- 'IID'
  environmental_factors_to_keep <- colnames(pca)[colnames(pca) %in% environmental_factors]
  environmental_dataset <- merge(environmental_dataset,pca[,c("IID",environmental_factors_to_keep)],by='IID')
  
  # diet covariate data
  diet_covariate <- fread('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/diet_covariates.txt',data.table = F,stringsAsFactors = F)
  col_names <- c('Cooked vegetable intake','Salad / raw vegetable intake',
                 'Fresh fruit intake','Dried fruit intake'); diet_covariate[,'produce_intake'] <- apply(diet_covariate[,col_names],1,sum)
  col_names <-   c('Cooked vegetable intake',
                   'Salad / raw vegetable intake'); diet_covariate[,'veggie_intake'] <- apply(diet_covariate[,col_names],1,sum)
  col_names <- c('Fresh fruit intake',
    'Dried fruit intake'); diet_covariate[,'fruit_intake'] <- apply(diet_covariate[,col_names],1,sum)
  col_names <- c('Tea intake',
                 'Coffee intake'); diet_covariate[,'tea_coffee_intake'] <- apply(diet_covariate[,col_names],1,sum)
  col_names <- c('Oily fish intake',
                 'Non-oily fish intake'); diet_covariate[,'fish_intake'] <- apply(diet_covariate[,col_names],1,sum)
  col_names <- c('Processed meat intake',
                 'Poultry intake',
                 'Beef intake',
                 'Lamb/mutton intake',
                 'Pork intake'); diet_covariate[,'meat_intake'] <- apply(diet_covariate[,col_names],1,sum)
  col_names <- c('Processed meat intake',
                 'Beef intake',
                 'Lamb/mutton intake',
                 'Pork intake'); diet_covariate[,'red_meat_intake'] <- apply(diet_covariate[,col_names],1,sum)
  col_names <- c('Oily fish intake',
                 'Non-oily fish intake',
                 'Processed meat intake',
                 'Poultry intake',
                 'Beef intake',
                 'Lamb/mutton intake',
                 'Pork intake'); diet_covariate[,'fish_meat_intake'] <- apply(diet_covariate[,col_names],1,sum)
  
  colnames(diet_covariate)[1] <- 'IID'
  environmental_factors_to_keep <- colnames(diet_covariate)[colnames(diet_covariate) %in% environmental_factors]
  environmental_dataset <- merge(environmental_dataset,diet_covariate[,c("IID",environmental_factors_to_keep)],by='IID')
  
  # rearrange and finalize
  environmental_factors_to_keep <- colnames(environmental_dataset)[colnames(environmental_dataset) %in% environmental_factors]
  environmental_dataset <- environmental_dataset[,c('IID',environmental_factors_to_keep)]
  colnames(environmental_dataset) <- gsub("-| |/",'_',c('IID',environmental_factors_to_keep))
  environmental_factors <- gsub("-| |/",'_',c(environmental_factors))
  
  # lastly: want diet score:
  # ????
  
  # merge covariate & environmental
  covariate_environmental_dataset <- merge(covariate_dataset,environmental_dataset,by='IID')
  
  # merge phenotypic data
  phenotype_dataset <- fam[,c('IID',phenoName)] # read in all phenotypes at once? 
  full_dataset <- merge(covariate_environmental_dataset,phenotype_dataset,by='IID')

  # genetic data
  # f.geno <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.',pheno,'.merged_subset2')
  # geno <- BEDMatrix(f.geno)
  f.geno <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.',pheno,'.QTL_matched_snp.merged_subset')
  geno <- BEDMatrix(f.geno)
  
  geno_names <- unlist(lapply(strsplit(rownames(geno),'_'),function(x) {return(x[2])}))
  
  # merge stuff
  ind <- which(geno_names %in% full_dataset$IID)
  full_dataset <- full_dataset[match(geno_names[ind],full_dataset$IID),]
  fwrite(full_dataset,paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/full_data.',pheno,'.GxE.',s,'.txt'),quote = F,sep = '\t',na = 'NA',col.names = T,row.names = F)
  
  environmental_factors <- c(
                            # paste0('DIET_PC',1:10),
                             # 'DIET_SCORE',
                             'age','Alcohol_intake_frequency',
                             'PA','SB','sex','Smoking.E')
  
  # can loop through and call covariates if necessary
  covariate_adjusted_phenotype <- FALSE
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
    df.save <- mclapply(start:p,GxE,mc.cores=16)
    df.save <- do.call(rbind,df.save)
    df.save <- as.data.frame(df.save)
    df.save$SNP <- colnames(geno)[start:p]
    colnames(df.save)[1:2] <- c('Estimate','Pr(>|t|)')
    df.save <- df.save[,c(3,1:2)]
    return(df.save)
  }
  
  for (k in 1:length(environmental_factors)) {
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
  
  # fwrite(df.results.save,paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.',s,'.ext.more_snp.txt'),quote = F,sep = '\t',na = 'NA',row.names = F,col.names = T)
  fwrite(df.results.save,paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.',s,'.ext.more_snp.QTL_matched_snp.txt'),quote = F,sep = '\t',na = 'NA',row.names = F,col.names = T)
}





