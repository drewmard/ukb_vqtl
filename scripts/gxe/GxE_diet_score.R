# diet score
# Sys.sleep(600)
use_rint <- FALSE
use_log <- TRUE

# recommend pre-processing using GxE_updated.R to get bmi raw or rint
library(data.table)
s='80'
full_dataset <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/full_data.bmi.GxE.',s,'.txt'),data.table = F,stringsAsFactors = F)
colnames(full_dataset) <- gsub("-| |/",'_',colnames(full_dataset))

pheno <- 'bmi'
if (use_rint) {
  # full_dataset[,paste0(pheno,'.na.RINT')] <- rntransform(full_dataset[,paste0(pheno,'.na')])
  phenoName <- paste0(pheno,'.na.RINT')
} else if (use_log) {
  phenoName <- paste0(pheno,'.na.log')
} else {
  phenoName <- paste0(pheno,'.na')
}

set.seed(031995)
i.20 <- sample(1:nrow(full_dataset),floor(0.25*nrow(full_dataset)),replace = F)
dataf.20 <- full_dataset[i.20,]
dataf.60 <- full_dataset[-i.20,]

# food_covariates <- c('Cooked vegetable intake',
#                      'Salad / raw vegetable intake',
#                      'Fresh fruit intake',
#                      'Dried fruit intake',
#                      'Bread intake',
#                      'Cereal intake',
#                      'Tea intake',
#                      'Coffee intake',
#                      'Water intake',
#                      'Oily fish intake',
#                      'Non-oily fish intake',
#                      'Processed meat intake',
#                      'Poultry intake',
#                      'Beef intake',
#                      'Lamb/mutton intake',
#                      'Pork intake',
#                      'Cheese intake',
#                      'Salt added to food')
# food_covariates <- gsub("-| |/",'_',c(food_covariates))

food_covariates <- c('cooked.veggie','raw.veggie','fresh.fruit','dried.fruit','oily.fish','non.oily.fish',
                      'processed.meat','poultry','beef','lamb','pork',
                      'cheese','milk','bread.intake','cereal.intake',
                      'salt.added.to.food','tea','coffee')

mod.formula.2 <- formula(paste(
  phenoName,
  ' ~ age+age2+genotyping.array+sex+age*sex+age2*sex+',
  paste(paste0('PC',1:20),collapse = '+'),'+',
  paste(food_covariates,collapse = '+')
))
mod.fit <- lm(mod.formula.2,data=dataf.20)
summary(mod.fit)

# hold out set
BETA <- as.matrix(coef(mod.fit)[paste0(food_covariates)])
DIET_SCORE <- scale(as.matrix(dataf.60[,paste0(food_covariates)]) %*% BETA)
dataf.60$DIET_SCORE <- DIET_SCORE

# genetic data
library(BEDMatrix)
f.geno <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.',pheno,'.merged_subset2')
geno <- BEDMatrix(f.geno)
geno_names <- unlist(lapply(strsplit(rownames(geno),'_'),function(x) {return(x[2])}))

# merge stuff
ind <- which(geno_names %in% dataf.60$IID)
dataf.60 <- dataf.60[match(geno_names[ind],dataf.60$IID),]
s='80'
if (use_rint) {
  f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/diet_data.',pheno,'.GxE.',s,'.RINT.txt')
} else if (use_log) {
  f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/diet_data.',pheno,'.GxE.',s,'.log.txt')
} else {
  f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/diet_data.',pheno,'.GxE.',s,'.txt')
}
fwrite(dataf.60[,c('IID','DIET_SCORE')],f,quote = F,na='NA',row.names = F,col.names = T,sep = '\t')

environmental_factors <- c('DIET_SCORE')
# can loop through and call covariates if necessary
covariate_adjusted_phenotype <- FALSE
library(parallel)
GxE <- function(i) {
  print(i)
  dataf.60$SNP <- geno[ind,i]
  mod.formula <- formula(paste(phenoName,' ~ age+age2+genotyping.array+sex+age*sex+age2*sex+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
               PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+',
                               envir_name,'*SNP'))
  mod <- lm(mod.formula,data=dataf.60)
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

if (use_rint) {
  f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.',s,'.diet_score.more_snp.RINT.txt')
} else if (use_log) {
  f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.',s,'.diet_score.more_snp.log.txt')
} else {
  f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.',s,'.diet_score.more_snp.txt')
}
fwrite(df.results.save,f,quote = F,sep = '\t',na = 'NA',row.names = F,col.names = T)

###### replication testing set

library(data.table)
s='20'
full_dataset <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/full_data.bmi.GxE.',s,'.txt'),data.table = F,stringsAsFactors = F)
colnames(full_dataset) <- gsub("-| |/",'_',colnames(full_dataset))

dataf.20_test <- full_dataset

# hold out set
BETA <- as.matrix(coef(mod.fit)[paste0(food_covariates)])
DIET_SCORE <- scale(as.matrix(dataf.20_test[,paste0(food_covariates)]) %*% BETA)
dataf.20_test$DIET_SCORE <- DIET_SCORE

# merge stuff
ind <- which(geno_names %in% dataf.20_test$IID)
dataf.20_test <- dataf.20_test[match(geno_names[ind],dataf.20_test$IID),]
s='20'; f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/diet_data.',pheno,'.GxE.',s,'.txt')
fwrite(dataf.20_test[,c('IID','DIET_SCORE')],f,quote = F,na='NA',row.names = F,col.names = T,sep = '\t')

library(parallel)
GxE <- function(i) {
  print(i)
  dataf.20_test$SNP <- geno[ind,i]
  mod.formula <- formula(paste(phenoName,' ~ age+age2+genotyping.array+sex+age*sex+age2*sex+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
               PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+',
                               envir_name,'*SNP'))
  mod <- lm(mod.formula,data=dataf.20_test)
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

fwrite(df.results.save,paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.',s,'.diet_score.more_snp.txt'),quote = F,sep = '\t',na = 'NA',row.names = F,col.names = T)






