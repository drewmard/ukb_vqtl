library(data.table)
s='80'

f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/full_data.',s,'.txt')
fam <- fread(f,data.table = F,stringsAsFactors = F)

full_dataset <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/full_data.bmi.GxE.',s,'.txt'),data.table = F,stringsAsFactors = F)
colnames(full_dataset) <- gsub("-| |/",'_',colnames(full_dataset))

pheno <- 'bmi'
phenoName <- paste0(pheno,'.na')

fam[,paste0(pheno,'.na')] <- fam[,pheno]
fam[,paste0(pheno,'.na')][which(fam$In==0)] <- NA
fam[,paste0(pheno,'.na')][which(fam$QC_In==0)] <- NA
i.outlier <- which(abs(scale(fam[,phenoName])) > 6)
fam[,paste0(phenoName,'.outlier')] <- fam[,phenoName]
fam[i.outlier,paste0(phenoName,'.outlier')] <- NA
phenoName <- paste0(phenoName,'.outlier')
full_dataset <- merge(full_dataset,fam[,c('IID',phenoName)],by='IID')

# genetic data
f.geno <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.',pheno,'.merged_subset')
geno <- BEDMatrix(f.geno)
geno_names <- unlist(lapply(strsplit(rownames(geno),'_'),function(x) {return(x[2])}))

# merge stuff
ind <- which(geno_names %in% full_dataset$IID)
full_dataset <- full_dataset[match(geno_names[ind],full_dataset$IID),]


# 
# # snp
i <- which(colnames(geno) %in% "rs11642015_T")
envir_name <- "PA"
# envir_name <- "Alcohol_intake_frequency"
full_dataset$SNP <- geno[ind,i]
mod.formula <- formula(paste(phenoName,' ~ age+age2+genotyping.array+sex+age*sex+age2*sex+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
               PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+',
                             envir_name,'*SNP'))
mod <- lm(mod.formula,data=full_dataset)
res <- summary(mod)$coef[paste0(envir_name,':SNP'),c("Estimate","Pr(>|t|)")]
res







# # environmental factors
environmental_factors <- c(paste0('DIET_PC',1:10),
                           # 'DIET_SCORE',
                           'age','Alcohol_intake_frequency',
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

for (k in 1:length(environmental_factors)) {
  envir_name <- environmental_factors[k]
  df.results <- Fit_Model(60,60)
  # df.results <- Fit_Model(60,60)
  df.results$E <- envir_name
  if (k==1) {
    df.results.save <- df.results
  } else {
    df.results.save <- rbind(df.results.save,df.results)
  }
}

df.results.save.remove_outliers <- df.results.save

df.results.save.mg <- cbind(df.results.save.remove_outliers,df.results.save)

df.results.save.mg$ratio <- -log10(df.results.save.remove_outliers[,3])/-log10(df.results.save[,3])
cbind(df.results.save.remove_outliers[,3],df.results.save[,3])

df.results.save.mg[,c(8,3,7,9)]

