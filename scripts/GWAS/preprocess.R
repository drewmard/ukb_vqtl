library(data.table)


# read in data
fam2 <- fread('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/full_data.txt',data.table = F,stringsAsFactors = F)

# Create phenotype file
phenotypeDataFile <- fam2[,c('FID','IID')]
colnames(phenotypeDataFile) <- c('FID','IID')

#initialize
PHENOTYPE_NAMES <- c('lymphocyte.count','monocyte.count','neutrophil.count','neutrophil.percentage','wbc.leukocyte.count')

############

# split in 80% train, 20% test
# first: need to only measure outliers in the 80%
fam2$bmi2 <- fam2$bmi2.with_outliers
fam2$time.since.period2 <- fam2$time.since.period2.with_outliers
fam3 <- fam2

set.seed(99)
i.20 <- sample(1:nrow(fam2),0.2*nrow(fam2),replace = F)
fam2.80 <- fam2; 
fam2.20 <- fam2; 
fam3.80 <- fam3
fam3.20 <- fam3
# 20 percent
fam2.80[i.20,] <- NA
fam2.20[(1:nrow(fam2))[-i.20],] <- NA
fam3.80[i.20,] <- NA
fam3.20[(1:nrow(fam2))[-i.20],] <- NA
fwrite(data.frame(row=i.20,IID=fam2$IID[i.20]),'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/phenotypes_blood.indiv_id.txt',sep='\t',quote = F,col.names = F,row.names = F)


#######################

# other
for (k in 1:length(PHENOTYPE_NAMES)) {
  
  phenoName <- PHENOTYPE_NAMES[k]
  i.outlier <- which(abs(scale(fam2.80[,paste0(phenoName,'.na')])) > 5)
  fam2.80[,paste0(phenoName,'.na')][i.outlier] <- NA
  
  i <- which(fam2$sex==1); j <- which(fam2$sex==0); 
  mod.formula.2 <- formula(paste(paste0(phenoName,'.na'),' ~ age+age2+genotyping.array+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
                 PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+
                               Smoking+Smoking.dummy+time.since.period2+time.since.period2.dummy+
                               menopause2+alcohol.freq2+alcohol.freq2.dummy+bmi2.dummy+bmi2+bmi2*age'))
  mod.formula.1 <- formula(paste(paste0(phenoName,'.na'),' ~ age+age2+genotyping.array+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
                                 PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+
                                 Smoking+Smoking.dummy+
                                 alcohol.freq2+alcohol.freq2.dummy+bmi2.dummy+bmi2+bmi2*age'))
  mod1 <- lm(mod.formula.1,
             data=fam2[i,],na.action=na.exclude)
  mod2 <- lm(mod.formula.2,
             data=fam2[j,],na.action=na.exclude)
  mod1 <- lm(mod.formula.1,
             data=fam2.80[i,],na.action=na.exclude)
  mod2 <- lm(mod.formula.2,
             data=fam2[j,],na.action=na.exclude)
  
  resid1 <- residuals(mod1)
  resid2 <- residuals(mod2)
  resid1 <- scale(resid1)
  resid2 <- scale(resid2)
  
  # Supplement to phenotype file
  phenotypeDataFile[,phenoName] <- NA
  phenotypeDataFile[,phenoName][i] <- resid1
  phenotypeDataFile[,phenoName][j] <- resid2
  
  # switch NA to -9
  phenotypeDataFile[,phenoName][which(is.na(phenotypeDataFile[,phenoName]))] <- -9
  
}

fwrite(phenotypeDataFile,'/athena/elementolab/scratch/anm2868/vQTL/UKB/phenotypes_blood.txt',sep='\t',quote = F,col.names = T,row.names = F)


