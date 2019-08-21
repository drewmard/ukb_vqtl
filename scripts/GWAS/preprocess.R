df2 <- fread('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/ukb_sample_qc.txt',data.table = F,stringsAsFactors = F)
df.disease <- fread('/athena/elementolab/scratch/anm2868/vQTL/UKB/blood_disease.indiv_id.txt',data.table = F,stringsAsFactors = F,header = T)

pheno <- fread('/home/kulmsc/athena/ukbiobank/phenotypes/ukb26867.csv.gz',data.table=F,stringsAsFactors = F)
pheno2 <- pheno[,c('eid','22001-0.0','21022-0.0','22000-0.0','22007-0.0','22008-0.0','22009-0.1',
                   '21001-0.0')]
colnames(pheno2)[2:ncol(pheno2)] <- c('sex','age','batch','plate','well','p.c.1',
                                      'bmi')
pheno2 <- subset(pheno2, !(eid %in% df.disease$eid)) # remove individuals w/ disease
pheno2$sex[which(is.na(pheno2$sex))] <- median(pheno2$sex,na.rm=T)
pheno2$age[which(is.na(pheno2$age))] <- median(pheno2$age,na.rm=T)
pheno2$age2 <- pheno2$age^2

tmp <- pheno2$bmi
tmp[which(abs(scale(tmp)) > 5)] <- NA
val <- median(tmp,na.rm=T)
pheno2$bmi2 <- pheno2$bmi;
pheno2$bmi2.dummy <- as.numeric(is.na(pheno2$bmi))
i <- which(is.na(pheno2$bmi2))
j <- which(abs(scale(pheno2$bmi)) > 5)
pheno2$bmi2[i] <- val
pheno2$bmi2.with_outliers <- pheno2$bmi2
pheno2$bmi2[j] <- NA

pheno.new <- fread('/home/kulmsc/athena/ukbiobank/setup_morePhenos/ukb33822.csv.gz',data.table=F,stringsAsFactors = F)
pheno.new2 <- pheno.new[,c('eid',paste0('X',c('3700.0.0','2724.0.0','1558.0.0','1239.0.0','1249.0.0')))]
colnames(pheno.new2)[2:ncol(pheno.new2)] <- c('time.since.period','menopause','alcohol.freq','current.smoking','past.smoking')
pheno.new2 <- subset(pheno.new2, !(eid %in% df.disease$eid)) # remove indiv w/ disease

pheno.new2$Smoking <- NA
pheno.new2$Smoking[which( (pheno.new2$current.smoking==1 | pheno.new2$current.smoking==2) | (pheno.new2$past.smoking==1 | pheno.new2$past.smoking==2) )] <- 1
pheno.new2$Smoking[which( pheno.new2$current.smoking==0 & (pheno.new2$past.smoking==3 | pheno.new2$past.smoking==4) )] <- 0
pheno.new2$Smoking.dummy <- as.numeric(is.na(pheno.new2$Smoking))
pheno.new2$Smoking[which(is.na(pheno.new2$Smoking))] <- median(pheno.new2$Smoking,na.rm=T)

tmp <- pheno.new2$time.since.period
tmp[which(abs(scale(tmp)) > 5)] <- NA
val <- median(tmp,na.rm=T)
pheno.new2$time.since.period2 <- pheno.new2$time.since.period;
pheno.new2$time.since.period2[which(pheno.new2$time.since.period2 %in% c(-1,-3))] <- NA
pheno.new2$time.since.period2.dummy <- as.numeric(is.na(pheno.new2$time.since.period2))
i <- which(is.na(pheno.new2$time.since.period2))
j <- which(abs(scale(pheno.new2$time.since.period2)) > 5)
pheno.new2$time.since.period2[i] <- val
pheno.new2$time.since.period2.with_outliers <- pheno.new2$time.since.period2
pheno.new2$time.since.period2[j] <- NA

pheno.new2$alcohol.freq2 <- pheno.new2$alcohol.freq;
pheno.new2$alcohol.freq2[which(pheno.new2$alcohol.freq %in% c(-3))] <- NA
pheno.new2$alcohol.freq2.dummy <- as.numeric(is.na(pheno.new2$alcohol.freq2))
pheno.new2$alcohol.freq2[which(is.na(pheno.new2$alcohol.freq2))] <- median(pheno.new2$alcohol.freq2,na.rm=T)

pheno.new2$menopause2 <- pheno.new2$menopause
pheno.new2$menopause2[which(pheno.new2$menopause==-3)] <- NA
pheno.new2$menopause2[which(pheno.new2$menopause==3)] <- NA
pheno.new2$menopause2[which(is.na(pheno.new2$menopause2))] <- 3
pheno.new2$menopause2 <- as.factor(pheno.new2$menopause2)

# phenotypes
pheno3 <- pheno[,c('eid','30120-0.0','30130-0.0','30140-0.0','30200-0.0','30000-0.0')]
pheno3 <- subset(pheno3, !(eid %in% df.disease$eid))  # remove indiv w/ disease
PHENOTYPE_NAMES <- c('lymphocyte.count','monocyte.count','neutrophil.count','neutrophil.percentage','wbc.leukocyte.count')
colnames(pheno3)[2:ncol(pheno3)] <- PHENOTYPE_NAMES
pheno2 <- merge(pheno2,pheno3,by='eid')
pheno2 <- merge(pheno2,pheno.new2,by='eid')

# Read in fam & merge w/ covariate & phenotype data
fam <- fread('/home/kulmsc/athena/ukbiobank/calls/ukbb.1.fam',data.table = F,stringsAsFactors = F)
fam2 <- cbind(fam,df2)
fam2 <- merge(fam2,pheno2,by.x='V1',by.y='eid',all.x=TRUE)
for (i in 1:length(PHENOTYPE_NAMES)) {
  phenoName <- PHENOTYPE_NAMES[i]
  fam2[,paste0(phenoName,'.na')] <- fam2[,phenoName]
}


# Neale subset:
Neale_subset <- read.table('/athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/samples.both_sexes.tsv.bgz',header=T,stringsAsFactors = F)
fam3 <- merge(fam2,Neale_subset,all.x=TRUE,by.x=c('plate.x','well.x'),by.y=c('plate_name','well'))

Neale_subset$In <- 1
fam3$In[which(is.na(fam3$In))] <- 0
for (i in 1:length(PHENOTYPE_NAMES)) {
  phenoName <- PHENOTYPE_NAMES[i]
  
  # creates British-only phenotype that has NA for non-European individuals
  fam3[,paste0(phenoName,'.na')][which(fam3$In==0)] <- NA
  fam3[,paste0(phenoName,'.na')][which(fam3$QC_In==0)] <- NA
# which(fam3$In==1 & fam3$QC_In==0) # only 1 individual
}
fam2 <- fam3

# phenotypeDataFile <- fam2[,c('V1','V2',paste0(PHENOTYPE_NAMES,'.na'),
#                             'sex','age','age2','genotyping.array',
#                             paste0('PC',1:20))]
phenotypeDataFile <- fam2[,-c(1:2)]
colnames(phenotypeDataFile)[1:2] <- c('FID','IID')

fwrite(phenotypeDataFile,'/athena/elementolab/scratch/anm2868/vQTL/UKB/files/pheno.blood.txt',quote=F,col.names = T,row.names = F,na='NA',sep='\t')

# Covariate model & residuals

# Create phenotype file
phenotypeDataFile <- fam2[,c('V1','V2')]
colnames(phenotypeDataFile) <- c('FID','IID')

fam3 <- fam2
fam3$bmi2 <- fam3$bmi2.with_outliers
fam3$time.since.period2 <- fam3$time.since.period2.with_outliers

############

set.seed(99)
i.20 <- sample(1:nrow(fam2),0.2*nrow(fam2),replace = F)
fam2.80 <- fam2; fam3.80 <- fam3
fam2.20 <- fam2; fam3.20 <- fam3
# 20 percent
fam2.80[i.20,] <- NA
fam2.20[(1:nrow(fam2))[-i.20],] <- NA
fam3.80[i.20,] <- NA
fam3.20[(1:nrow(fam2))[-i.20],] <- NA

# other
for (k in 1:length(PHENOTYPE_NAMES)) {
  
  phenoName <- PHENOTYPE_NAMES[k]
  i.outlier <- which(abs(scale(fam2[,phenoName])) > 5)
  fam2[,paste0(phenoName,'.na')][i.outlier] <- NA
  
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
fwrite(as.data.frame(i.20),'/athena/elementolab/scratch/anm2868/vQTL/UKB/phenotypes_blood.indiv_id.txt',sep='\t',quote = F,col.names = F,row.names = F)


