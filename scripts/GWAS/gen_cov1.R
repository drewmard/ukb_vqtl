library(data.table)
pheno <- fread('/home/kulmsc/athena/ukbiobank/phenotypes/ukb26867.csv.gz',data.table=F,stringsAsFactors = F)

# extract certain data
pheno2 <- pheno[,c('eid','22001-0.0','21022-0.0','22000-0.0','22007-0.0','22008-0.0','22009-0.1',
                   '21001-0.0')]
colnames(pheno2)[2:ncol(pheno2)] <- c('sex','age','batch','plate','well','p.c.1',
                                      'bmi')

# remove individuals w/ disease
df.disease <- fread('/athena/elementolab/scratch/anm2868/vQTL/UKB/blood_disease.indiv_id.txt',data.table = F,stringsAsFactors = F,header = T)
pheno2 <- subset(pheno2, !(eid %in% df.disease$eid)) 

# split 80/20 train/test
df.20 <- fread('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/phenotypes_blood.indiv_id.txt',data.table = F,stringsAsFactors = F)
pheno2.80 <- subset(pheno2,!(eid %in% df.20$IID))
pheno2.20 <- subset(pheno2,(eid %in% df.20$IID))

for (iter in 1:2) {
  
  if (iter==1) {pheno2 <- pheno2.80}
  if (iter==2) {pheno2 <- pheno2.20}
  
  # impute sex & age
  pheno2$sex[which(is.na(pheno2$sex))] <- median(pheno2$sex,na.rm=T)
  pheno2$age[which(is.na(pheno2$age))] <- median(pheno2$age,na.rm=T)

  # squared age
  pheno2$age2 <- pheno2$age^2
  
  # goal: remove outliers & impute median
  pheno2$bmi2 <- pheno2$bmi; # new variable: bmi2
  pheno2$bmi2.dummy <- as.numeric(is.na(pheno2$bmi)) #create new dummy variable 
  i <- which(is.na(pheno2$bmi2)) # i = indiv w/ missing data
  j <- which(abs(scale(pheno2$bmi)) > 5) # j = indiv w/ outliers
  tmp <- pheno2$bmi
  tmp[which(abs(scale(tmp)) > 5)] <- NA # make outliers equal to NA
  val <- median(tmp,na.rm=T) # median value
  pheno2$bmi2[i] <- val # assign median value to indiv w/ missing data
  pheno2$bmi2.with_outliers <- pheno2$bmi2
  pheno2$bmi2[j] <- NA
  # bmi: original
  # bmi2: imputed values for indiv missing, and outliers removed after
  # bmi2.dummy: 0, 1 for whether indiv missing
  # bmi2.with_outliers: bmi2 but outliers not removed
  
  if (iter == 1) {
    f <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/covariates1.80.txt'
  } else if (iter==2) {
    f <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/covariates1.20.txt'
  }
  fwrite(pheno2,f,na='NA',row.names = F,col.names = T,quote = F,sep = '\t')
  
}


