library(data.table)

# read in data
pheno.new <- fread('/home/kulmsc/athena/ukbiobank/setup_morePhenos/ukb33822.csv.gz',data.table=F,stringsAsFactors = F)
pheno.new2 <- pheno.new[,c('eid',c('3700-0.0','2724-0.0','1558-0.0','1239-0.0','1249-0.0'))]
colnames(pheno.new2)[2:ncol(pheno.new2)] <- c('time.since.period','menopause','alcohol.freq','current.smoking','past.smoking')

# remove individuals w/ disease
df.disease <- fread('/athena/elementolab/scratch/anm2868/vQTL/UKB/blood_disease.indiv_id.txt',data.table = F,stringsAsFactors = F,header = T)
pheno.new2 <- subset(pheno.new2, !(eid %in% df.disease$eid)) # remove indiv w/ disease

# split 80/20 train/test
df.20 <- fread('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/phenotypes_blood.indiv_id.txt',data.table = F,stringsAsFactors = F)
pheno.new2.80 <- subset(pheno.new2,!(eid %in% df.20$IID))
pheno.new2.20 <- subset(pheno.new2,(eid %in% df.20$IID))

for (iter in 1:2) {
  
  if (iter==1) {pheno.new2 <- pheno.new2.80}
  if (iter==2) {pheno.new2 <- pheno.new2.20}
  
  # smoking
  pheno.new2$Smoking <- NA
  pheno.new2$Smoking[which( (pheno.new2$current.smoking==1 | pheno.new2$current.smoking==2) | (pheno.new2$past.smoking==1 | pheno.new2$past.smoking==2) )] <- 1
  pheno.new2$Smoking[which( pheno.new2$current.smoking==0 & (pheno.new2$past.smoking==3 | pheno.new2$past.smoking==4) )] <- 0
  pheno.new2$Smoking.dummy <- as.numeric(is.na(pheno.new2$Smoking))
  pheno.new2$Smoking[which(is.na(pheno.new2$Smoking))] <- median(pheno.new2$Smoking,na.rm=T)
  
  # time since period
  pheno.new2$time.since.period2 <- pheno.new2$time.since.period;
  pheno.new2$time.since.period2[which(pheno.new2$time.since.period2 %in% c(-1,-3))] <- NA # missing
  pheno.new2$time.since.period2.dummy <- as.numeric(is.na(pheno.new2$time.since.period2)) # dummy
  i <- which(is.na(pheno.new2$time.since.period2)) # missing
  j <- which(abs(scale(pheno.new2$time.since.period2)) > 5) # outliers
  tmp <- pheno.new2$time.since.period
  tmp[which(abs(scale(tmp)) > 5)] <- NA
  val <- median(tmp,na.rm=T) # median value w/o outliers
  pheno.new2$time.since.period2[i] <- val
  pheno.new2$time.since.period2.with_outliers <- pheno.new2$time.since.period2
  pheno.new2$time.since.period2[j] <- NA # remove outliers
  
  # alcohol freq
  pheno.new2$alcohol.freq2 <- pheno.new2$alcohol.freq;
  pheno.new2$alcohol.freq2[which(pheno.new2$alcohol.freq %in% c(-3))] <- NA
  pheno.new2$alcohol.freq2.dummy <- as.numeric(is.na(pheno.new2$alcohol.freq2))
  pheno.new2$alcohol.freq2[which(is.na(pheno.new2$alcohol.freq2))] <- median(pheno.new2$alcohol.freq2,na.rm=T)
  
  # menopause
  pheno.new2$menopause2 <- pheno.new2$menopause
  pheno.new2$menopause2[which(pheno.new2$menopause==-3)] <- NA
  pheno.new2$menopause2[which(pheno.new2$menopause==3)] <- NA
  pheno.new2$menopause2[which(is.na(pheno.new2$menopause2))] <- 3
  
  if (iter == 1) {
    f <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/covariates2.80.txt'
  } else if (iter==2) {
    f <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/covariates2.20.txt'
  }
  fwrite(pheno.new2,f,na='NA',row.names = F,col.names = T,quote = F,sep = '\t')
}