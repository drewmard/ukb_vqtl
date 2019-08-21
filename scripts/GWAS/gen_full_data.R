library(data.table)

# sample_qc.R
df2 <- fread('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/ukb_sample_qc.txt',data.table = F,stringsAsFactors = F)

# gen_cov1.R
pheno2 <- fread('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/covariates1.txt',data.table = F,stringsAsFactors = F)

# gen_cov2.R
pheno.new2 <- fread('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/covariates2.txt',data.table = F,stringsAsFactors = F)
pheno.new2$menopause2 <- as.factor(pheno.new2$menopause2)

# gen_pheno.R
pheno3 <- fread('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/phenotypes.txt',data.table = F,stringsAsFactors = F)

# identify_indiv_blood_disorders.R
# df.disease <- fread('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/blood_disease.indiv_id.txt',data.table = F,stringsAsFactors = F,header = T)

# merge files
pheno2 <- merge(pheno2,pheno3,by='eid')
pheno2 <- merge(pheno2,pheno.new2,by='eid')

# Read in fam & merge w/ covariate & phenotype data
fam <- fread('/home/kulmsc/athena/ukbiobank/calls/ukbb.1.fam',data.table = F,stringsAsFactors = F)
fam2 <- as.data.frame(cbind(fam,df2))
fam2 <- merge(x=fam2,y=pheno2,by.x='V1',by.y='eid',all.x=TRUE)
# for (i in 1:length(PHENOTYPE_NAMES)) {
#   phenoName <- PHENOTYPE_NAMES[i]
# }

# Neale subset:
Neale_subset <- read.table('/athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/samples.both_sexes.tsv.bgz',header=T,stringsAsFactors = F)
Neale_subset$In <- 1
fam3 <- merge(fam2,Neale_subset,all.x=TRUE,by.x=c('plate.x','well.x'),by.y=c('plate_name','well'))
fam3$In[which(is.na(fam3$In))] <- 0
for (i in 1:length(PHENOTYPE_NAMES)) {
  phenoName <- PHENOTYPE_NAMES[i]
  
  # creates British-only phenotype that has NA for non-European individuals by using Neale QC
  fam3[,paste0(phenoName,'.na')] <- fam3[,phenoName]
  fam3[,paste0(phenoName,'.na')][which(fam3$In==0)] <- NA
  fam3[,paste0(phenoName,'.na')][which(fam3$QC_In==0)] <- NA
  # which(fam3$In==1 & fam3$QC_In==0) # only 1 individual
}
fam2 <- fam3

phenotypeDataFile <- fam2[,-c(1:2)] # removes plate.x and well.x; still have *.y
colnames(phenotypeDataFile)[1:2] <- c('FID','IID')

fwrite(phenotypeDataFile,'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/full_data.txt',quote=F,col.names = T,row.names = F,na='NA',sep='\t')
