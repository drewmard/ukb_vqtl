# CAD ascertainment was based on a composite of myocardial infarction or coronary revascularization. 
# Myocardial infarction was based on self-report or hospital admission diagnosis, as performed centrally. 
# This included individuals with 
# ICD-9 codes of 410.X, 411.0, 412.X, or 429.79, or 
# ICD-10 codes of I21.X, I22.X, I23.X, I24.1, or I25.2 in hospitalization records
# Coronary revascularization was assessed based on an OPCS-4 coded procedure for coronary artery bypass grafting 
# (K40.1–40.4, K41.1–41.4, or K45.1–45.5), or 
# coronary angioplasty with or without stenting 
# (K49.1–49.2, K49.8–49.9, K50.2, K75.1–75.4, or K75.8–75.9).

library('data.table')
library('pROC')

read2=TRUE
if (read2==TRUE) {
  pheno <- fread('/home/kulmsc/athena/ukbiobank/phenotypes/ukb26867.csv.gz',data.table=F,stringsAsFactors = F)
}

i1 <- which(pheno[,'20002-0.0'] %in% c(1075))

i2 <- which(pheno[,'41203-0.0'] %in% c('410', '4109','411','4119','412','4129'))
i3 <- which(pheno[,'41205-0.0'] %in% c('410', '4109','411','4119','412','4129'))

i4 <- which(pheno[,'41202-0.0'] %in% c('I21','I210','I211','I212','I213','I214','I219','I21X', 
                                       'I22','I220','I221','I228','I229',
                                       'I23','I230','I231','I232','I233','I234','I235','I236' ,
                                       'I241','I252')) # icd10
i5 <- which(pheno[,'41204-0.0'] %in% c('I21','I210','I211','I212','I213','I214','I219','I21X', 
                                       'I22','I220','I221','I228','I229',
                                       'I23','I230','I231','I232','I233','I234','I235','I236' ,
                                       'I241','I252')) # icd10

i6 <- which(pheno[,'41200-0.0'] %in% c('K40','K401','K402','K403','K404',
                                       'K41','K411','K412','K413','K414',
                                       'K451','K452','K453','K454','K455',
                                       'K491','K492',
                                       'K498','K499', 
                                       'K502', 
                                       'K751','K752','K753','K754', 
                                       'K758','K759'))
i7 <- which(pheno[,'41210-0.0'] %in% c('K40','K401','K402','K403','K404',
                                       'K41','K411','K412','K413','K414',
                                       'K451','K452','K453','K454','K455',
                                       'K491','K492',
                                       'K498','K499', 
                                       'K502', 
                                       'K751','K752','K753','K754', 
                                       'K758','K759'))

x <- unique(c(i1,i2,i3,i4,i5,i6,i7))

CAD <- data.frame(eid=pheno$eid,CAD=as.numeric(c(1:nrow(pheno) %in% x)))
f.out<-'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/CAD.txt'
fwrite(CAD,f.out,quote = F,na='NA',sep = '\t',col.names = T,row.names = F)

Diabetes <- data.frame(eid=pheno$eid,Diabetes=as.numeric(pheno[,'2443-0.0']))
Diabetes$Diabetes[Diabetes$Diabetes < 0] <- NA
f.out<-'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/Diabetes.txt'
fwrite(Diabetes,f.out,quote = F,na='NA',sep = '\t',col.names = T,row.names = F)

HBP <- data.frame(eid=pheno$eid,HBP=as.numeric(pheno[,'6150-0.0']))
HBP$HBP[HBP$HBP < 0] <-0
HBP$HBP[HBP$HBP %in% c(1,2,3)] <- NA
HBP$HBP[HBP$HBP == 4] <- 1
f.out<-'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/HBP.txt'
fwrite(HBP,f.out,quote = F,na='NA',sep = '\t',col.names = T,row.names = F)



