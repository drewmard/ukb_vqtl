library(data.table)

pheno <- fread('/home/kulmsc/athena/ukbiobank/phenotypes/ukb26867.csv.gz',data.table=F,stringsAsFactors = F)
df.disease <- fread('/athena/elementolab/scratch/anm2868/vQTL/UKB/blood_disease.indiv_id.txt',data.table = F,stringsAsFactors = F,header = T)

# phenotypes
pheno3 <- pheno[,c('eid','21022-0.0','20002-0.0')]
colnames(pheno3)[2] <- 'age'
pheno3$psoriasis <- as.numeric(pheno3[,'20002-0.0']==1453)
pheno3$psoriasis[which(is.na(pheno3$psoriasis))] <- 0
pheno3$coeliac_disease <- as.numeric(pheno3[,'20002-0.0']==1456)
pheno3$coeliac_disease[which(is.na(pheno3$coeliac_disease))] <- 0
pheno3$hyperthyroidism <- as.numeric(pheno3[,'20002-0.0']==1225)
pheno3$hyperthyroidism[which(is.na(pheno3$hyperthyroidism))] <- 0


pheno3 <- subset(pheno3, !(eid %in% df.disease$eid))  # remove indiv w/ disease

pheno3 <- pheno3[,-(which(colnames(pheno3) %in% '20002-0.0'))]

# fwrite(pheno3,'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/phenotypes.txt',na='NA',row.names = F,col.names = T,quote = F,sep = '\t')


user_direc <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl'
phenoName <- 'lymphocyte.count.rint'
f.geno <- paste0(user_direc,'/output/GWAS/subset/',phenoName,'/ukbb.ALL_vQTL.raw')
df.geno <- fread(f.geno,data.table = F,stringsAsFactors = F)

df <- merge(df.geno,pheno3,by.x='IID',by.y='eid')

summary(lm(psoriasis~rs887468+age,data=df))
summary(lm(psoriasis~rs887468*age,data=df)) # significant
summary(lm(coeliac_disease~rs887468+age,data=df)) 
summary(lm(coeliac_disease~rs887468*age,data=df)) # loses lots of significance?
summary(lm(hyperthyroidism~rs887468*age,data=df))
# case control imbalance, maybe not the best data set to try this out on?

summary(lm(psoriasis~rs887468*age,data=df))


