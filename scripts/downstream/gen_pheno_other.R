library(data.table)

pheno <- fread('/home/kulmsc/athena/ukbiobank/phenotypes/ukb26867.csv.gz',data.table=F,stringsAsFactors = F)
df.disease <- fread('/athena/elementolab/scratch/anm2868/vQTL/UKB/blood_disease.indiv_id.txt',data.table = F,stringsAsFactors = F,header = T)
# phenotypes
pheno3 <- pheno[,c('eid','21022-0.0','20002-0.0','30120-0.0')]
colnames(pheno3)[2] <- 'age'
colnames(pheno3)[4] <- 'lymphocyte.count'
pheno3$psoriasis <- as.numeric(pheno3[,'20002-0.0']==1453)
pheno3$psoriasis[which(is.na(pheno3$psoriasis))] <- 0
pheno3$coeliac_disease <- as.numeric(pheno3[,'20002-0.0']==1456)
pheno3$coeliac_disease[which(is.na(pheno3$coeliac_disease))] <- 0
pheno3$hyperthyroidism <- as.numeric(pheno3[,'20002-0.0']==1225)
pheno3$hyperthyroidism[which(is.na(pheno3$hyperthyroidism))] <- 0
# pheno3 <- subset(pheno3, !(eid %in% df.disease$eid))  # remove indiv w/ disease
pheno3 <- pheno3[,-(which(colnames(pheno3) %in% c('20002-0.0','age')))]

f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/results/full_data_gxe.',s,'.txt')
df2 <- fread(f,data.table = F,stringsAsFactors = F)
df2$menopause2 <- as.factor(df2$menopause2)
fam3 <- df2
fam3$bmi2 <- df2$bmi2.with_outliers
fam3$time.since.period2 <- df2$time.since.period2.with_outliers
fam3 <- merge(fam3,pheno3,by.x='IID',by.y='eid')

DISEASE='psoriasis'
# vQTL='rs887468'
vQTL='rs3132506'
mod.formula.1 <- paste0(DISEASE,' ~ genotyping.array+
       PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
       PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+
       menopause2+
       bmi2.dummy+bmi2+
    Smoking+Smoking.dummy+
    alcohol.freq2+alcohol.freq2.dummy')

mod.formula.1.ext <- paste0(mod.formula.1,'+age*',vQTL)
mod.formula.1.ext <- formula(mod.formula.1.ext)
mod <- glm(mod.formula.1.ext,data=fam3,family=binomial(link='logit'))
summary(mod)


# mod <- glm(psoriasis~age*rs887468,data=fam3,family=binomial(link='logit'))
# summary(mod)

# fwrite(pheno3,'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/phenotypes.txt',na='NA',row.names = F,col.names = T,quote = F,sep = '\t')


user_direc <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl'
phenoName <- 'lymphocyte.count.rint'
f.geno <- paste0(user_direc,'/output/GWAS/subset/',phenoName,'/ukbb.ALL_vQTL.raw')
df.geno <- fread(f.geno,data.table = F,stringsAsFactors = F)

df <- merge(df.geno,pheno3,by.x='IID',by.y='eid')


snp <- 'rs377763'
summary(lm(psoriasis~rs887468+age,data=df))
summary(lm(psoriasis~rs887468*age,data=df)) # significant
summary(lm(coeliac_disease~rs887468+age,data=df)) 
summary(lm(coeliac_disease~rs887468*age,data=df)) # loses lots of significance?
summary(lm(hyperthyroidism~rs887468*age,data=df))
summary(lm(lymphocyte.count~rs887468*age,data=df))
library(car)
leveneTest(lymphocyte.count~as.factor(rs887468),data=df)
# case control imbalance, maybe not the best data set to try this out on?

summary(lm(psoriasis~rs887468*age,data=df))


