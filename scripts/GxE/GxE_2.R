library(data.table)
library(sandwich)
library(lmtest)
library(parallel)


# initialize
user_direc <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl'
phenoName <- 'lymphocyte.count.rint.ALL'

# var hits
f <- paste0(user_direc,'/output/GWAS/results2/',phenoName,'/ukbb.',phenoName,'.results.var.clumped.cut.txt')
index.SNP <- fread(f,data.table = F,stringsAsFactors = F)

# mean hits
# f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/UKB/results/ukbb.',phenoName,'.results.mean.clumped.cut.txt')
# index <- fread(f,data.table = F,stringsAsFactors = F)

# index <- fread('/athena/elementolab/scratch/anm2868/vQTL/UKB/results/ukbb.bmi.results.custom.txt',data.table = F,stringsAsFactors = F)
print('Creating genotype file...')
for (i in 1:nrow(index.SNP)) {
  # for (i in 1:10) {
  print(paste0(i,'/',nrow(index.SNP),' SNPs...'))
  vQTL=index.SNP[i,2]
  CHR_vQTL=index.SNP[i,1]
  tryCatch(
    {
      f.vQTL <- paste0(user_direc,'/output/GWAS/subset/',phenoName,'/ukbb.',CHR_vQTL,'.',vQTL,'.raw')
      df.vQTL <- fread(f.vQTL,data.table = F,stringsAsFactors = F)
      
      df.vQTL <- df.vQTL[,c(2,grep(vQTL,colnames(df.vQTL)))]
      colnames(df.vQTL)[2] <- vQTL
      
      if (i != 1) {
        df.geno <- merge(df.geno,df.vQTL,by='IID')
      } else {
        df.geno <- df.vQTL
      }
    },error=function(e) {
      print(paste0('ERROR: ',i))
      break
    }
  )
  
}

f.out <- paste0(user_direc,'/output/GWAS/subset/',phenoName,'/ukbb.ALL_vQTL.raw')
fwrite(df.geno,f.out,sep='\t',quote=F,row.names = F,col.names = T,na='NA')

source('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/scripts/GxE/gen_datafiles.R')

# initialize
# phenoName <- 'lymphocyte.count.rint'
phenoName <- 'lymphocyte.count.rint.ALL'
# phenoName <- 'lymphocyte.count.ALL'
phenoName2 <- 'lymphocyte.count'

PHENOTYPE_NAMES <- phenoName

# PHENOTYPE_NAMES <- c('lymphocyte.count','monocyte.count','neutrophil.count','neutrophil.percentage','wbc.leukocyte.count')
# ENVIR_NAMES <- c("PA","SB","Smoking.E","sleep.duration","getting.up.morning","nap.during.day",
#                  "time.spent.outdoors.summer","time.spent.outdoors.winter","time.spent.outdoors",
#                  "tobacco.smoke.exposure","alcohol.freq.E","childhood.sunburn.occasions",
#                  "age.started.smoking","hormone.replacement.therapy",
#                  'age')#,
# 'sex')
ENVIR_NAMES <- c("PA","SB","Smoking.E","sleep.duration","alcohol.freq.E",
                 'age')


#,'bmi2','menopause','time.since.period2')
# ENVIR_NAMES <- c('PA','SB','Smoking.E','alcohol.freq.E','childhood.sunburn.occasions','age','sex') #,'bmi2','menopause','time.since.period2')

# loading:
source('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/scripts/GxE/gen_datafiles.R');
gen_datafiles(phenoName,phenoName2)
# PHENOTYPE_NAMES <- paste0(rep(PHENOTYPE_NAMES,each=3),c('','.log','.rint'),'.na')
# phenoName2 <- paste0(phenoName,'.na')


s <- '20'

f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/results/full_data_gxe.',s,'.txt')
df2 <- fread(f,data.table = F,stringsAsFactors = F)
df2$menopause2 <- as.factor(df2$menopause2)

index <- colnames(df2)[2:(which(colnames(df2)=='FID')-1)]

fam3 <- df2
fam3$bmi2 <- df2$bmi2.with_outliers
fam3$time.since.period2 <- df2$time.since.period2.with_outliers

# source('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/scripts/GxE/GxE_acrossPhenotypes_2.R')
source('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/scripts/GxE/GxE_acrossEnvirFactor.R')

print('GxE...')
START=TRUE
y2 <- mclapply(1:length(index),GxE_acrossSNPs,mc.cores = 16)
y2.df <- as.data.frame(do.call(rbind,y2))
for (l in 4:length(y2.df)) {
  y2.df[,l] <- as.numeric(levels(y2.df[,l]))[y2.df[,l]]
}

results <- y2.df
print(results[order(results$P_GxE),][1:10,])

# SAVE # # # #
f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/results/ukbb.gxe.',phenoName,'.',s,'.txt')
fwrite(results,f.out,sep='\t',na='NA',row.names = F,col.names = T,quote = F)
