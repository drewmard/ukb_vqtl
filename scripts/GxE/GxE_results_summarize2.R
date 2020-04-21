library(data.table)
pheno <- 'bmi'
# pheno <- 'lymphocyte.count'

f1 <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_subset/ukbb.',paste0(pheno,'.ALL'),'.vGWAS.txt')
f2 <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_subset/ukbb.',paste0(pheno,'.rint.ALL'),'.vGWAS.txt')
f3 <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/imputed/results/ukbb.',paste0(pheno,'.ALL'),'.results.txt')

var.raw <- fread(f1,data.table = F,stringsAsFactors = F)
var.rint <- fread(f2,data.table = F,stringsAsFactors = F)
mean.raw <- fread(f3,data.table = F,stringsAsFactors = F); colnames(mean.raw)[1] <- 'rs'

f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.',pheno,'.merged_subset.bim')
bim <- fread(f,data.table = F,stringsAsFactors = F)
var.raw.sub <- subset(var.raw,rs %in% bim[,2])
var.rint.sub <- subset(var.rint,rs %in% bim[,2])
mean.raw.sub <- subset(mean.raw,rs %in% bim[,2])

mean.raw.sub <- mean.raw.sub[,c('rs','CHR','BP','A1','A2','MAF','BETA','P')]
colnames(mean.raw.sub)[(ncol(mean.raw.sub)-1):ncol(mean.raw.sub)] <- paste0(colnames(mean.raw.sub)[(ncol(mean.raw.sub)-1):ncol(mean.raw.sub)],'.MEAN')
var.raw.sub <- var.raw.sub[,c('rs','BETA','P')]; colnames(var.raw.sub)[2:3] <- paste0(colnames(var.raw.sub)[2:3],'.VAR.RAW')
var.rint.sub <- var.rint.sub[,c('rs','BETA','P')]; colnames(var.rint.sub)[2:3] <- paste0(colnames(var.rint.sub)[2:3],'.VAR.RINT')

df.mg <- merge(merge(mean.raw.sub,var.raw.sub,by='rs'),var.rint.sub,by='rs')

#################################################################

f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.',pheno,'.merged_subset.LD.ld')
ld <- fread(f,data.table = F,stringsAsFactors = F)
ld <- subset(ld,R2 > 0.1)
df.mg.ld <- subset(df.mg, rs %in% c(ld$SNP_A,ld$SNP_B))
df.mg.ld <- df.mg.ld[order(df.mg.ld$CHR),]

df.mg2 <- df.mg[,c('rs','CHR','BP','P.MEAN','P.VAR.RAW','P.VAR.RINT')]
#################################################################
# GxE

s='20';results.20 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.',s,'.diet_score.more_snp.txt'),data.table = F,stringsAsFactors = F)
s='80';results.80 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.',s,'.diet_score.more_snp.txt'),data.table = F,stringsAsFactors = F)
results.mg.diet <- merge(results.80,results.20,by=c('SNP','E'))

s='20';results.20 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.',s,'.ext.more_snp.txt'),data.table = F,stringsAsFactors = F)
s='80';results.80 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.',s,'.ext.more_snp.txt'),data.table = F,stringsAsFactors = F)
results.mg.all <- merge(results.80,results.20,by=c('SNP','E'))

results.mg <- results.mg.all
# results.mg <- rbind(results.mg.diet,results.mg.all)
g <- regexpr("_[^_]*$", results.mg$SNP)-1
results.mg$SNP <- substring(results.mg$SNP,1,g)

results = merge(results.mg,df.mg2,by.x='SNP',by.y='rs')

environmental_factors <- c(
                           'DIET_SCORE',
                           'age','Alcohol_intake_frequency',
                           'PA','SB','sex','Smoking.E')
results <- subset(results,E %in% environmental_factors)
subset(results,results[,4] < 0.05 / nrow(results))


f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.',pheno,'.merged_subset.LD.ld')
ld <- fread(f,data.table = F,stringsAsFactors = F)
ld <- subset(ld,R2 > 0.1)
results.aggre <-  aggregate(results[,'Pr(>|t|).x'],by=list(results$SNP),min)
results.ld <- subset(results.aggre, results.aggre[,1] %in% c(ld$SNP_A,ld$SNP_B))
results.ld <- results.ld[order(results.ld[,2],decreasing = F),]

i=1
while (i <= nrow(results.ld)) {
  snp = results.ld[i,1]
  tmp <- subset(ld,SNP_A==snp | SNP_B==snp)
  tmp <- subset(c(tmp$SNP_A,tmp$SNP_B),c(tmp$SNP_A,tmp$SNP_B)!=snp)
  results.ld <- subset(results.ld,!(results.ld[,1] %in% c(tmp)))
  i=i+1
}

snp_to_remove <- subset(c(ld$SNP_A,ld$SNP_B), !(c(ld$SNP_A,ld$SNP_B) %in% results.ld[,1]))
f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/',pheno,'.snp_to_remove.txt')
fwrite(data.frame(snp_to_remove=snp_to_remove),f.out,quote = F,sep = '\t',na = 'NA',row.names = F,col.names = T)
snp_to_remove <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/',pheno,'.snp_to_remove.txt'),data.table = F,stringsAsFactors = F)
results.ld_sub <- subset(results,!(SNP %in% snp_to_remove[,1]))
results.mg <- results.ld_sub

subset(results.mg,results.mg[,4] < 0.05 / nrow(results.mg))
results.mg[order(results.mg[,4],decreasing = F)[1:5],]
results.trim <- subset(results.mg,E %in% environmental_factors)
subset(results.trim,results.trim[,4] < 0.05 / nrow(results.trim))
f.out<-paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.ALL_RESULTS.trim.txt')
fwrite(results.trim,f.out,quote = F,sep = '\t',na = 'NA',row.names = F,col.names = T)

validation_statistics <- function(results.mg,thres.vec=10^-(seq(0,3,by=0.5)),flip=T,PVAL.THRES=NULL) {
  start <- T
  for (i in 1:length(thres.vec)) {
    thres <- thres.vec[i]
    df.sub <- subset(results.mg,results.mg[,4] < thres); 
    df.sub.sub <- subset(df.sub,sign(Estimate.x)==sign(Estimate.y))
    if (!is.null(PVAL.THRES)) {
      df.sub.sub <- subset(df.sub.sub,df.sub.sub[,6] < PVAL.THRES)
      X=nrow(df.sub.sub); N=nrow(df.sub)
      binomial.test.results <- binom.test(x=X,n=N,PVAL.THRES / 2)
    } else {
      X=nrow(df.sub.sub); N=nrow(df.sub)
      binomial.test.results <- binom.test(x=X,n=N,0.5)
    }
    p=as.numeric(binomial.test.results$estimate)
    ci=as.numeric(binomial.test.results$conf.int)
    lower=ci[1]; upper=ci[2]
    pval=binomial.test.results$p.value
    winners_curse <- nrow(subset(df.sub.sub,abs(Estimate.x) > abs(Estimate.y)))/nrow(df.sub.sub)
    
    # se=sqrt(p*(1-p)/(nrow(df.sub)))
    # p+1.96*se
    # same_sign_prop <- nrow(df.sub.sub <- subset(df.sub,sign(Estimate.x)==sign(Estimate.y)))/nrow(df.sub); 
    # winners_curse <- nrow(subset(df.sub.sub,abs(Estimate.x) > abs(Estimate.y)))/nrow(df.sub.sub)
    # n=nrow(df.sub)
    
    if (flip) {
      df.sub <- subset(results.mg,results.mg[,6] < thres); 
      df.sub.sub <- subset(df.sub,sign(Estimate.x)==sign(Estimate.y))
      if (!is.null(PVAL.THRES)) {
        df.sub.sub <- subset(df.sub.sub,df.sub.sub[,4] < PVAL.THRES)
        X2=nrow(df.sub.sub); N2=nrow(df.sub)
        binomial.test.results <- binom.test(x=X2,n=N2,PVAL.THRES / 2)
      } else {
        X2=nrow(df.sub.sub); N2=nrow(df.sub)
        binomial.test.results <- binom.test(x=X2,n=N2,0.5)
      }
      p2=as.numeric(binomial.test.results$estimate)
      ci=as.numeric(binomial.test.results$conf.int)
      lower2=ci[1]; upper2=ci[2]
      pval2=binomial.test.results$p.value
      winners_curse2 <- nrow(subset(df.sub.sub,abs(Estimate.x) < abs(Estimate.y)))/nrow(df.sub.sub)
    } else {
      N2=NA
      p2=NA
      lower2=NA
      upper2=NA
      pval2=NA
      winners_curse2=NA
    }
    
    
    # df.sub <- subset(results.mg,results.mg[,4] > thres); 
    # same_sign_prop2 <- nrow(df.sub.sub <- subset(df.sub,sign(Estimate.x)==sign(Estimate.y)))/nrow(df.sub); 
    # winners_curse2 <- nrow(subset(df.sub.sub,abs(Estimate.x) > abs(Estimate.y)))/nrow(df.sub.sub)
    # n2=nrow(df.sub)
    
    df.tmp <- data.frame(thres=thres,
                         N=N,p=p,lower=lower,upper=upper,pval=pval,winner=winners_curse,
                         N2=N2,p2=p2,lower2=lower2,upper2=upper2,pval2=pval2,winner2=winners_curse2)
    if (start) {
      df.save <- df.tmp
      start <- F
    } else {
      df.save <- rbind(df.save,df.tmp)
    }
  }
  return(df.save)
}

read_data=FALSE
if (read_data) {
  library(data.table)
  pheno='bmi'
  f<-paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.ALL_RESULTS.trim.txt')
  results.mg <- fread(f,data.table = F,stringsAsFactors = F)
  f <- "/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.80.var.raw_matched_snp.merged_subset.txt"
  results.80 <- fread(f,data.table = F,stringsAsFactors = F)
  f <- "/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.20.var.raw_matched_snp.merged_subset.txt"
  results.20 <- fread(f,data.table = F,stringsAsFactors = F)
  f <- "/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.80.mean.raw_matched_snp.merged_subset.txt"
  results.80 <- fread(f,data.table = F,stringsAsFactors = F)
  f <- "/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.20.mean.raw_matched_snp.merged_subset.txt"
  results.20 <- fread(f,data.table = F,stringsAsFactors = F)
  results.mg <- merge(results.80,results.20,by=c('SNP','E'))
}

thres=NULL; suff <- ifelse(is.null(thres),'','PVAL.')
df.save <- validation_statistics(results.mg,PVAL.THRES = thres)
f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.validation.mean.raw_matched_snp.full.',suff,'txt')
fwrite(df.save,f.out,quote = F,na='NA',sep = '\t',row.names = F,col.names = T)
thres=0.05; suff <- ifelse(is.null(thres),'','PVAL.')
df.save <- validation_statistics(results.mg,PVAL.THRES = thres)
f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.validation.mean.raw_matched_snp.full.',suff,'txt')
fwrite(df.save,f.out,quote = F,na='NA',sep = '\t',row.names = F,col.names = T)

# using results.trim
thres=0.05; suff <- ifelse(is.null(thres),'','PVAL.')
df.save <- validation_statistics(results.mg,PVAL.THRES = thres)
df.save.var_raw <- validation_statistics(subset(results.mg,P.VAR.RAW < 5e-8),PVAL.THRES = thres)#,thres.vec=10^-(seq(0,1.5,by=0.5)))
df.save.only_mean <- validation_statistics(subset(results.mg,P.MEAN < 5e-8 & P.VAR.RAW > 5e-8 & P.VAR.RINT > 1e-5),PVAL.THRES = thres)
df.save.mean.plus_var <- validation_statistics(subset(results.mg,P.MEAN < 5e-8 & (P.VAR.RAW < 5e-8 | P.VAR.RINT < 1e-5)),PVAL.THRES = thres)
df.save.mean <- validation_statistics(subset(results.mg,P.MEAN < 5e-8),PVAL.THRES = thres)
df.save.strict_mean <- validation_statistics(subset(results.mg,P.MEAN < 1e-15),PVAL.THRES = thres)
# df.save.var_raw <- validation_statistics(subset(results.mg,P.VAR.RAW < 5e-8),flip = F,thres.vec=10^-(seq(0,2.5,by=0.5)),PVAL.THRES = thres)
df.save.var_raw <- validation_statistics(subset(results.mg,P.VAR.RAW < 5e-8),PVAL.THRES = thres)#,thres.vec=10^-(seq(0,1.5,by=0.5)))
df.save.var_rint <- validation_statistics(subset(results.mg,P.VAR.RINT < 1e-5),PVAL.THRES = thres)
# df.save.var_rint <- validation_statistics(subset(results.mg,P.VAR.RINT < 1e-5),flip = F,thres.vec=10^-(seq(0,2,by=0.5)),PVAL.THRES = thres)

# df.save.var_rint <- validation_statistics(subset(results.mg,P.VAR.RINT < 5e-8)) # rs13198716 interesting. all interactions replicate
# df.save.var_rint <- validation_statistics(subset(results.mg,P.VAR.RINT < 1e-5 & P.VAR.RINT > 5e-8),thres.vec=10^-(seq(0,2,by=0.5)))
df.save.mean.var_raw <- validation_statistics(subset(results.mg,P.MEAN < 5e-8 & P.VAR.RAW < 5e-8),PVAL.THRES = thres)
df.save.mean.var_rint <- validation_statistics(subset(results.mg,P.MEAN < 5e-8 & P.VAR.RINT < 1e-5),PVAL.THRES = thres)
df.save.var_raw.var_rint <- validation_statistics(subset(results.mg,P.VAR.RAW < 5e-8 & P.VAR.RINT < 1e-5),PVAL.THRES = thres)
df.save.mean.var_raw.var_rint <- validation_statistics(subset(results.mg,P.MEAN < 5e-8 & P.VAR.RAW < 5e-8 & P.VAR.RINT < 1e-5),PVAL.THRES = thres)
df.save.only_mean <- validation_statistics(subset(results.mg,P.MEAN < 5e-8 & P.VAR.RAW > 5e-8 & P.VAR.RINT > 1e-5),PVAL.THRES = thres)
df.save.mean.plus_var <- validation_statistics(subset(results.mg,P.MEAN < 5e-8 & (P.VAR.RAW < 5e-8 | P.VAR.RINT < 1e-5)),PVAL.THRES = thres)
df.save.some_var <- validation_statistics(subset(results.mg,(P.VAR.RAW < 5e-8 | P.VAR.RINT < 1e-5)),PVAL.THRES = thres)

df.save.strict_mean.no_var <- validation_statistics(subset(results.mg,P.MEAN < 1e-15 & P.VAR.RAW > 5e-8 & P.VAR.RINT > 1e-5),PVAL.THRES = thres)
df.save.strict_mean.plus_var <- validation_statistics(subset(results.mg,P.MEAN < 1e-15 & (P.VAR.RAW < 5e-8 | P.VAR.RINT < 1e-5)),PVAL.THRES = thres)

f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.validation.summary.full.',suff,'txt')
fwrite(df.save,f.out,quote = F,na='NA',sep = '\t',row.names = F,col.names = T)
f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.validation.summary.var_raw.',suff,'txt')
fwrite(df.save.var_raw,f.out,quote = F,na='NA',sep = '\t',row.names = F,col.names = T)
f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.validation.summary.var_rint.',suff,'txt')
fwrite(df.save.var_rint,f.out,quote = F,na='NA',sep = '\t',row.names = F,col.names = T)
f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.validation.summary.mean.var_raw.',suff,'txt')
fwrite(df.save.mean.var_raw,f.out,quote = F,na='NA',sep = '\t',row.names = F,col.names = T)
f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.validation.summary.mean.var_rint.',suff,'txt')
fwrite(df.save.mean.var_rint,f.out,quote = F,na='NA',sep = '\t',row.names = F,col.names = T)
f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.validation.summary.var_raw.var_rint.',suff,'txt')
fwrite(df.save.var_raw.var_rint,f.out,quote = F,na='NA',sep = '\t',row.names = F,col.names = T)
f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.validation.summary.only_mean.',suff,'txt')
fwrite(df.save.only_mean,f.out,quote = F,na='NA',sep = '\t',row.names = F,col.names = T)
f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.validation.summary.mean.plus_var.',suff,'txt')
fwrite(df.save.mean.plus_var,f.out,quote = F,na='NA',sep = '\t',row.names = F,col.names = T)
f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.validation.summary.some_var.',suff,'txt')
fwrite(df.save.some_var,f.out,quote = F,na='NA',sep = '\t',row.names = F,col.names = T)
f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.validation.summary.strict_mean.no_var.',suff,'txt')
fwrite(df.save.strict_mean.no_var,f.out,quote = F,na='NA',sep = '\t',row.names = F,col.names = T)
f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.validation.summary.strict_mean.plus_var.',suff,'txt')
fwrite(df.save.strict_mean.plus_var,f.out,quote = F,na='NA',sep = '\t',row.names = F,col.names = T)



df.save.full <- data.frame(
  thres=df.save$thres,
  p.all=df.save$p,
  p.mean=df.save.mean$p,
  p.strict_mean=df.save.strict_mean$p,
  p.var_raw=df.save.var_raw$p,
  p.var_rint=df.save.var_rint$p,
  p.mean.var_raw=df.save.mean.var_raw$p,
  p.mean.var_rint=df.save.mean.var_rint$p,
  p.var_raw.var_rint=df.save.var_raw.var_rint$p,
  p.mean.var_raw.var_rint=df.save.mean.var_raw.var_rint$p
)

f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.validation.summary.',suff,'txt')
fwrite(df.save.full,f.out,quote = F,na='NA',sep = '\t',row.names = F,col.names = T)

# f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.','ALL','.ext.more_snp.',suff,'txt')
# fwrite(results.mg,f.out,quote = F,na='NA',sep = '\t',row.names = F,col.names = T)



