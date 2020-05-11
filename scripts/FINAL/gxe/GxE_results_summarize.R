library(data.table)
pheno <- 'bmi'

#################################################################
# GxE

s='20';results.20 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.',s,'.diet_score.more_snp.txt'),data.table = F,stringsAsFactors = F)
s='80';results.80 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.',s,'.diet_score.more_snp.txt'),data.table = F,stringsAsFactors = F)
results.mg.diet <- merge(results.80,results.20,by=c('SNP','E'))

s='20';results.20 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.',s,'.ext.more_snp.txt'),data.table = F,stringsAsFactors = F)
s='80';results.80 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.',s,'.ext.more_snp.txt'),data.table = F,stringsAsFactors = F)
results.mg.all <- merge(results.80,results.20,by=c('SNP','E'))

# results.mg <- results.mg.all
results.mg <- rbind(results.mg.diet,results.mg.all)
g <- regexpr("_[^_]*$", results.mg$SNP)-1
results.mg$SNP <- substring(results.mg$SNP,1,g)

f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/sig_results/',pheno,'.sig.txt')
df.mg2 <- fread(f,data.table = F,stringsAsFactors = F)
results.mg = merge(results.mg,df.mg2,by.x='SNP',by.y='rs')

f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/sig_results/',pheno,'.sig.GxE.txt')
fwrite(results.mg,f.out,quote = F,na = 'NA',sep = '\t',row.names = F,col.names = T)

library(plyr)
ddply(df.mg2,.(Mean.QTL,Raw.vQTL,Rint.vQTL,dQTL),nrow)
ddply(df.mg2,.(Rint.vQTL,Rint.log_vQTL),nrow)
ddply(df.mg2,.(Rint.vQTL,dQTL),nrow)
ddply(df.mg2,.(Rint.log_vQTL,dQTL),nrow)
ddply(df.mg2,.(Mean.QTL,Raw.vQTL),nrow)
ddply(df.mg2,.(Mean.QTL,Rint.vQTL),nrow)
ddply(df.mg2,.(Mean.QTL,dQTL),nrow)
df.mg2[order(df.mg2$dispersion_pval)[1:5],]

ddply(df.mg2,.(Mean.QTL,Raw.vQTL,Rint.vQTL),nrow)
ddply(df.mg2,.(Mean.QTL,Raw.vQTL,dQTL),nrow)

apply(df.mg2[,c('Mean.QTL','Raw.vQTL','Rint.vQTL','dQTL')],2,sum)

environmental_factors <- c(
                           'DIET_SCORE',
                           'age','Alcohol_intake_frequency',
                           'PA','SB','sex','Smoking.E')
results <- subset(results,E %in% environmental_factors)
subset(results,results[,4] < 0.05 / nrow(results))
subset(results,SNP=='rs4743930')
subset(results,SNP=='rs17451107')
subset(df.mg2,rs=='rs17451107')
##################################################


validation_statistics <- function(results.mg,thres.vec=10^-(seq(0,3,by=0.5)),flip=T,PVAL.THRES=NULL) {
  start <- T
  for (i in 1:length(thres.vec)) {
    thres <- thres.vec[i]
    df.sub <- subset(results.mg,results.mg[,4] <= thres); 
    df.sub.sub <- subset(df.sub,sign(Estimate.x)==sign(Estimate.y))
    if (nrow(df.sub) > 0) {
      if (!is.null(PVAL.THRES)) {
        df.sub.sub <- subset(df.sub.sub,df.sub.sub[,6] <= PVAL.THRES)
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
    } else {
      N=NA
      p=NA
      lower=NA
      upper=NA
      pval=NA
      winners_curse=NA
    }
    
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

# using results.trim
results.mg$FDR <- p.adjust(results.mg[,4],method = 'fdr')
tmp <- subset(results.mg,FDR<0.1); P.FDR_0.1 <- tmp[which.max(tmp$FDR),4]; nrow(tmp)
tmp <- subset(results.mg,FDR<0.05); P.FDR_0.05 <- tmp[which.max(tmp$FDR),4]; nrow(tmp)
tmp <- subset(results.mg,FDR<0.01); P.FDR_0.01 <- tmp[which.max(tmp$FDR),4]; nrow(tmp)

# NEED NULL FOR WINNER'S CURSE
# thres=0.05; suff <- ifelse(is.null(thres),'','PVAL.')
thres=NULL; suff <- ifelse(is.null(thres),'','PVAL.')
# THRESHOLD_VECTOR <- sort(c(1,0.1,0.05,0.01,0.005,0.001,0.05/length(unique(results.mg$SNP)),0.05/nrow(results.mg),P.FDR_0.1,P.FDR_0.05,P.FDR_0.01),decreasing = T)
THRESHOLD_VECTOR <- sort(c(1,0.1,0.05,0.01,0.005,0.001,P.FDR_0.1,P.FDR_0.05,P.FDR_0.01),decreasing = T)
df.save <- validation_statistics(results.mg,PVAL.THRES = thres,thres.vec = THRESHOLD_VECTOR,flip = F)
# df.save <- validation_statistics(results.mg,PVAL.THRES = thres,thres.vec = 6e-5,flip = F) 
df.save.var_raw <- validation_statistics(subset(results.mg,Raw.vQTL==1),PVAL.THRES = thres,thres.vec = THRESHOLD_VECTOR,flip = F)#,thres.vec=10^-(seq(0,1.5,by=0.5)))
df.save.only_mean <- validation_statistics(subset(results.mg,Mean.QTL==1 & Criteria==1),PVAL.THRES = thres,thres.vec = THRESHOLD_VECTOR,flip = F)
df.save.var_rint <- validation_statistics(subset(results.mg,Rint.vQTL==1 & Raw.vQTL==0 & Mean.QTL==0),PVAL.THRES = thres,thres.vec = THRESHOLD_VECTOR,flip = F)
df.save.dispersion <- validation_statistics(subset(results.mg,dQTL==1  & Raw.vQTL==0 & Mean.QTL==0),PVAL.THRES = thres,thres.vec = THRESHOLD_VECTOR,flip = F)
df.save.dispersion <- validation_statistics(subset(results.mg,(dQTL==1 | Rint.vQTL==1)  & Raw.vQTL==0 & Mean.QTL==0),PVAL.THRES = thres,thres.vec = THRESHOLD_VECTOR,flip = F)

df.save.young <- validation_statistics(subset(results.mg,Mean.QTL==1 & dispersion_pval < 0.05),PVAL.THRES = thres,thres.vec = THRESHOLD_VECTOR,flip = F)







df.save <- validation_statistics(results.mg,PVAL.THRES = thres,thres.vec = THRESHOLD_VECTOR,flip = F)
df.save.var_raw <- validation_statistics(subset(results.mg,P.VAR.RAW < 5e-8),PVAL.THRES = thres,thres.vec = THRESHOLD_VECTOR,flip = F)#,thres.vec=10^-(seq(0,1.5,by=0.5)))
df.save.only_mean <- validation_statistics(subset(results.mg,P.MEAN < 5e-8 & P.VAR.RAW > 5e-8 & P.VAR.RINT > 1e-5),PVAL.THRES = thres,thres.vec = THRESHOLD_VECTOR,flip = F)
df.save.mean.plus_var <- validation_statistics(subset(results.mg,P.MEAN < 5e-8 & (P.VAR.RAW < 5e-8 | P.VAR.RINT < 1e-5)),PVAL.THRES = thres,thres.vec = THRESHOLD_VECTOR,flip = F)
df.save.mean <- validation_statistics(subset(results.mg,P.MEAN < 5e-8),PVAL.THRES = thres,thres.vec = THRESHOLD_VECTOR,flip = F)
df.save.strict_mean <- validation_statistics(subset(results.mg,P.MEAN < 1e-15),PVAL.THRES = thres,thres.vec = THRESHOLD_VECTOR,flip = F)
# df.save.var_raw <- validation_statistics(subset(results.mg,P.VAR.RAW < 5e-8),flip = F,thres.vec=10^-(seq(0,2.5,by=0.5)),PVAL.THRES = thres)
df.save.var_raw <- validation_statistics(subset(results.mg,P.VAR.RAW < 5e-8),PVAL.THRES = thres,thres.vec = THRESHOLD_VECTOR,flip = F)#,thres.vec=10^-(seq(0,1.5,by=0.5)))
df.save.var_rint <- validation_statistics(subset(results.mg,P.VAR.RINT < 1e-5),PVAL.THRES = thres,thres.vec = THRESHOLD_VECTOR,flip = F)
df.save.strict_var_rint <- validation_statistics(subset(results.mg,P.VAR.RINT < 5e-8),PVAL.THRES = thres,thres.vec = THRESHOLD_VECTOR,flip = F)
df.save.only_var_rint <- validation_statistics(subset(results.mg,P.VAR.RINT < 1e-5 & P.VAR.RAW > 5e-8 & P.MEAN > 5e-8),PVAL.THRES = thres,thres.vec = THRESHOLD_VECTOR,flip = F)
df.save.var_raw <- validation_statistics(subset(results.mg,P.VAR.RAW < 5e-8),PVAL.THRES = thres,thres.vec = THRESHOLD_VECTOR,flip = F)#,thres.vec=10^-(seq(0,1.5,by=0.5)))
df.save.mean.var_raw_matched_mean <- validation_statistics(subset(results.mg,P.VAR.RAW > 5e-8 & P.MEAN < 5e-17),PVAL.THRES = thres,thres.vec = THRESHOLD_VECTOR,flip = F)#,thres.vec=10^-(seq(0,1.5,by=0.5)))

x <- subset(results.mg,P.VAR.RAW > 5e-8 & P.MEAN < 5e-17); x[order(x[,4])[1:5],]
x <- subset(results.mg,P.VAR.RAW > 5e-8); x[order(x[,4])[1:5],]

nrow(subset(results.mg,P.MEAN < 5e-8 & FDR < 0.1))
nrow(subset(results.mg,P.VAR.RAW < 5e-8 & FDR < 0.1))
nrow(subset(results.mg,P.MEAN < 5e-8 & P.VAR.RAW > 5e-8 & FDR < 0.1))

# median(subset(results.mg,P.VAR.RAW < 5e-8)$P.MEAN)
# median(subset(results.mg,P.VAR.RAW > 5e-8 & P.MEAN < 5e-17)$P.MEAN)

# df.save.var_rint <- validation_statistics(subset(results.mg,P.VAR.RINT < 1e-5),flip = F,thres.vec=10^-(seq(0,2,by=0.5)),PVAL.THRES = thres)

# df.save.var_rint <- validation_statistics(subset(results.mg,P.VAR.RINT < 5e-8)) # rs13198716 interesting. all interactions replicate
# df.save.var_rint <- validation_statistics(subset(results.mg,P.VAR.RINT < 1e-5 & P.VAR.RINT > 5e-8),thres.vec=10^-(seq(0,2,by=0.5)))
df.save.mean.var_raw <- validation_statistics(subset(results.mg,P.MEAN < 5e-8 & P.VAR.RAW < 5e-8),PVAL.THRES = thres,thres.vec = THRESHOLD_VECTOR,flip = F)
df.save.mean.var_rint <- validation_statistics(subset(results.mg,P.MEAN < 5e-8 & P.VAR.RINT < 1e-5),PVAL.THRES = thres,thres.vec = THRESHOLD_VECTOR,flip = F)
df.save.var_raw.var_rint <- validation_statistics(subset(results.mg,P.VAR.RAW < 5e-8 & P.VAR.RINT < 1e-5),PVAL.THRES = thres,thres.vec = THRESHOLD_VECTOR,flip = F)
df.save.mean.var_raw.var_rint <- validation_statistics(subset(results.mg,P.MEAN < 5e-8 & P.VAR.RAW < 5e-8 & P.VAR.RINT < 1e-5),PVAL.THRES = thres,thres.vec = THRESHOLD_VECTOR,flip = F)
df.save.only_mean <- validation_statistics(subset(results.mg,P.MEAN < 5e-8 & P.VAR.RAW > 5e-8 & P.VAR.RINT > 1e-5),PVAL.THRES = thres,thres.vec = THRESHOLD_VECTOR,flip = F)
df.save.mean.plus_var <- validation_statistics(subset(results.mg,P.MEAN < 5e-8 & (P.VAR.RAW < 5e-8 | P.VAR.RINT < 1e-5)),PVAL.THRES = thres,thres.vec = THRESHOLD_VECTOR,flip = F)
df.save.some_var <- validation_statistics(subset(results.mg,(P.VAR.RAW < 5e-8 | P.VAR.RINT < 1e-5)),PVAL.THRES = thres,thres.vec = THRESHOLD_VECTOR,flip = F)

df.save.strict_mean.no_var <- validation_statistics(subset(results.mg,P.MEAN < 1e-15 & P.VAR.RAW > 5e-8 & P.VAR.RINT > 1e-5),PVAL.THRES = thres,thres.vec = THRESHOLD_VECTOR,flip = F)
df.save.strict_mean.plus_var <- validation_statistics(subset(results.mg,P.MEAN < 1e-15 & (P.VAR.RAW < 5e-8 | P.VAR.RINT < 1e-5)),PVAL.THRES = thres,thres.vec = THRESHOLD_VECTOR,flip = F)

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
f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.validation.summary.mean.var_raw_matched_mean.',suff,'txt')
fwrite(df.save.mean.var_raw_matched_mean,f.out,quote = F,na='NA',sep = '\t',row.names = F,col.names = T)



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



