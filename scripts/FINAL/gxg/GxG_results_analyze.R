library(data.table)
phenoName <- 'bmi'

s <- '80'
f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.',phenoName,'.merged_subset2.GxG.',s,'.epi.qt')
df.80 <- fread(f,data.table=F,stringsAsFactors = F)
colnames(df.80)[5:7] <- paste0(colnames(df.80)[5:7],'.',s)
df.80[order(df.80$P.80)[1:5],] 

s <- '20'
f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.',phenoName,'.merged_subset2.GxG.',s,'.epi.qt')
df.20 <- fread(f,data.table = F,stringsAsFactors = F)
colnames(df.20)[5:7] <- paste0(colnames(df.20)[5:7],'.',s)

df <- merge(df.80,df.20,by=c('CHR1','SNP1','CHR2','SNP2'))
df$FDR <- p.adjust(df$P.80,method = 'fdr')

validation_statistics <- function(results.mg,thres.vec=10^-(seq(0,3,by=0.5)),flip=T,PVAL.THRES=NULL) {
  start <- T
  for (i in 1:length(thres.vec)) {
    thres <- thres.vec[i]
    df.sub <- subset(results.mg,P.80 <= thres); 
    df.sub.sub <- subset(df.sub,sign(BETA_INT.80)==sign(BETA_INT.20))
    if (nrow(df.sub) > 0) {
      if (!is.null(PVAL.THRES)) {
        df.sub.sub <- subset(df.sub.sub,P.20 <= PVAL.THRES)
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
      winners_curse <- nrow(subset(df.sub.sub,abs(BETA_INT.80) > abs(BETA_INT.20)))/nrow(df.sub.sub)
    } else {
      N=NA
      p=NA
      lower=NA
      upper=NA
      pval=NA
      winners_curse=NA
    }
    
    df.tmp <- data.frame(thres=thres,
                         N=N,p=p,lower=lower,upper=upper,pval=pval,winner=winners_curse)
    
    if (start) {
      df.save <- df.tmp
      start <- F
    } else {
      df.save <- rbind(df.save,df.tmp)
    }
  }
  return(df.save)
}

thres=0.05
df.save <- validation_statistics(df,thres.vec = 10^-(seq(0,5,by=0.5)),flip=F,PVAL.THRES=thres)
# validation_statistics(subset(df,FDR < .9),thres.vec = 10^-(seq(0,5,by=0.5)),flip=F,PVAL.THRES=thres)

cor.test(df$BETA_INT.80,df$BETA_INT.20)
cor.test(subset(df,P.80<1e-3)$BETA_INT.80,subset(df,P.80<1e-3)$BETA_INT.20)

f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.',phenoName,'.merged_subset2.GxG.FULL.txt')
fwrite(df,f,quote = F,na = 'NA',sep = '\t',row.names = F,col.names = T)
f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.',phenoName,'.merged_subset2.GxG.FULL.valid.txt')
fwrite(df.save,f,quote = F,na = 'NA',sep = '\t',row.names = F,col.names = T)



# f.freq <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.',phenoName,'.merged_subset.GxG.frq')
# # f.freq <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.lymphocyte.count.merged_subset.GxG.frq'
# df.freq <- fread(f.freq,data.table = F,stringsAsFactors = F)
# df2 <- merge(df,df.freq[,c('SNP','MAF')],by.x='SNP1',by.y='SNP')
# df2 <- merge(df2,df.freq[,c('SNP','MAF')],by.x='SNP2',by.y='SNP')
# colnames(df2)[which(colnames(df2) %in% c('MAF.x','MAF.y'))] <- c('MAF1','MAF2')
# df2 <- subset(df2,MAF1 > 0.05 & MAF2 > 0.05)
# df2$Sign.Validate <- as.numeric(sign(df2$BETA_INT.80) == sign(df2$BETA_INT.20))
# 
# f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.',phenoName,'.merged_subset.GxG.','ALL_TRIMMED','.epi.qt')
# fwrite(df2,f.out,quote = F,na='NA',sep = '\t',row.names = F,col.names = T)
# 
