library(data.table)
# phenoName <- 'lymphocyte.count'
# phenoName <- 'monocyte.count.rint.ALL'
# phenoName <- 'neutrophil.count.rint.ALL'
# phenoName <- 'wbc.leukocyte.count.rint.ALL'
# phenoName <- 'rbc.erythrocyte.count.rint.ALL'
phenoName <- 'bmi'

s <- '80'
f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.',phenoName,'.merged_subset.GxG.',s,'.epi.qt')
# f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.',phenoName,'.ALL.sub.GxG.',s,'.epi.qt')
# f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_GxG/ukbb.',phenoName,'.ALL.sub.GxG.epi.qt')
df.80 <- fread(f,data.table=F,stringsAsFactors = F)
colnames(df.80)[5:7] <- paste0(colnames(df.80)[5:7],'.',s)
df.80[order(df.80$P.80)[1:5],] 
# unique(subset(df.80,P.80 < 0.05/nrow(df.80)/5)[,'CHR1'])
# monocyte: 4
# lymphocyte: 28
# neutrophil: 1

subset(df.80,P.80 < 5e-8)

s <- '20'
f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.',phenoName,'.merged_subset.GxG.',s,'.epi.qt')
df.20 <- fread(f,data.table = F,stringsAsFactors = F)
# df.20 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.',phenoName,'.ALL.sub.GxG.',s,'.epi.qt'),data.table=F,stringsAsFactors = F)
colnames(df.20)[5:7] <- paste0(colnames(df.20)[5:7],'.',s)

df <- merge(df.80,df.20,by=c('CHR1','SNP1','CHR2','SNP2'))
df.sub <- df
df.sub[order(df.sub$P.80)[1:5],] 

f.freq <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.bmi.merged_subset.GxG.frq'
df.freq <- fread(f.freq,data.table = F,stringsAsFactors = F)
df2 <- merge(df,df.freq[,c('SNP','MAF')],by.x='SNP1',by.y='SNP')
df2 <- merge(df2,df.freq[,c('SNP','MAF')],by.x='SNP2',by.y='SNP')
colnames(df2)[which(colnames(df2) %in% c('MAF.x','MAF.y'))] <- c('MAF1','MAF2')

df.sub <- subset(df2,CHR1!=CHR2 & MAF1 > 0.05 & MAF2 > 0.05)
df.sub$Sign.Validate <- as.numeric(sign(df.sub$BETA_INT.80) == sign(df.sub$BETA_INT.20))
mean(subset(df.sub,P.20 < 1e-3)$Sign.Validate)

# results.mg <- df.sub
gxg_validation_statistics <- function(results.mg) {
  thres.vec <- 10^-(seq(0,3,by=0.1))
  start <- T
  for (i in 1:length(thres.vec)) {
    thres <- thres.vec[i]
    
    df.sub <- subset(results.mg,P.20 < thres); 
    same_sign_prop <- nrow(df.sub.sub <- subset(df.sub,sign(BETA_INT.80)==sign(BETA_INT.20)))/nrow(df.sub); 
    winners_curse <- nrow(subset(df.sub.sub,abs(BETA_INT.80) < abs(BETA_INT.20)))/nrow(df.sub.sub)
    n=nrow(df.sub)
    
    df.sub <- subset(results.mg,P.80 < thres); 
    same_sign_prop2 <- nrow(df.sub.sub <- subset(df.sub,sign(BETA_INT.80)==sign(BETA_INT.20)))/nrow(df.sub); 
    winners_curse2 <- nrow(subset(df.sub.sub,abs(BETA_INT.80) > abs(BETA_INT.20)))/nrow(df.sub.sub)
    n2=nrow(df.sub)
    
    df.tmp <- data.frame(thres=thres,n=n,
                         p=same_sign_prop,winner=winners_curse,
                         n2=n2,
                         p2=same_sign_prop2,winner2=winners_curse2)
    if (start) {
      df.save <- df.tmp
      start <- F
    } else {
      df.save <- rbind(df.save,df.tmp)
    }
  }
  return(df.save)
}

# df.save <- gxg_validation_statistics(results.mg)
# f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/',pheno,'.GxG.validation.summary.txt')
# fwrite(df.save,f.out,quote = F,na='NA',sep = '\t',row.names = F,col.names = T)

f.snp_to_remove <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/bmi.snp_to_remove.txt'
snp_to_remove <- fread(f.snp_to_remove,data.table = F,stringsAsFactors = F)

results.mg.2 <- subset(results.mg,!(SNP1 %in% snp_to_remove[,1]) & !(SNP2 %in% snp_to_remove[,1]))
df.save <- gxg_validation_statistics(results.mg)
f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/',pheno,'.GxG.validation.summary.txt')
fwrite(df.save,f.out,quote = F,na='NA',sep = '\t',row.names = F,col.names = T)

# working off of GxE_results_summarize2.R
results.mg.3 <- merge(results.mg.2,df.mg2,by.x='SNP1',by.y='rs')
colnames(results.mg.3)[14:18] <- paste0(colnames(results.mg.3)[14:18],'.1')
results.mg.3 <- merge(results.mg.3,df.mg2,by.x='SNP2',by.y='rs')
colnames(results.mg.3)[19:23] <- paste0(colnames(results.mg.3)[19:23],'.2')

df.save <- gxg_validation_statistics(results.mg)
df.save.mean <- gxg_validation_statistics(subset(results.mg.3,P.MEAN.1 < 5e-8 | P.MEAN.2 < 5e-8))
df.save.var.raw <- gxg_validation_statistics(subset(results.mg.3,P.VAR.RAW.1 < 5e-8 | P.VAR.RAW.2 < 5e-8))
df.save.var.rint <- gxg_validation_statistics(subset(results.mg.3,P.VAR.RINT.1 < 1e-5 | P.VAR.RINT.2 < 1e-5))








# 
# df.sub <- subset(df2,CHR1!=CHR2 & MAF1 > 0.2 & MAF2 > 0.2 & P.80 < 0.05 & (sign(BETA_INT.80) == sign(BETA_INT.20)))
# df.sub[order(df.sub$P.80)[1:5],] 
# 
# 
# df <- df.80
# # df.sub <- subset(df,CHR1!=6)
# df.sub <- subset(df,CHR1!=CHR2)
# df.sub[order(df.sub$P.80)[1:5],] 
# 
# df.sub <- subset(df,P.80 < 1e-3 & P.20 < 0.05 & CHR1 != CHR2 & sign(BETA_INT.80) == sign(BETA_INT.20))
# df.sub[order(df.sub$P.80)[1:5],] 
# 
# df.sub <- subset(df,P.80 < 1e-3 & P.20 < 0.05 & CHR1 == CHR2 & sign(BETA_INT.80) == sign(BETA_INT.20))
# df.sub[order(df.sub$P.80)[1:5],] 
# 
# 
# df.sub <- df
# # df.sub <- subset(df,SNP2=='rs146125856' & CHR1!=6)
# # df.sub <- subset(df,CHR1!=CHR2)
# # df.sub <- subset(df,CHR1!=CHR2 & SNP2!='rs146125856')
# # df.sub <- subset(df,CHR1!=CHR2 & CHR1!=6 & CHR2!=6)
# df.sub[order(df.sub$P)[1:5],] 
# 
# df.sub <- subset(df,CHR1!=CHR2)
# df.sub[order(df.sub$P)[1:5],] 
# 
# df.sub <- subset(df,CHR1==15 | CHR2==15)
