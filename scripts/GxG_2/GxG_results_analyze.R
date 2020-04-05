library(data.table)
# phenoName <- 'lymphocyte.count'
# phenoName <- 'monocyte.count.rint.ALL'
# phenoName <- 'neutrophil.count.rint.ALL'
# phenoName <- 'wbc.leukocyte.count.rint.ALL'
# phenoName <- 'rbc.erythrocyte.count.rint.ALL'
phenoName <- 'bmi'
# phenoName <- 'lymphocyte.count'

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

s <- '20'
f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.',phenoName,'.merged_subset.GxG.',s,'.epi.qt')
df.20 <- fread(f,data.table = F,stringsAsFactors = F)
# df.20 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.',phenoName,'.ALL.sub.GxG.',s,'.epi.qt'),data.table=F,stringsAsFactors = F)
colnames(df.20)[5:7] <- paste0(colnames(df.20)[5:7],'.',s)

df <- merge(df.80,df.20,by=c('CHR1','SNP1','CHR2','SNP2'))

f.freq <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.',phenoName,'.merged_subset.GxG.frq')
# f.freq <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.lymphocyte.count.merged_subset.GxG.frq'
df.freq <- fread(f.freq,data.table = F,stringsAsFactors = F)
df2 <- merge(df,df.freq[,c('SNP','MAF')],by.x='SNP1',by.y='SNP')
df2 <- merge(df2,df.freq[,c('SNP','MAF')],by.x='SNP2',by.y='SNP')
colnames(df2)[which(colnames(df2) %in% c('MAF.x','MAF.y'))] <- c('MAF1','MAF2')

df.sub <- subset(df2,CHR1!=CHR2 & MAF1 > 0.05 & MAF2 > 0.05)
df.sub$Sign.Validate <- as.numeric(sign(df.sub$BETA_INT.80) == sign(df.sub$BETA_INT.20))
df.sub[order(df.sub$P.80)[1:5],] 
.05/nrow(df.sub)
mean(subset(df.sub,P.20 < 1e-3)$Sign.Validate)

# results.mg <- df.sub
gxg_validation_statistics <- function(results.mg,thres.vec=(10^-(seq(0,3,by=0.5))),flip=TRUE) {
  start <- T
  for (i in 1:length(thres.vec)) {
    thres <- thres.vec[i]
    
    # df.sub <- subset(results.mg,P.20 < thres); 
    # same_sign_prop <- nrow(df.sub.sub <- subset(df.sub,sign(BETA_INT.80)==sign(BETA_INT.20)))/nrow(df.sub); 
    # winners_curse <- nrow(subset(df.sub.sub,abs(BETA_INT.80) < abs(BETA_INT.20)))/nrow(df.sub.sub)
    # n=nrow(df.sub)
    # 
    # df.sub <- subset(results.mg,P.80 < thres); 
    # same_sign_prop2 <- nrow(df.sub.sub <- subset(df.sub,sign(BETA_INT.80)==sign(BETA_INT.20)))/nrow(df.sub); 
    # winners_curse2 <- nrow(subset(df.sub.sub,abs(BETA_INT.80) > abs(BETA_INT.20)))/nrow(df.sub.sub)
    # n2=nrow(df.sub)
    
    df.sub <- subset(results.mg,P.80 < thres); 
    df.sub.sub <- subset(df.sub,sign(BETA_INT.80)==sign(BETA_INT.20))
    X=nrow(df.sub.sub); N=nrow(df.sub)
    binomial.test.results <- binom.test(x=X,n=N,0.5)
    p=as.numeric(binomial.test.results$estimate)
    ci=as.numeric(binomial.test.results$conf.int)
    lower=ci[1]; upper=ci[2]
    pval=binomial.test.results$p.value
    winners_curse <- nrow(subset(df.sub.sub,abs(BETA_INT.80) > abs(BETA_INT.20)))/nrow(df.sub.sub)
    
    if (flip) {
      df.sub <- subset(results.mg,P.20 < thres); 
      df.sub.sub <- subset(df.sub,sign(BETA_INT.80)==sign(BETA_INT.20))
      X=nrow(df.sub.sub); N2=nrow(df.sub)
      binomial.test.results <- binom.test(x=X,n=N2,0.5)
      p2=as.numeric(binomial.test.results$estimate)
      ci=as.numeric(binomial.test.results$conf.int)
      lower2=ci[1]; upper2=ci[2]
      pval2=binomial.test.results$p.value
      winners_curse2 <- nrow(subset(df.sub.sub,abs(BETA_INT.80) < abs(BETA_INT.20)))/nrow(df.sub.sub)
    } else {
      N2=NA
      p2=NA
      lower2=NA
      upper2=NA
      pval2=NA
      winners_curse2=NA
    }

    df.tmp <- data.frame(thres=thres,
                         N=N,p=p,lower=lower,upper=upper,pval=pval,winner=winners_curse,
                         N2=N2,p2=p2,lower2=lower2,upper2=upper2,pval2=pval2,winner2=winners_curse2)

    # df.tmp <- data.frame(thres=thres,n=n,
    #                      p=same_sign_prop,winner=winners_curse,
    #                      n2=n2,
    #                      p2=same_sign_prop2,winner2=winners_curse2)
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

results <- df.sub
f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.',phenoName,'.merged_subset.LD.ld')
ld <- fread(f,data.table = F,stringsAsFactors = F)
ld <- subset(ld,R2 > 0.1)
results.tmp <- results
results.tmp$SNP1 <- results$SNP2; results.tmp$SNP2 <- results$SNP1
results <- rbind(results,results.tmp)
results.aggre <-  aggregate(results[,'P.80'],by=list(results$SNP1),min)
results.ld <- subset(results.aggre, results.aggre[,1] %in% c(ld$SNP_A,ld$SNP_B))
results.ld <- results.ld[order(results.ld[,2],decreasing = F),]

# do instead based on gxe
# i=1
# while (i <= nrow(results.ld)) {
#   snp = results.ld[i,1]
#   tmp <- subset(ld,SNP_A==snp | SNP_B==snp)
#   tmp <- subset(c(tmp$SNP_A,tmp$SNP_B),c(tmp$SNP_A,tmp$SNP_B)!=snp)
#   results.ld <- subset(results.ld,!(results.ld[,1] %in% c(tmp)))
#   i=i+1
# }
# snp_to_remove <- subset(c(ld$SNP_A,ld$SNP_B), !(c(ld$SNP_A,ld$SNP_B) %in% results.ld[,1]))
# f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/',phenoName,'.snp_to_remove.txt')
# fwrite(data.frame(snp_to_remove=snp_to_remove),f.out,quote = F,sep = '\t',na = 'NA',row.names = F,col.names = T)

f.snp_to_remove <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/bmi.snp_to_remove.txt'
snp_to_remove <- fread(f.snp_to_remove,data.table = F,stringsAsFactors = F)

results.mg.2 <- subset(results,!((SNP1 %in% snp_to_remove[,1]) | (SNP2 %in% snp_to_remove[,1])))
results.mg.2 <- results.mg.2[1:(nrow(results.mg.2)/2),]
df.save <- gxg_validation_statistics(results.mg.2)
f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/',phenoName,'.GxG.validation.summary.txt')
fwrite(df.save,f.out,quote = F,na='NA',sep = '\t',row.names = F,col.names = T)

# working off of GxE_results_summarize2.R
pheno <- phenoName
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



######3
results.mg.3 <- merge(results.mg.2,df.mg2,by.x='SNP1',by.y='rs')
colnames(results.mg.3)[14:18] <- paste0(colnames(results.mg.3)[14:18],'.1')
results.mg.3 <- merge(results.mg.3,df.mg2,by.x='SNP2',by.y='rs')
colnames(results.mg.3)[19:23] <- paste0(colnames(results.mg.3)[19:23],'.2')

# df.save <- gxg_validation_statistics(results.mg.3)
# df.save.mean <- gxg_validation_statistics(subset(results.mg.3,P.MEAN.1 < 5e-8 | P.MEAN.2 < 5e-8))
gxg_validation_statistics.df <- gxg_validation_statistics(results.mg.3)
gxg_validation_statistics.df.var_raw <- gxg_validation_statistics(subset(results.mg.3,P.VAR.RAW.2 < 5e-8 | P.VAR.RAW.1 < 5e-8))
gxg_validation_statistics.df.var_rint <- gxg_validation_statistics(subset(results.mg.3,P.VAR.RINT.2 < 1e-5 | P.VAR.RINT.1 < 1e-5))
gxg_validation_statistics.df.mean <- gxg_validation_statistics(subset(results.mg.3,P.MEAN.1 < 5e-8 | P.MEAN.2 < 5e-8))
f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/',pheno,'.GxG.validation.summary.full.txt')
fwrite(gxg_validation_statistics.df,f.out,quote = F,na='NA',sep = '\t',row.names = F,col.names = T)
f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/',pheno,'.GxG.validation.summary.var_raw.txt')
fwrite(gxg_validation_statistics.df.var_raw,f.out,quote = F,na='NA',sep = '\t',row.names = F,col.names = T)
f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/',pheno,'.GxG.validation.summary.var_rint.txt')
fwrite(gxg_validation_statistics.df.var_rint,f.out,quote = F,na='NA',sep = '\t',row.names = F,col.names = T)
f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/',pheno,'.GxG.validation.summary.mean.txt')
fwrite(gxg_validation_statistics.df.mean,f.out,quote = F,na='NA',sep = '\t',row.names = F,col.names = T)

gxg_validation_statistics(subset(results.mg.3,P.VAR.RINT.2 < 1e-5 & P.VAR.RINT.1 < 1e-5),flip=F,thres.vec=(10^-(seq(0,2,by=0.5))))


gxg_validation_statistics(subset(results.mg.3,P.VAR.RINT.2 < 5e-8 | P.VAR.RINT.1 < 5e-8),flip=F)

gxg_validation_statistics(subset(results.mg.3,(P.MEAN.1 < 5e-8 & P.VAR.RAW.2 < 5e-8) | (P.VAR.RAW.1 < 5e-8 & P.MEAN.2 < 5e-8)))



# subsetting specific gxg:
# seems bull shit and I am skeptical
results.mg.3[order(results.mg.3$P.MEAN.2,decreasing = F),][1:5,]
gxg_validation_statistics(subset(results.mg.3,SNP1=='rs56094641' | SNP2=='rs56094641'),thres = (10^-(seq(0,3,by=0.5))),flip=F)
gxg_validation_statistics(subset(results.mg.3,SNP1=='rs35608615' | SNP2=='rs35608615'),thres = (10^-(seq(0,3,by=0.5))),flip=F)
gxg_validation_statistics(subset(results.mg.3,SNP1=='rs2256752' | SNP2=='rs2256752'),thres = (10^-(seq(0,3,by=0.5))),flip=F)

gxg_validation_statistics(subset(results.mg.3,SNP1=='rs35608615' | SNP2=='rs35608615'),thres = (10^-(seq(0,2.5,by=0.5))))
gxg_validation_statistics(subset(results.mg.3,SNP1=='rs2256752' | SNP2=='rs2256752'),thres = (10^-(seq(0,2.5,by=0.5))))
# 
# SNP_list <- unique(c(results.mg.3$SNP1,results.mg.3$SNP2))
# save <- T
# for (i in 1:length(SNP_list)) {
#   SNP=SNP_list[i]
#   df.tmp <- gxg_validation_statistics(subset(results.mg.3,SNP1==SNP | SNP2==SNP),thres = 0.1)
#   df.tmp$SNP <- SNP
#   if (save) {
#     df.save <- df.tmp
#     save=F
#   } else {
#     df.save <- rbind(df.save,df.tmp)
#   }
# } 
# 
# df.save.mg <- merge(df.save,df.mg2,by.x='SNP',by.y='rs')
# df.save.mg[order(df.save.mg$pval,decreasing = F),][1:5,]







