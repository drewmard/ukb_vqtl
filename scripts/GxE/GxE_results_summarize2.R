library(data.table)
pheno <- 'bmi'

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

s='20';results.20 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.',s,'.diet_score.more_snp.txt'),data.table = F,stringsAsFactors = F)
s='80';results.80 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.',s,'.diet_score.more_snp.txt'),data.table = F,stringsAsFactors = F)
results.mg.diet <- merge(results.80,results.20,by=c('SNP','E'))

s='20';results.20 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.',s,'.ext.more_snp.txt'),data.table = F,stringsAsFactors = F)
s='80';results.80 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.',s,'.ext.more_snp.txt'),data.table = F,stringsAsFactors = F)
results.mg.all <- merge(results.80,results.20,by=c('SNP','E'))

results.mg <- rbind(results.mg.diet,results.mg.all)
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
results.ld_sub <- subset(results,!(SNP %in% snp_to_remove))
results.mg <- results.ld_sub

subset(results.mg,results.mg[,4] < 0.05 / nrow(results.mg))

validation_statistics <- function(results.mg) {
  thres.vec <- 10^-(seq(0,-log10(min(results.mg[,4]))+0.1,by=0.1))
  start <- T
  for (i in 1:length(thres.vec)) {
    thres <- thres.vec[i]
    
    df.sub <- subset(results.mg,results.mg[,4] < thres); 
    same_sign_prop <- nrow(df.sub.sub <- subset(df.sub,sign(Estimate.x)==sign(Estimate.y)))/nrow(df.sub); 
    winners_curse <- nrow(subset(df.sub.sub,abs(Estimate.x) > abs(Estimate.y)))/nrow(df.sub.sub)
    n=nrow(df.sub)
    
    df.sub <- subset(results.mg,results.mg[,4] > thres); 
    same_sign_prop2 <- nrow(df.sub.sub <- subset(df.sub,sign(Estimate.x)==sign(Estimate.y)))/nrow(df.sub); 
    winners_curse2 <- nrow(subset(df.sub.sub,abs(Estimate.x) > abs(Estimate.y)))/nrow(df.sub.sub)
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
df.save <- validation_statistics(results.mg)
df.save.mean <- validation_statistics(subset(results.mg,P.MEAN < 5e-8))
df.save.strict_mean <- validation_statistics(subset(results.mg,P.MEAN < 1e-15))
df.save.var_raw <- validation_statistics(subset(results.mg,P.VAR.RAW < 5e-8))
df.save.var_rint <- validation_statistics(subset(results.mg,P.VAR.RINT < 1e-5))
df.save.mean.var_raw <- validation_statistics(subset(results.mg,P.MEAN < 5e-8 & P.VAR.RAW < 5e-8))
df.save.mean.var_rint <- validation_statistics(subset(results.mg,P.MEAN < 5e-8 & P.VAR.RINT < 1e-5))
df.save.var_raw.var_rint <- validation_statistics(subset(results.mg,P.VAR.RAW < 5e-8 & P.VAR.RINT < 1e-5))
df.save.mean.var_raw.var_rint <- validation_statistics(subset(results.mg,P.MEAN < 5e-8 & P.VAR.RAW < 5e-8 & P.VAR.RINT < 1e-5))

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

f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.validation.summary.txt')
fwrite(df.save.full,f.out,quote = F,na='NA',sep = '\t',row.names = F,col.names = T)

f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.','ALL','.ext.more_snp.txt')
fwrite(results.mg,f.out,quote = F,na='NA',sep = '\t',row.names = F,col.names = T)



