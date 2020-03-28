library(data.table)

# original
s='20';results.20 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.',s,'.diet_score.txt'),data.table = F,stringsAsFactors = F)
s='80';results.80 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.',s,'.diet_score.txt'),data.table = F,stringsAsFactors = F)
results.mg.diet <- merge(results.80,results.20,by=c('SNP','E'))

s='20';results.20 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.',s,'.ext.txt'),data.table = F,stringsAsFactors = F)
s='80';results.80 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.',s,'.ext.txt'),data.table = F,stringsAsFactors = F)
results.mg.all <- merge(results.80,results.20,by=c('SNP','E'))

# fake
# s='80';results.20 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.',s,'.diet_score.txt'),data.table = F,stringsAsFactors = F)
# s='20';results.80 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.',s,'.diet_score.txt'),data.table = F,stringsAsFactors = F)
# results.mg.diet <- merge(results.80,results.20,by=c('SNP','E'))
# 
# s='80';results.20 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.',s,'.ext.txt'),data.table = F,stringsAsFactors = F)
# s='20';results.80 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.',s,'.ext.txt'),data.table = F,stringsAsFactors = F)
# results.mg.all <- merge(results.80,results.20,by=c('SNP','E'))

#############

results.mg <- rbind(results.mg.diet,results.mg.all)
results.mg[order(results.mg[,4],decreasing = F),][1:5,]
subset(results.mg,results.mg[,4] < 0.05 / nrow(results.mg))

# aggregate(results.mg[,4],by=list(results.mg$E),min)

environmental_factors <- c(paste0('DIET_PC',1:10),
                           'DIET_SCORE',
                           'age','Alcohol_intake_frequency',
                           'PA','SB','sex','Smoking.E')
results.mg <- subset(results.mg,E %in% environmental_factors)
#
# subset(results.mg,results.mg[,6] < 0.05 / nrow(results.mg)) # how many significant in 20% cohort ?

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

# subset(results.mg,results.mg[,4] < 0.05 / nrow(results.mg))

fwrite(df.save,paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.validation.reversed.txt'),quote = F,na='NA',sep = '\t',row.names = F,col.names = T)


