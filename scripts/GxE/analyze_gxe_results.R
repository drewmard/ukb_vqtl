library(data.table)
phenoName <- 'monocyte.count.rint.ALL'
phenoName <- 'lymphocyte.count.rint.ALL'
phenoName <- 'neutrophil.count.rint.ALL'
envirFactor <- c('PA' ,'SB' ,'Smoking.E','sleep.duration', 'time.spent.outdoors' ,'tobacco.smoke.exposure', 'alcohol.freq.E', 'age' ,'sex')
for (i in 1:length(envirFactor)) {
  f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_res/ukbb.',phenoName,'.ALL.sub.',envirFactor[i],'.GxE.80.assoc.linear')
  df <- fread(f,data.table=F,stringsAsFactors=F)
  df.sub <- df[seq(3,nrow(df),by=3),]
  df.sub$envirFactor <- envirFactor[i]
  if (i==1) {
    df.save <- df.sub
  } else {
    df.save <- rbind(df.save,df.sub)
  }
}

df.save[order(df.save$P)[1:5],]


phenoName <- 'lymphocyte.count.rint.ALL'
f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/results/ukbb.gxe.',phenoName,'.',s,'.txt')
df <- fread(f,data.table = F,stringsAsFactors = F)
df[order(df$P_GxE)[1:5],]


