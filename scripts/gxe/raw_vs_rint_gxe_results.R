library(data.table)
f.raw <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.80.ext.more_snp.txt'
f.loge <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.80.ext.more_snp.log.txt'
f.rint <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.80.ext.more_snp.RINT.txt'
f.diet_raw <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.80.diet_score.more_snp.txt'
f.diet_rint <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.80.diet_score.more_snp.RINT.txt'
raw <- fread(f.raw,data.table = F,stringsAsFactors = F)
rint <- fread(f.rint,data.table=F,stringsAsFactors = F)
loge <- fread(f.loge,data.table=F,stringsAsFactors = F)

# diet_raw <-  fread(f.diet_raw,data.table = F,stringsAsFactors = F)
# diet_rint <-  fread(f.diet_rint,data.table = F,stringsAsFactors = F)
# diet_log <-  fread(f.diet_rint,data.table = F,stringsAsFactors = F)
# raw <- rbind(raw,diet_raw)
# rint <- rbind(rint,diet_rint)
# loge <- rbind(loge,diet_log)

colnames(raw)[2:3] <- c('BETA.RAW','P.RAW')
colnames(loge)[2:3] <- c('BETA.LOG','P.LOG')
colnames(rint)[2:3] <- c('BETA.RINT','P.RINT')

results.mg <- merge(raw,rint,by=c('SNP','E'))
results.mg <- merge(results.mg,loge)
results.mg$FDR.RAW <- p.adjust(results.mg$P.RAW,method = 'fdr')
results.mg$FDR.RINT <- p.adjust(results.mg$P.RINT,method = 'fdr')
results.mg$FDR.LOG <- p.adjust(results.mg$P.LOG,method = 'fdr')

g <- regexpr("_[^_]*$", results.mg$SNP)-1
results.mg$SNP <- substring(results.mg$SNP,1,g)
cor(results.mg[,c('BETA.RAW','BETA.RINT','BETA.LOG')])
cor(-log10(results.mg[,c('P.RAW','P.RINT','P.LOG')]))
mean(results.mg$FDR.RAW<0.1)
mean(results.mg$FDR.RINT<0.1)
mean(results.mg$FDR.LOG<0.1)
aggregate(results.mg[,c('FDR.RAW','FDR.RINT','FDR.LOG')],by=list(results.mg$E),function(x) mean(x<0.1))

pheno <- 'bmi'
f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/sig_results/',pheno,'.sig.txt')
df.mg2 <- fread(f,data.table = F,stringsAsFactors = F)
results.mg <- merge(results.mg,df.mg2,by.x='SNP',by.y='rs')

results.mg[order(results.mg$FDR.RAW,decreasing = F),][1:5,]

res.sub <- subset(results.mg,Mean.QTL==1 & Raw.vQTL==0)
sum(res.sub$FDR.RAW < 0.1)/nrow(res.sub)
sum(res.sub$FDR.RINT < 0.1)/nrow(res.sub)
sum(res.sub$FDR.LOG < 0.1)/nrow(res.sub)
res.sub <- subset(results.mg,Raw.vQTL==1)
sum(res.sub$FDR.RAW < 0.1)/nrow(res.sub)
sum(res.sub$FDR.RINT < 0.1)/nrow(res.sub)
sum(res.sub$FDR.LOG < 0.1)/nrow(res.sub)
res.sub <- subset(results.mg,Rint.vQTL==1 & Raw.vQTL==0)
sum(res.sub$FDR.RAW < 0.1)/nrow(res.sub)
sum(res.sub$FDR.RINT < 0.1)/nrow(res.sub)
sum(res.sub$FDR.LOG < 0.1)/nrow(res.sub)
res.sub <- subset(results.mg,dQTL==1 &  Raw.vQTL==0)
sum(res.sub$FDR.RAW < 0.1)/nrow(res.sub)
sum(res.sub$FDR.RINT < 0.1)/nrow(res.sub)
sum(res.sub$FDR.LOG < 0.1)/nrow(res.sub)

f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/','bmi','.GxE.','80','.ext.more_snp.multiplicative.txt')
multi <- fread(f,data.table = F,stringsAsFactors = F)
colnames(multi)[2:3] <- c('BETA','P')
multi$FDR <- p.adjust(multi$P,method = 'fdr')
mean(multi$FDR<0.1)



