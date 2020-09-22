library(data.table)
library(ggplot2)
library(cowplot)

s='80'; results.80 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.',s,'.pleiotropy.txt'),data.table = F,stringsAsFactors = F)
s='20'; results.20 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.',s,'.pleiotropy.txt'),data.table = F,stringsAsFactors = F)
df.results.save <- merge(results.80,results.20,by=c('SNP','E','pheno'))

df.results.save <- fread('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.ALL.pleiotropy.txt',data.table = F,stringsAsFactors = F)
# df.results.save <- fread('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.ALL.pleiotropy.used_bmi.txt',data.table = F,stringsAsFactors = F)
colnames(df.results.save)[4:7] <- c('BETA.80.DISEASE','P.80.DISEASE','BETA.20.DISEASE','P.20.DISEASE')
f='~/Documents/Research/vQTL/ukb_vqtl/output/sig_results/bmi.sig.GxE.txt'
results.mg <- fread(f,data.table = F,stringsAsFactors = F)
results.mg.mg <- merge(results.mg,df.results.save,by=c('SNP','E'))

# 1
results.mg.mg.diab <- subset(results.mg.mg,pheno=='Diabetes')
cor.test(results.mg.mg.diab$BETA.80.DISEASE,results.mg.mg.diab$BETA.20.DISEASE)
cor.test(results.mg.mg.diab$Estimate.x,results.mg.mg.diab$BETA.80.DISEASE)
cor.test(results.mg.mg.diab$Estimate.y,results.mg.mg.diab$BETA.20.DISEASE)
cor.test(results.mg.mg.diab$Estimate.x,results.mg.mg.diab$BETA.20.DISEASE)

g4 <- ggplot(results.mg.mg.diab,aes(x=Estimate.x,y=BETA.20.DISEASE)) + geom_point() + 
  geom_smooth(method='lm',se=F,col='red') +
  theme_bw() + theme(panel.grid=element_blank()) + 
  labs(x=expression(BMI: beta['D']),y=expression(Diabetes: beta['R']))

x=2
# png('~/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/diabetes_bmi_gxe_effects.used_bmi.png',width=2800*x,height=3000*x,res=1700)
png('~/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/diabetes_bmi_gxe_effects.png',width=2800*x,height=3000*x,res=1700)
g4
dev.off()


g1 <- ggplot(results.mg.mg.diab,aes(x=BETA.80.DISEASE,y=BETA.20.DISEASE)) + geom_point() + 
  geom_smooth(method='lm',se=F,col='red') +
  theme_bw() + theme(panel.grid=element_blank()) + 
  labs(x=expression(Diabetes: beta['80']),y=expression(Diabetes: beta['20']))
x=2
png('~/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/diabetes_gxe_effects.png',width=2800*x,height=3000*x,res=1400)
g1
dev.off()



