library(data.table)
phenoName='lymphocyte.count'
f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/results/ukbb.',phenoName,'.results.txt')
df <- fread(f,stringsAsFactors = F,data.table = F)

# phenoName='monocyte.count.rint'
f.rint <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/results/ukbb.',paste0(phenoName,'.rint'),'.results.txt')
df.rint <- fread(f.rint,stringsAsFactors = F,data.table = F)

df <- subset(df,MAF > 0.05 & (NMISS > max(df$NMISS)*0.9))
df.rint <- subset(df.rint,MAF > 0.05 & (NMISS > max(df.rint$NMISS)*0.9))

df2 <- df[,c('SNP','CHR','BP','MAF','BETA.x','P.x','BETA.y','P.y')]
colnames(df2)[(ncol(df2)-3):ncol(df2)] <- c('BETA.MEAN','P.MEAN','BETA.VAR','P.VAR')
df2.rint <- df.rint[,c('SNP','CHR','BP','BETA.x','P.x','BETA.y','P.y')]
colnames(df2.rint)[(ncol(df2.rint)-3):ncol(df2.rint)] <- c('BETA.MEAN.RINT','P.MEAN.RINT','BETA.VAR.RINT','P.VAR.RINT')

df.mg <- merge(df2,df2.rint,by=c('SNP','CHR','BP'))

subset(df.mg,P.VAR.RINT < 1e-3 & SNP %in% results$vQTL)$SNP
subset(results,vQTL %in% subset(df.mg,P.VAR.RINT < 1e-3 & SNP %in% results$vQTL)$SNP)
x <- subset(results,vQTL %in% subset(df.mg,P.VAR.RINT < 1e-3 & SNP %in% results$vQTL)$SNP)
print(x[order(x$P_GxE),][1:10,])

cor.test(df.mg$BETA.MEAN,df.mg$BETA.VAR,use = 'pairwise.complete.obs')
cor.test(df.mg$BETA.MEAN.RINT,df.mg$BETA.VAR.RINT,use = 'pairwise.complete.obs')

head(df.mg[order(df.mg$P.VAR),],30)

x <- subset(df.mg,P.VAR < 5e-8 & P.VAR.RINT < 1e-3)


subset(df.mg,SNP%in%c('rs34890930','rs8026803'))
head(df.mg[order(df.mg$P.VAR.RINT),],10)

phenoName <- 'lymphocyte.count.transform_all'
f.save <- paste0('/athena/elementolab/scratch/anm2868/vQTL/UKB/results/ukbb.',phenoName,'.results.txt')
fwrite(df.mg,f.save,col.names = T,row.names = F,na='NA',quote=F,sep='\t')


user_direc <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl'
f.geno <- paste0(user_direc,'/output/GWAS/subset/',phenoName,'/ukbb.ALL_vQTL.raw')
df.geno <- fread(f.geno,data.table = F,stringsAsFactors = F,check.names = T)



