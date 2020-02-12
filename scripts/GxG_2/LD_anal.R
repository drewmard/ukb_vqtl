library(data.table)
phenoName <- 'lymphocyte.count'

s <- '80'
f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.',phenoName,'.ALL.sub.GxG.',s,'.epi.qt')
df.80 <- fread(f,data.table=F,stringsAsFactors = F)
colnames(df.80)[5:7] <- paste0(colnames(df.80)[5:7],'.',s)
df.80.sub <- subset(df.80,CHR1==CHR2)
subset(df.80,CHR1!=CHR2)[order(subset(df.80,CHR1!=CHR2)$P.80)[1:5],]

s <- '20'
df.20 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.',phenoName,'.ALL.sub.GxG.',s,'.epi.qt'),data.table=F,stringsAsFactors = F)
colnames(df.20)[5:7] <- paste0(colnames(df.20)[5:7],'.',s)
df <- merge(df.80,df.20,by=c('CHR1','SNP1','CHR2','SNP2'))
df.sub <- subset(df,CHR1!=CHR2)
df.sub[order(df.sub$P.80)[1:5],] 


f <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.lymphocyte.count.ALL.sub.LD.ld'
LD <- fread(f,data.table = F,stringsAsFactors = F)
LD.old <- LD
x <- LD[,1:3]
LD[,1:3] <- LD[,4:6]
LD[,4:6] <- x
LD.full <- data.frame(rbind(LD.old,LD))


df.80.sub.mg <- merge(df.80.sub,LD.full,by.x=c('SNP1','SNP2'),by.y=c('SNP_A','SNP_B'),all.x=TRUE)
df.80.sub.mg[order(df.80.sub.mg$P.80)[1:5],]

df.80.sub.mg.sub <- subset(df.80.sub.mg,CHR1!=6)
df.80.sub.mg.sub[order(df.80.sub.mg.sub$P.80)[1:5],] 


df.80.sub.mg.sub <- subset(df.80.sub.mg,CHR1==19)
df.80.sub.mg.sub[order(df.80.sub.mg.sub$P.80)[1:30],] 


df.80.sub.mg.sub$SNP1 


