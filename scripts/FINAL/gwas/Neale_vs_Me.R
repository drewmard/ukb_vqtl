df.Neale <- read.table('/athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/21001_raw.gwas.imputed_v3.both_sexes.tsv.bgz',stringsAsFactors = F,header=T)
x=df.Neale[1:1045000,]
x[order(x$pval,decreasing=F),][1:10,]#RELEVANT_COLUMNS]

library(data.table)
df.Neale <- fread('~/Documents/Research/vQTL/UKB/21001_raw.gwas.imputed_v3.both_sexes.chr1.txt',data.table = F,stringsAsFactors = F)
df.mean <- fread('~/Documents/Research/vQTL/ukb_vqtl/output/imputed/results/ukbb.bmi.ALL.results.txt',data.table = F,stringsAsFactors = F)
df.mean.tmp <- df.mean

df.mean$Neale_SNP <- apply(df.mean[,c('CHR','BP','A2','A1')],1,function(x) paste(as.character(x),collapse = ':'))
df.mean.tmp$Neale_SNP <- apply(df.mean.tmp[,c('CHR','BP','A1','A2')],1,function(x) paste(as.character(x),collapse = ':'))

df.mean <- rbind(df.mean,df.mean.tmp)
df.mean$Neale_SNP <- str_replace_all(df.mean$Neale_SNP,' ','')
df.mg=merge(df.mean,df.Neale,by.x='Neale_SNP',by.y='variant')

i=which(df.mg$A1!=df.mg$minor_allele)
df.mg$beta[i] <- -df.mg$beta[i]

library(ggplot2)
ggplot(df.mg[sample(1:nrow(df.mg),10000,replace=F),],aes(x=BETA,y=beta)) + geom_point() + geom_smooth(col='red',method='lm')
g <- ggplot(df.mg[sample(1:nrow(df.mg),10000,replace=F),],aes(x=-log10(P),y=-log10(pval))) + geom_point() + geom_smooth(col='red',method='lm',lty='dashed',se=F) + theme_bw() + theme(panel.grid=element_blank()) + labs(x='Current study',y='Neale')
g <- ggplot(df.mg,aes(x=-log10(P),y=-log10(pval))) + geom_point() + geom_smooth(col='red',method='lm',lty='dashed',se=F) + theme_bw() + theme(panel.grid=element_blank()) + labs(x='Current study',y='Neale')
g
f <- '~/Documents/Research/vQTL/Neale_vs_Me.png'
png(f,width=6000,height=6000,res=1000)
g
dev.off()

df.mg$P.INFL <- -log10(df.mg$P/df.mg$pval)
df.mg[order(df.mg$P.INFL,decreasing = T)[1:5],]


median(df.mg$pval)
median(df.mg$P)
cor(-log10(df.mg$P),-log10(df.mg$pval))


x <- subset(df.mean,!(SNP %in% df.mg$SNP) & CHR==1)
x[order(x$P)[1:5],]

subset(df.Neale,variant=='1:177864136:A:T')

nrow(df.mg)
nrow(df.mean)
nrow(df.Neale)