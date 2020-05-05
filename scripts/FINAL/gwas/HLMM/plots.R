library(ggplot2)
library(data.table)
library(cowplot)
mean <- fread('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/imputed/results/ukbb.bmi.ALL.results.txt',data.table = F,stringsAsFactors = F)
mean.rint <- fread('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/imputed/results/ukbb.bmi.rint.ALL.results.txt',data.table = F,stringsAsFactors = F)
var.rint <- fread('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/vGWAS_subset/ukbb.bmi.rint.ALL.vGWAS.txt',data.table = F,stringsAsFactors = F)
var.raw <- fread('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/vGWAS_subset/ukbb.bmi.ALL.vGWAS.txt',data.table = F,stringsAsFactors = F)

colnames(mean)[1] <- 'rs'
mean <- mean[,c('rs','CHR','BP','A1','A2','MAF','BETA','P')]
colnames(mean)[(ncol(mean)-1):ncol(mean)] <- paste0(colnames(mean)[(ncol(mean)-1):ncol(mean)],'.MEAN')
var.raw <- var.raw[,c('rs','BETA','P')]; colnames(var.raw)[2:3] <- paste0(colnames(var.raw)[2:3],'.VAR.RAW')
var.rint <- var.rint[,c('rs','BETA','P')]; colnames(var.rint)[2:3] <- paste0(colnames(var.rint)[2:3],'.VAR.RINT')
colnames(mean.rint)[1] <- 'rs'
mean.rint <- mean.rint[,c('rs','BETA','P')]
colnames(mean.rint)[2:3] <- paste0(colnames(mean.rint)[2:3],'.MEAN.RINT')

df.mg <- merge(merge(mean,var.raw,by='rs'),var.rint,by='rs')
df.mg <- merge(df.mg,mean.rint,by='rs')
df.mg.sub <- subset(df.mg,MAF > 0.05)

HLMM <- fread('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GWAS/HLMM_results/ukbb.bmi.rint.ALL.HLMM.dispersion_nochr6.txt',data.table = F,stringsAsFactors = F)
HLMM.sub <- HLMM[,c('SNP','add','var','dispersion',
               'add_pval','var_pval','av_pval','dispersion_pval')]

df.mg.sub <- merge(df.mg.sub,HLMM.sub,by.x='rs',by.y='SNP')
df.mg.sub[,c('add','var','dispersion')] <- -1*df.mg.sub[,c('add','var','dispersion')]
cor.mat <- cor(df.mg.sub[,c('add','var','dispersion','BETA.MEAN','BETA.VAR.RAW','BETA.VAR.RINT')])
cor.mat.melt <- melt(cor.mat)
cor.mat.melt[order(cor.mat.melt$value),]

df.mg.sub[,c('P.MEAN','P.VAR.RAW','P.VAR.RINT')] <- -log10(df.mg.sub[,c('P.MEAN','P.VAR.RAW','P.VAR.RINT')])
cor.mat <- cor(df.mg.sub[,c('add_pval','var_pval','av_pval','dispersion_pval','P.MEAN','P.VAR.RAW','P.VAR.RINT')])
cor.mat.melt <- melt(cor.mat)
cor.mat.melt[order(cor.mat.melt$value),]

df.mg.sub$dispersion_col <- floor(df.mg.sub$dispersion_pval)
df.mg.sub$dispersion_col[df.mg.sub$dispersion_col > -log10(1e-5)] <- 15
g1 <- ggplot(df.mg.sub,aes(x=P.MEAN,y=dispersion_pval,col=dispersion_col)) + geom_point(size=rel(0.75)) + geom_smooth(col='red',method='lm',se=F,size=rel(0.6))+ labs(x='-log10 raw muQTL p-values',y='-log10 DET p-values') + theme_bw() + theme(panel.grid=element_blank(),legend.position = 'none') + scale_colour_gradient2(mid='grey83')
g2 <- ggplot(df.mg.sub,aes(x=P.VAR.RINT,y=dispersion_pval,col=dispersion_col)) + geom_point(size=rel(0.75)) + geom_smooth(col='red',method='lm',se=F,size=rel(0.6))+ labs(x='-log10 RINT vQTL p-values',y='-log10 DET p-values') + theme_bw() + theme(panel.grid=element_blank(),legend.position = 'none') + scale_colour_gradient2(mid='grey83')
g3 <- ggplot(df.mg.sub,aes(x=BETA.MEAN,y=dispersion,col=dispersion_col)) + geom_point(size=rel(0.75)) + geom_smooth(col='red',method='lm',se=F,size=rel(0.6))+ labs(x='Raw muQTL effects',y='Dispersion effects') + theme_bw() + theme(panel.grid=element_blank(),legend.position = 'none') + scale_colour_gradient2(mid='grey83')
g4 <- ggplot(df.mg.sub,aes(x=BETA.VAR.RINT,y=dispersion,col=dispersion_col)) + geom_point(size=rel(0.75)) + geom_smooth(col='red',method='lm',se=F,size=rel(0.6))+ labs(x='RINT vQTL effects',y='Dispersion effects') + theme_bw() + theme(panel.grid=element_blank(),legend.position = 'none') + scale_colour_gradient2(mid='grey83')

x=1.5
png('~/Documents/Research/vQTL/ukb_vqtl/output/GWAS/HLMM_results/HLMM_correlation_plots.png',width=4000*x,height=4000*x,res=550*x)
plot_grid(g3,g1,g4,g2,ncol=2)
dev.off()

