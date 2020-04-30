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

df.mg.sub$var.raw_col <- floor(-log10(df.mg.sub$P.VAR.RAW))
df.mg.sub$var.raw_col[df.mg.sub$var.raw_col > -log10(5e-8)] <- 15
g1 <- ggplot(df.mg.sub,aes(x=BETA.MEAN,y=BETA.VAR.RAW,col=var.raw_col)) + geom_point(size=rel(0.75)) + geom_smooth(col='red',method='lm',se=F,size=rel(0.6)) + labs(x='muQTL effects',y='Raw vQTL effects') + theme_bw() + theme(panel.grid=element_blank(),legend.position = 'none') + scale_colour_gradient2(mid='grey83')
g2 <- ggplot(df.mg.sub,aes(x=-log10(P.MEAN),y=-log10(P.VAR.RAW),col=(var.raw_col))) + geom_point(size=rel(0.75)) + geom_smooth(col='red',method='lm',se=F,size=rel(0.6)) + labs(x='-log10 muQTL p-values',y='-log10 raw vQTL p-values') + theme_bw() + theme(panel.grid=element_blank(),legend.position = 'none') + scale_colour_gradient2(mid='grey83')
df.mg.sub$var.rint_col <- floor(-log10(df.mg.sub$P.VAR.RINT))
df.mg.sub$var.rint_col[df.mg.sub$var.rint_col > -log10(1e-5)] <- 15
g3 <- ggplot(df.mg.sub,aes(x=-log10(P.MEAN),y=-log10(P.VAR.RINT),col=(var.rint_col))) + geom_point(size=rel(0.75)) + geom_smooth(col='red',method='lm',se=F,size=rel(0.6)) + labs(x='-log10 muQTL p-values',y='-log10 RINT vQTL p-values') + theme_bw() + theme(panel.grid=element_blank(),legend.position = 'none') + scale_colour_gradient2(mid='grey83')

x=2
png('~/Documents/Research/vQTL/ukb_vqtl/output/GWAS/mean_var_PANEL.png',width=5000*x,height=1500*x,res=450*x)
plot_grid(g1,g2,g3,nrow=1)
dev.off()

#########################

g1 <- ggplot(df.mg.sub,aes(x=BETA.MEAN,y=BETA.VAR.RAW,col=var.rint_col)) + geom_point(size=rel(0.75)) + geom_smooth(col='red',method='lm',se=F,size=rel(0.6)) + labs(x='muQTL effects',y='Raw vQTL effects') + theme_bw() + theme(panel.grid=element_blank(),legend.position = 'none') + scale_colour_gradient2(mid='grey83')
g2 <- ggplot(df.mg.sub,aes(x=BETA.MEAN,y=BETA.VAR.RINT,col=var.rint_col)) + geom_point(size=rel(0.75)) + geom_smooth(col='red',method='lm',se=F,size=rel(0.6)) + labs(x='muQTL effects',y='RINT vQTL effects') + theme_bw() + theme(panel.grid=element_blank(),legend.position = 'none') + scale_colour_gradient2(mid='grey83')
g3 <- ggplot(df.mg.sub,aes(x=-log10(P.MEAN),y=-log10(P.VAR.RAW),col=(var.rint_col))) + geom_point(size=rel(0.75)) + geom_smooth(col='red',method='lm',se=F,size=rel(0.6)) + labs(x='-log10 muQTL p-values',y='-log10 raw vQTL p-values') + theme_bw() + theme(panel.grid=element_blank(),legend.position = 'none') + scale_colour_gradient2(mid='grey83')

x=2
png('~/Documents/Research/vQTL/ukb_vqtl/output/GWAS/mean_var_PANEL_2.png',width=5000*x,height=1500*x,res=450*x)
plot_grid(g1,g2,g3,nrow=1)
dev.off()

###################

df.mg.sub$mean.raw_col <- floor(-log10(df.mg.sub$P.MEAN))
df.mg.sub$mean.raw_col[df.mg.sub$P.MEAN > -log10(5e-8)] <- 15

g1 <- ggplot(df.mg.sub,aes(x=BETA.VAR.RAW,y=BETA.VAR.RINT,col=mean.raw_col)) + geom_point(size=rel(0.75)) + geom_smooth(col='red',method='lm',se=F,size=rel(0.6)) + labs(x='Raw vQTL effects',y='RINT vQTL effects') + theme_bw() + theme(panel.grid=element_blank(),legend.position = 'none') + scale_colour_gradient2(mid='grey83')
g2 <- ggplot(df.mg.sub,aes(x=-log10(P.VAR.RAW),y=-log10(P.VAR.RINT),col=mean.raw_col)) + geom_point(size=rel(0.75)) + geom_smooth(col='red',method='lm',se=F,size=rel(0.6)) + labs(x='-log10 raw vQTL p-values',y='-log10 RINT vQTL p-values') + theme_bw() + theme(panel.grid=element_blank(),legend.position = 'none') + scale_colour_gradient2(mid='grey83')

g3 <- ggplot(df.mg.sub,aes(x=BETA.MEAN,y=BETA.MEAN.RINT,col=var.raw_col)) + geom_point(size=rel(0.75)) + geom_smooth(col='red',method='lm',se=F,size=rel(0.6)) +labs(x='Raw muQTL effects',y='RINT muQTL effects') + theme_bw() + theme(panel.grid=element_blank(),legend.position = 'none') + scale_colour_gradient2(mid='grey83')
g4 <- ggplot(df.mg.sub,aes(x=-log10(P.MEAN),y=-log10(P.MEAN.RINT),col=var.raw_col)) + geom_point(size=rel(0.75)) + geom_smooth(col='red',method='lm',se=F,size=rel(0.6)) +  labs(x='-log10 raw muQTL p-values',y='-log10 RINT muQTL p-values') + theme_bw() + theme(panel.grid=element_blank(),legend.position = 'none') + scale_colour_gradient2(mid='grey83')


x=2
png('~/Documents/Research/vQTL/ukb_vqtl/output/GWAS/mean_var_PANEL_3.png',width=3000*x,height=3000*x,res=450*x)
plot_grid(g1,g2,g3,g4,nrow=2)
dev.off()

