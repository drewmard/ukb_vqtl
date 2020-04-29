library(data.table)
library(ggplot2)
library(cowplot)
library(reshape2)

# parameters
MAF1 <- 0.4
# MAF1 <- 0.1
MAF2 <- 0.4
nsim <- 1000; 
nindiv <- 10000
simulation_type='gxg'
genetic_variance_explained.vec <- seq(0.01,0.1,by=0.01);

simulation_type <- 'gxg'
phenotype_noise <- 'NORMAL' 
dir='/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/simulation/'
f=paste0(dir,'MAF1_',MAF1,'.MAF2_',MAF2,'.NSIM_',nsim,'.NINDIV_',nindiv,'.TYPE_',simulation_type,'.NOISE_',phenotype_noise,'.txt')
df1 <- fread(f,data.table = F,stringsAsFactors = F)
phenotype_noise <- 'CHISQ4' 
dir='/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/simulation/'
f=paste0(dir,'MAF1_',MAF1,'.MAF2_',MAF2,'.NSIM_',nsim,'.NINDIV_',nindiv,'.TYPE_',simulation_type,'.NOISE_',phenotype_noise,'.txt')
df2 <- fread(f,data.table = F,stringsAsFactors = F)
df.mg <- rbind(df1,df2)
simulation_type <- 'mean'
phenotype_noise <- 'NORMAL'
dir='/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/simulation/'
f=paste0(dir,'MAF1_',MAF1,'.MAF2_',MAF2,'.NSIM_',nsim,'.NINDIV_',nindiv,'.TYPE_',simulation_type,'.NOISE_',phenotype_noise,'.txt')
df1 <- fread(f,data.table = F,stringsAsFactors = F)
df.mg <- rbind(df.mg,df1)
phenotype_noise <- 'CHISQ4' 
dir='/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/simulation/'
f=paste0(dir,'MAF1_',MAF1,'.MAF2_',MAF2,'.NSIM_',nsim,'.NINDIV_',nindiv,'.TYPE_',simulation_type,'.NOISE_',phenotype_noise,'.txt')
df2 <- fread(f,data.table = F,stringsAsFactors = F)
df.mg <- rbind(df.mg,df2)
rm(df1);rm(df2)
# phenotype_noise <- 'CHISQ4' 
# dir='/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/simulation/'
# f=paste0(dir,'MAF1_',MAF1,'.MAF2_',MAF2,'.NSIM_',nsim,'.NINDIV_',nindiv,'.TYPE_',simulation_type,'.NOISE_',phenotype_noise,'.txt')
# df2 <- fread(f,data.table = F,stringsAsFactors = F)
# df.mg <- rbind(df1,df2)


# ggplot(df.mg,aes(x=as.factor(h),y=beta,fill=noise)) + geom_boxplot(alpha=0.5) + theme_bw() + theme(panel.grid = element_blank()) + labs(x='Variance Explained by GxG',y=expression(beta['var']))# +
gA <- ggplot(subset(df.mg,noise=='NORMAL' & type=='gxg'),aes(x=as.factor(h),y=beta,fill=h)) + geom_boxplot(alpha=0.5) + theme_bw() + theme(panel.grid = element_blank(),legend.position = 'none') + labs(x='Variance explained by simulated GxG',y=expression(beta['var'])) +
  scale_fill_continuous(low='yellow',high='red')
gB <- ggplot(subset(df.mg,noise=='NORMAL' & type=='mean'),aes(x=as.factor(h),y=beta,fill=h)) + geom_boxplot(alpha=0.5) + theme_bw() + theme(panel.grid = element_blank(),legend.position = 'none') + labs(x='Variance explained by simulated mean QTL',y=expression(beta['var'])) +
  scale_fill_continuous(low='yellow',high='red')
f='~/Documents/Research/vQTL/ukb_vqtl/output/simulation/variance_explained_by_vqtl_beta.png'
png(f,width=4000,height=4000,res=500)
plot_grid(gA,gB,ncol=1)
dev.off()

cor(subset(df.mg,noise=='NORMAL' & type=='gxg')[,c('h','beta')])

df.aggre <- aggregate(.~MAF1+MAF2+N+noise+h+type,df.mg[,-c(7)],function(x) {mean(x < 0.05)})
df.aggre.melt <- melt(df.aggre[,c(4:ncol(df.aggre))],id.vars=c('h','noise','type'))
g2 <- ggplot(subset(df.aggre.melt,noise=='NORMAL' & type=='gxg'),aes(x=h,y=value,col=variable)) + geom_line() + geom_point() + theme_bw() + theme(panel.grid = element_blank()) + labs(x='Variance Explained by GxG',y='Power') + scale_color_brewer(palette='Dark2') + theme(legend.title = element_blank())
g3 <- ggplot(subset(df.aggre.melt,noise=='CHISQ4' & type=='gxg'),aes(x=h,y=value,col=variable)) + geom_line() + geom_point() + theme_bw() + theme(panel.grid = element_blank()) + labs(x='Variance Explained by GxG',y='Power') + scale_color_brewer(palette='Dark2') + theme(legend.title = element_blank())
g4 <- ggplot(subset(df.aggre.melt,noise=='NORMAL' & type=='mean'),aes(x=h,y=value,col=variable)) + geom_line() + geom_point() + theme_bw() + theme(panel.grid = element_blank()) + labs(x='Variance Explained by mean QTL',y='False Positive Rate') + geom_abline(slope=0,intercept=0.05,col='red',linetype='dashed') + scale_color_brewer(palette='Dark2') + theme(legend.title = element_blank())
g5 <- ggplot(subset(df.aggre.melt,noise=='CHISQ4' & type=='mean'),aes(x=h,y=value,col=variable)) + geom_line() + geom_point() + theme_bw() + theme(panel.grid = element_blank()) + labs(x='Variance Explained by mean QTL',y='False Positive Rate') + geom_abline(slope=0,intercept=0.05,col='red',linetype='dashed') + scale_color_brewer(palette='Dark2')+ theme(legend.title = element_blank())
plot_grid(g2,g3,g4,g5,ncol=2)
aggregate(.~MAF1+MAF2+N+noise+type,df.mg[,-c(7)],function(x) {mean(x < 0.05)})

f='~/Documents/Research/vQTL/ukb_vqtl/output/simulation/power_fpr_plots.png'
png(f,width=5000,height=5000,res=500)
plot_grid(g2,g3,g4,g5,ncol=2)
dev.off()

g2 <- ggplot(subset(df.aggre.melt,noise=='NORMAL' & type=='gxg' & variable%in%c('LT','DRM')),aes(x=h,y=value,col=variable)) + geom_line() + geom_point() + theme_bw() + theme(panel.grid = element_blank()) + labs(x='Variance Explained by GxG',y='Power') + scale_color_brewer(palette='Dark2')+ theme(legend.title = element_blank())
g3 <- ggplot(subset(df.aggre.melt,noise=='CHISQ4' & type=='gxg' & variable%in%c('LT','DRM')),aes(x=h,y=value,col=variable)) + geom_line() + geom_point() + theme_bw() + theme(panel.grid = element_blank()) + labs(x='Variance Explained by GxG',y='Power') + scale_color_brewer(palette='Dark2')+ theme(legend.title = element_blank())
f='~/Documents/Research/vQTL/ukb_vqtl/output/simulation/power_fpr_plots2.png'
png(f,width=5000,height=1000,res=500)
plot_grid(g2,g3,ncol=2)
dev.off()



