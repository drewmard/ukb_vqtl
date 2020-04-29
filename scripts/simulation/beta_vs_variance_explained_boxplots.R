library(data.table)
library(ggplot2)
library(cowplot)
library(reshape2)

# parameters
MAF1 <- 0.4
MAF2 <- 0.4
nsim <- 1000; 
nindiv <- 10000
simulation_type='gxg'
genetic_variance_explained.vec <- seq(0.01,0.1,by=0.01);
phenotype_noise <- 'NORMAL' 
dir='/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/simulation/'
f=paste0(dir,'MAF1_',MAF1,'.MAF2_',MAF2,'.NSIM_',nsim,'.NINDIV_',nindiv,'.TYPE_',simulation_type,'.NOISE_',phenotype_noise,'.txt')
df <- fread(f,data.table = F,stringsAsFactors = F)
ggplot(df,aes(x=as.factor(h),y=beta,fill=h)) + geom_boxplot(alpha=0.5) + theme_bw() + theme(panel.grid = element_blank(),legend.position = 'none') + labs(x='Variance Explained by GxG',y=expression(beta['var'])) +
  # scale_color_continuous(low='orange',high='blue')
  scale_fill_continuous(low='yellow',high='red')

df.aggre <- aggregate(.~MAF1+MAF2+N+noise+h+type,df[,-c(7)],function(x) {mean(x < 0.05)})
df.aggre.melt <- melt(df.aggre[,c(5,7:ncol(df.aggre))],id.vars='h')
ggplot(df.aggre.melt,aes(x=h,y=value,col=variable)) + geom_line() + geom_point() + theme_bw() + theme(panel.grid = element_blank()) + labs(x='Variance Explained by GxG',y='Power')


phenotype_noise <- 'CHISQ4' 
dir='/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/simulation/'
f=paste0(dir,'MAF1_',MAF1,'.MAF2_',MAF2,'.NSIM_',nsim,'.NINDIV_',nindiv,'.TYPE_',simulation_type,'.NOISE_',phenotype_noise,'.txt')
df <- fread(f,data.table = F,stringsAsFactors = F)
ggplot(df,aes(x=as.factor(h),y=beta,fill=h)) + geom_boxplot(alpha=0.5) + theme_bw() + theme(panel.grid = element_blank(),legend.position = 'none') + labs(x='Variance Explained by GxG',y=expression(beta['var'])) +
  # scale_color_continuous(low='orange',high='blue')
  scale_fill_continuous(low='yellow',high='red')

df.aggre <- aggregate(.~MAF1+MAF2+N+noise+h+type,df[,-c(7)],function(x) {mean(x < 0.05)})
df.aggre.melt <- melt(df.aggre[,c(5,7:ncol(df.aggre))],id.vars='h')
ggplot(df.aggre.melt,aes(x=h,y=value,col=variable)) + geom_line() + geom_point() + theme_bw() + theme(panel.grid = element_blank()) + labs(x='Variance Explained by GxG',y='Power')



