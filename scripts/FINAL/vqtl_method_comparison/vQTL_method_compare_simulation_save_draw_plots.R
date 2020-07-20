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
gA <- ggplot(subset(df.mg,noise=='NORMAL' & type=='gxg'),aes(x=as.factor(h),y=beta,fill=h)) + geom_boxplot(alpha=0.5) + theme_bw() + theme(panel.grid = element_blank(),legend.position = 'none') + labs(x='Variance explained by mean SNP effect',y=expression(beta['var'])) +
  scale_fill_continuous(low='yellow',high='red')
gB <- ggplot(subset(df.mg,noise=='NORMAL' & type=='mean'),aes(x=as.factor(h),y=beta,fill=h)) + geom_boxplot(alpha=0.5) + theme_bw() + theme(panel.grid = element_blank(),legend.position = 'none') + labs(x='Variance explained by interaction effect',y=expression(beta['var'])) +
  scale_fill_continuous(low='yellow',high='red')
f='~/Documents/Research/vQTL/ukb_vqtl/output/simulation/variance_explained_by_vqtl_beta.png'
png(f,width=4000,height=4000,res=500)
plot_grid(gA,gB,ncol=1)
dev.off()

cor(subset(df.mg,noise=='NORMAL' & type=='gxg')[,c('h','beta')])

df.aggre <- aggregate(.~MAF1+MAF2+N+noise+h+type,df.mg[,-c(7)],function(x) {mean(x < 0.05)})
df.aggre.melt <- melt(df.aggre[,c(4:ncol(df.aggre))],id.vars=c('h','noise','type'))
df.aggre.melt$shape <- 0; df.aggre.melt$shape[df.aggre.melt$variable=='DRM'] <- 1; df.aggre.melt$shape <- as.factor(df.aggre.melt$shape)
# g2 <- ggplot(subset(df.aggre.melt,noise=='NORMAL' & type=='gxg'),aes(x=h,y=value,col=variable,shape=shape)) + geom_line() + geom_point() + theme_bw() + theme(panel.grid = element_blank()) + labs(x='Variance explained by interaction effect',y='Power') + scale_color_brewer(palette='Dark2') + theme(legend.title = element_blank()) + scale_shape_manual(labels=c('Other','DRM')) # scale_shape_discrete(labels=c('Other','DRM')) + guides(col=guide_legend(shape=c(2,1,1,1,1,1)))#guides(col = guide_legend(override.aes=list(shape=c(2,1,1,1,1,1))))
g2 <- ggplot(subset(df.aggre.melt,noise=='NORMAL' & type=='gxg'),aes(x=h,y=value,col=variable,shape=variable)) + geom_line() + geom_point() + theme_bw() + theme(panel.grid = element_blank()) + labs(x='Variance explained by interaction effect',y='Power') + scale_color_brewer(palette='Dark2') + theme(legend.title = element_blank()) + scale_shape_manual(labels=c('DRM','LT','BT','FK','DGLM','TSSR'),values=c(17,16,16,16,16,16)) 
g3 <- ggplot(subset(df.aggre.melt,noise=='CHISQ4' & type=='gxg'),aes(x=h,y=value,col=variable,shape=variable)) + geom_line() + geom_point() + theme_bw() + theme(panel.grid = element_blank()) + labs(x='Variance explained by interaction effect',y='Power') + scale_color_brewer(palette='Dark2') + theme(legend.title = element_blank()) + scale_shape_manual(labels=c('DRM','LT','BT','FK','DGLM','TSSR'),values=c(17,16,16,16,16,16)) 
g4 <- ggplot(subset(df.aggre.melt,noise=='NORMAL' & type=='mean'),aes(x=h,y=value,col=variable,shape=variable)) + geom_abline(slope=0,intercept=0.05,col='red',linetype='dashed') + geom_line() + geom_point() + theme_bw() + theme(panel.grid = element_blank()) + labs(x='Variance explained by mean SNP effect',y='False Positive Rate') + scale_color_brewer(palette='Dark2') + theme(legend.title = element_blank()) + scale_shape_manual(labels=c('DRM','LT','BT','FK','DGLM','TSSR'),values=c(17,16,16,16,16,16)) 
g5 <- ggplot(subset(df.aggre.melt,noise=='CHISQ4' & type=='mean'),aes(x=h,y=value,col=variable,shape=variable)) + geom_abline(slope=0,intercept=0.05,col='red',linetype='dashed') + geom_line() + geom_point() + theme_bw() + theme(panel.grid = element_blank()) + labs(x='Variance explained by mean SNP effect',y='False Positive Rate') + scale_color_brewer(palette='Dark2')+ theme(legend.title = element_blank()) + scale_shape_manual(labels=c('DRM','LT','BT','FK','DGLM','TSSR'),values=c(17,16,16,16,16,16)) 
plot_grid(g2,g3,g4,g5,ncol=2)
aggregate(.~MAF1+MAF2+N+noise+type,df.mg[,-c(7)],function(x) {mean(x < 0.05)})

f='~/Documents/Research/vQTL/ukb_vqtl/output/simulation/power_fpr_plots.png'
png(f,width=5000,height=5000,res=500)
plot_grid(g2,g3,g4,g5,ncol=2)
dev.off()

g2 <- ggplot(subset(df.aggre.melt,noise=='NORMAL' & type=='gxg' & variable%in%c('LT','DRM')),aes(x=h,y=value,col=variable)) + geom_line() + geom_point(size=rel(0.8)) + theme_bw() + theme(panel.grid = element_blank()) + labs(x='Variance explained by interaction effect',y='Power') + scale_color_brewer(palette='Dark2')+ theme(legend.title = element_blank())
g3 <- ggplot(subset(df.aggre.melt,noise=='CHISQ4' & type=='gxg' & variable%in%c('LT','DRM')),aes(x=h,y=value,col=variable)) + geom_line() + geom_point(size=rel(0.8)) + theme_bw() + theme(panel.grid = element_blank()) + labs(x='Variance explained by interaction effect',y='Power') + scale_color_brewer(palette='Dark2')+ theme(legend.title = element_blank())
f='~/Documents/Research/vQTL/ukb_vqtl/output/simulation/power_fpr_plots2.png'
# png(f,width=3750,height=1000,res=500)
# plot_grid(g2,g3,nrow=1)
png(f,width=2500,height=2500,res=500)
plot_grid(g2,g3,ncol=1)
dev.off()



