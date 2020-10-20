library(data.table)
library(ggplot2)
library(cowplot)
library(reshape2)
library(RColorBrewer)

# parameters
MAF1 <- 0.4
# MAF1 <- 0.1
MAF2 <- 0.4
nsim <- 1000; 
nindiv <- 250000
simulation_type='gxg'
genetic_variance_explained.vec <- seq(0.01,0.1,by=0.01);

# simulation_type <- 'gxg'
# phenotype_noise <- 'NORMAL' 
# dir='/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/simulation/'
# f=paste0(dir,'MAF1_',MAF1,'.MAF2_',MAF2,'.NSIM_',nsim,'.NINDIV_',nindiv,'.TYPE_',simulation_type,'.NOISE_',phenotype_noise,'.txt')
# df1 <- fread(f,data.table = F,stringsAsFactors = F)
# phenotype_noise <- 'CHISQ4' 
# dir='/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/simulation/'
# f=paste0(dir,'MAF1_',MAF1,'.MAF2_',MAF2,'.NSIM_',nsim,'.NINDIV_',nindiv,'.TYPE_',simulation_type,'.NOISE_',phenotype_noise,'.txt')
# df2 <- fread(f,data.table = F,stringsAsFactors = F)
# df.mg <- rbind(df1,df2)
# simulation_type <- 'mean'
# phenotype_noise <- 'NORMAL'
# dir='/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/simulation/'
# f=paste0(dir,'MAF1_',MAF1,'.MAF2_',MAF2,'.NSIM_',nsim,'.NINDIV_',nindiv,'.TYPE_',simulation_type,'.NOISE_',phenotype_noise,'.txt')
# df1 <- fread(f,data.table = F,stringsAsFactors = F)
# df.mg <- rbind(df.mg,df1)
# phenotype_noise <- 'CHISQ4' 
# dir='/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/simulation/'
# f=paste0(dir,'MAF1_',MAF1,'.MAF2_',MAF2,'.NSIM_',nsim,'.NINDIV_',nindiv,'.TYPE_',simulation_type,'.NOISE_',phenotype_noise,'.txt')
# df2 <- fread(f,data.table = F,stringsAsFactors = F)
# df.mg <- rbind(df.mg,df2)
# rm(df1);rm(df2)

simulation_type <- 'ALL'
phenotype_noise <- 'ALL'
dir='/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/simulation/'
f=paste0(dir,'MAF1_',MAF1,'.MAF2_',MAF2,'.NSIM_',nsim,'.NINDIV_',nindiv,'.TYPE_',simulation_type,'.NOISE_',phenotype_noise,'.txt')
df.mg <- fread(f,data.table = F,stringsAsFactors = F)
df.mg <- subset(df.mg,transformation=='orig')

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

# df.aggre <- aggregate(.~MAF1+MAF2+N+noise+h+type,df.mg[,-c(7)],function(x) {mean(x < 0.05)})
ind <- which(colnames(df.mg) %in% c('transformation','beta','gJLS',paste0(c('DRM','LT','BF','BT','FK','DGLM','TSSR','SVLM','gJLS'),'.time')))
df.aggre <- aggregate(.~MAF1+MAF2+N+noise+h+type,df.mg[,-ind],function(x) {mean(x < 0.05)})
df.aggre.melt <- melt(df.aggre[,c(4:ncol(df.aggre))],id.vars=c('h','noise','type'))
df.aggre.melt$shape <- 0; df.aggre.melt$shape[df.aggre.melt$variable=='DRM'] <- 1; df.aggre.melt$shape <- as.factor(df.aggre.melt$shape)
# g2 <- ggplot(subset(df.aggre.melt,noise=='NORMAL' & type=='gxg'),aes(x=h,y=value,col=variable,shape=shape)) + geom_line() + geom_point() + theme_bw() + theme(panel.grid = element_blank()) + labs(x='Variance explained by interaction effect',y='Power') + scale_color_brewer(palette='Dark2') + theme(legend.title = element_blank()) + scale_shape_manual(labels=c('Other','DRM')) # scale_shape_discrete(labels=c('Other','DRM')) + guides(col=guide_legend(shape=c(2,1,1,1,1,1)))#guides(col = guide_legend(override.aes=list(shape=c(2,1,1,1,1,1))))
mycolors = c(brewer.pal(name="Dark2", n = 6), brewer.pal(name="Paired", n = 3))
g2 <- ggplot(subset(df.aggre.melt,noise=='NORMAL' & type=='gxg'),aes(x=h,y=value,col=variable,shape=variable)) + geom_line() + geom_point() + theme_bw() + theme(panel.grid = element_blank(),legend.title = element_blank()) + labs(x='% variance explained by interaction effect',y='Power') + scale_color_manual(values=mycolors) + scale_shape_manual(labels=c('DRM','LT','BF','BT','FK','DGLM','TSSR','SVLM','gS'),values=c(17,16,16,16,16,16,16,16,16)) + scale_x_continuous(breaks=c(0,0.005,0.01,0.015,0.02),labels=c("0%","0.5%","1%","1.5%","2%"))
g3 <- ggplot(subset(df.aggre.melt,noise=='CHISQ4' & type=='gxg'),aes(x=h,y=value,col=variable,shape=variable)) + geom_line() + geom_point() + theme_bw() + theme(panel.grid = element_blank(),legend.title = element_blank()) + labs(x='% variance explained by interaction effect',y='Power') + scale_color_manual(values=mycolors) + scale_shape_manual(labels=c('DRM','LT','BF','BT','FK','DGLM','TSSR','SVLM','gS'),values=c(17,16,16,16,16,16,16,16,16)) + scale_x_continuous(breaks=c(0,0.005,0.01,0.015,0.02),labels=c("0%","0.5%","1%","1.5%","2%"))
g4 <- ggplot(subset(df.aggre.melt,noise=='NORMAL' & type=='mean'),aes(x=h,y=value,col=variable,shape=variable)) + geom_abline(slope=0,intercept=0.05,col='red',linetype='dashed') + geom_line() + geom_point() + theme_bw() + theme(panel.grid = element_blank(),legend.title = element_blank()) + labs(x='% variance explained by mean SNP effect',y='False Positive Rate') + scale_color_manual(values=mycolors) + scale_shape_manual(labels=c('DRM','LT','BF','BT','FK','DGLM','TSSR','SVLM','gS'),values=c(17,16,16,16,16,16,16,16,16)) + scale_x_continuous(breaks=c(0,0.005,0.01,0.015,0.02),labels=c("0%","0.5%","1%","1.5%","2%"))
g5 <- ggplot(subset(df.aggre.melt,noise=='CHISQ4' & type=='mean'),aes(x=h,y=value,col=variable,shape=variable)) + geom_abline(slope=0,intercept=0.05,col='red',linetype='dashed') + geom_line() + geom_point() + theme_bw() + theme(panel.grid = element_blank(),legend.title = element_blank()) + labs(x='% variance explained by mean SNP effect',y='False Positive Rate') + scale_color_manual(values=mycolors) + scale_shape_manual(labels=c('DRM','LT','BF','BT','FK','DGLM','TSSR','SVLM','gS'),values=c(17,16,16,16,16,16,16,16,16)) + scale_x_continuous(breaks=c(0,0.005,0.01,0.015,0.02),labels=c("0%","0.5%","1%","1.5%","2%"))
# plot_grid(g2,g3,g4,g5,ncol=2)
aggregate(.~MAF1+MAF2+N+noise+type,df.mg[,-ind],function(x) {mean(x < 0.05)})
x <- aggregate(.~MAF1+MAF2+N+noise+type,subset(df.mg[,-ind],h<=.01),function(x) {mean(x < 0.05)})
x <- subset(x,noise=='CHISQ4' & type=='gxg')
x$DRM/x$BF
x$DRM/x$SVLM
ind2 <- which(colnames(df.mg) %in% c('transformation','noise','beta','gJLS',paste0(c('DRM','LT','BF','BT','FK','DGLM','TSSR','SVLM','gJLS'),'.time')))
x <- aggregate(.~MAF1+MAF2+N+type,df.mg[,-ind2],function(x) {mean(x < 0.05)})
x <- subset(x,type=='gxg')
x$DRM;x$gS

# time to complete
df.mg <- subset(df.mg,transformation=='orig')
ind2 <- which(colnames(df.mg) %in% paste0(c('DRM','LT','BF','BT','FK','DGLM','TSSR','SVLM','gJLS'),'.time'))
x <- aggregate(df.mg[,ind2],list(df.mg$N),median)
x$gJLS.time/x$DRM


results.sub <- subset(df.mg,h<=0.01 & noise=='CHISQ4' & type=='gxg' & transformation=='orig')
t.test(results.sub$DRM < 0.05,results.sub$BF < 0.05,paired=T)$p.value
t.test(results.sub$DRM < 0.05,results.sub$SVLM < 0.05,paired=T)$p.value
t.test(results.sub$DRM < 0.05,results.sub$gS < 0.05,paired=T)$p.value
results.sub <- subset(df.mg,h<=0.01 & type=='gxg' & transformation=='orig')
t.test(results.sub$DRM < 0.05,results.sub$gS < 0.05,paired=T)$p.value

# just to see other half of simulations
aggregate(.~MAF1+MAF2+N+noise+type,subset(df.mg[,-ind],h>.01),function(x) {mean(x < 0.05)})

# aggregate(.~MAF1+MAF2+N+noise+type,df.mg[,-c(7)],function(x) {mean(x < 0.05)})

# g2 <- ggplot(subset(df.aggre.melt,noise=='NORMAL' & type=='gxg'),aes(x=h,y=value,col=variable,shape=variable)) + geom_line() + geom_point() + theme_bw() + theme(panel.grid = element_blank(),legend.title = element_blank()) + labs(x='Variance explained by interaction effect',y='Power') + scale_color_brewer(palette='Dark2') + scale_shape_manual(labels=c('DRM','LT','BF','BT','FK','DGLM','TSSR','SVLM','gS'),values=c(17,16,16,16,16,16,16,16,16)) 
# g3 <- ggplot(subset(df.aggre.melt,noise=='CHISQ4' & type=='gxg'),aes(x=h,y=value,col=variable,shape=variable)) + geom_line() + geom_point() + theme_bw() + theme(panel.grid = element_blank()) + labs(x='Variance explained by interaction effect',y='Power') + scale_color_brewer(palette='Dark2') + theme(legend.title = element_blank()) + scale_shape_manual(labels=c('DRM','LT','BT','FK','DGLM','TSSR'),values=c(17,16,16,16,16,16)) 
# g4 <- ggplot(subset(df.aggre.melt,noise=='NORMAL' & type=='mean'),aes(x=h,y=value,col=variable,shape=variable)) + geom_abline(slope=0,intercept=0.05,col='red',linetype='dashed') + geom_line() + geom_point() + theme_bw() + theme(panel.grid = element_blank()) + labs(x='Variance explained by mean SNP effect',y='False Positive Rate') + scale_color_brewer(palette='Dark2') + theme(legend.title = element_blank()) + scale_shape_manual(labels=c('DRM','LT','BT','FK','DGLM','TSSR'),values=c(17,16,16,16,16,16)) 
# g5 <- ggplot(subset(df.aggre.melt,noise=='CHISQ4' & type=='mean'),aes(x=h,y=value,col=variable,shape=variable)) + geom_abline(slope=0,intercept=0.05,col='red',linetype='dashed') + geom_line() + geom_point() + theme_bw() + theme(panel.grid = element_blank()) + labs(x='Variance explained by mean SNP effect',y='False Positive Rate') + scale_color_brewer(palette='Dark2')+ theme(legend.title = element_blank()) + scale_shape_manual(labels=c('DRM','LT','BT','FK','DGLM','TSSR'),values=c(17,16,16,16,16,16)) 

f='~/Documents/Research/vQTL/ukb_vqtl/output/simulation/power_fpr_plots.png'
png(f,width=5000,height=5000,res=500)
plot_grid(g2,g3,g4,g5,ncol=2)
plot_grid(g4,g2,g3,g5,ncol=2)
dev.off()

mycolors = c(brewer.pal(name="Dark2", n = 6), brewer.pal(name="Paired", n = 3))[c(1,3,8)]
g2 <- ggplot(subset(df.aggre.melt,noise=='NORMAL' & type=='gxg' & variable%in%c('BF','SVLM','DRM')),aes(x=h,y=value,col=variable)) + geom_line() + geom_point(size=rel(0.8)) + theme_bw() + theme(panel.grid = element_blank()) + labs(x='Variance explained by interaction effect',y='Power') + scale_color_manual(values=mycolors) + theme(legend.title = element_blank())

g2 <- ggplot(subset(df.aggre.melt,noise=='NORMAL' & type=='gxg' & variable%in%c('BF','SVLM','DRM')),aes(x=h,y=value,col=variable)) + geom_line() + geom_point(size=rel(0.8)) + theme_bw() + theme(panel.grid = element_blank()) + labs(x='Percent variance explained by interaction effect',y='Power') + scale_color_brewer(palette='Dark2')+ theme(legend.title = element_blank())+ scale_x_continuous(breaks=c(0,0.005,0.01,0.015,0.02),labels=c("0%","0.5%","1%","1.5%","2%"))
g3 <- ggplot(subset(df.aggre.melt,noise=='CHISQ4' & type=='gxg' & variable%in%c('BF','SVLM','DRM')),aes(x=h,y=value,col=variable)) + geom_line() + geom_point(size=rel(0.8)) + theme_bw() + theme(panel.grid = element_blank()) + labs(x='Percent variance explained by interaction effect',y='Power') + scale_color_brewer(palette='Dark2')+ theme(legend.title = element_blank())+ scale_x_continuous(breaks=c(0,0.005,0.01,0.015,0.02),labels=c("0%","0.5%","1%","1.5%","2%"))
f='~/Documents/Research/vQTL/ukb_vqtl/output/simulation/power_fpr_plots2.png'
# png(f,width=3750,height=1000,res=500)
# plot_grid(g2,g3,nrow=1)
png(f,width=2500,height=2500,res=500)
plot_grid(g2,g3,ncol=1)
dev.off()


ind <- which(colnames(df.mg) %in% paste0(c('DRM','LT','BF','BT','FK','DGLM','TSSR','SVLM','gJLS'),'.time'))
df.aggre <- aggregate(df.mg[,ind],list(rep(1,nrow(df.mg))),mean)

df.mg.melt <- melt(df.mg[,ind])
df.mg.melt[,1] <- substring(df.mg.melt[,1],1,nchar(as.character(df.mg.melt[,1]))-5)
df.mg.melt[df.mg.melt[,1]=='gJLS',1] <- 'gS'
mycolors = c(brewer.pal(name="Dark2", n = 6), brewer.pal(name="Paired", n = 3))
g2 <- ggplot(df.mg.melt,aes(x=variable,y=value,fill=variable)) + geom_boxplot(outlier.shape = NA,col='black') + theme_bw() + theme(panel.grid=element_blank(),legend.position = 'none') + labs(y='Elapsed time in seconds',x='Method')# + geom_jitter(alpha=0.1,col='black')
f='~/Documents/Research/vQTL/ukb_vqtl/output/simulation/power_fpr_plots2_v2.png'
png(f,width=2500,height=2500,res=500)
plot_grid(g3,g2,ncol=1)
dev.off()

