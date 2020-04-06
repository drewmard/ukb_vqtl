library(data.table)
library(ggplot2)
library(cowplot)
validation <- fread('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.validation.txt',data.table = F,stringsAsFactors = F)

g1 <- ggplot(validation,aes(x=-log10(thres),y=p)) + geom_point(col='orange') + theme_bw() +
  labs(x='-log10 p-value threshold cutoff',y='proportion w/ same sign in discovery & validation sets') +
  labs(title = 'All interactions w/ p-value < threshold')

g2 <- ggplot(validation,aes(x=-log10(thres),y=p2)) + geom_point(col='orange') + theme_bw() +
  labs(x='-log10 p-value cutoff',y='proportion w/ same sign in discovery & validation sets') +
  labs(title= 'All interactions w/ p-value > threshold')

plot_grid(g1,g2,ncol=2)

####################

f <- '/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.validation.summary.txt'
validation <- fread(f,data.table = F,stringsAsFactors = F)
library(reshape2)
validation.melt <- melt(validation,id.vars ='thres')

validation.melt.sub <- subset(validation.melt,thres > 1e-5)
validation.melt.sub <- subset(validation.melt,thres %in% unique(validation.melt$thres)[seq(1,60,by=5)])
g1 <- ggplot(validation.melt.sub,aes(x=-log10(thres),y=value)) + geom_line(aes(col=variable)) + theme_bw() +
  labs(x='-log10 p-value threshold cutoff',y='proportion w/ same sign in discovery & validation sets') +
  labs(title = 'All interactions w/ p-value < threshold') +
  scale_color_brewer(palette="Dark2")
g1

#################

f <- '/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxG_2/bmi.GxG.validation.summary.txt'
validation <- fread(f,data.table = F,stringsAsFactors = F)
validation <- validation[,c('thres','p','winner','p2','winner2')]
# colnames(validation) <- c('thres')
library(reshape2)
validation.melt <- melt(validation[,c('thres','p','winner','p2','winner2')],id.vars ='thres')

validation.melt.sub <- subset(validation.melt,thres %in% unique(validation.melt$thres)[seq(1,length(unique(validation.melt$thres)),by=5)])
g1 <- ggplot(subset(validation.melt.sub,variable %in% c('p','p2')),aes(x=-log10(thres),y=value)) + geom_line(aes(col=variable)) + theme_bw() +
  labs(x='-log10 p-value threshold cutoff',y='proportion w/ same sign in discovery & validation sets') +
  labs(title = 'All interactions w/ p-value < threshold') +
  scale_color_brewer(palette="Dark2") +
  geom_hline(yintercept =0.5,col='black',lty='dashed')
g1
g2 <- ggplot(subset(validation.melt.sub,variable %in% c('winner','winner2')),aes(x=-log10(thres),y=value)) + geom_line(aes(col=variable)) + theme_bw() +
  labs(x='-log10 p-value threshold cutoff',y='proportion w/ same sign in discovery & validation sets') +
  labs(title = 'All interactions w/ p-value < threshold') +
  scale_color_brewer(palette="Dark2") +
  geom_hline(yintercept =0.5,col='black',lty='dashed')
g2


##########################

# validation <- validation[,c('thres','p','lower','upper')]
# validation.melt <- melt(validation,id.vars ='thres')
# labs(x='-log10 p-value threshold cutoff',y='proportion w/ same sign in discovery & validation sets') +
#   labs(title = 'All interactions w/ p-value < threshold') +

library(data.table)
library(ggplot2)
library(reshape2)
library(cowplot)
  
f <- '/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.validation.summary.full.txt'
df.save <- fread(f,data.table = F,stringsAsFactors = F)
g1 <- ggplot(df.save,aes(x=-log10(thres),y=p)) + geom_line() + geom_point() + theme_bw() +
  labs(title='All') +
  geom_ribbon(aes(ymin=lower,ymax=upper),linetype=2,alpha=0.1) + ylim(0,1) + geom_hline(yintercept=0.5,linetype='dashed',col='red') +
  theme(axis.title=element_blank())

f <- '/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.validation.summary.var_raw.txt'
df.save.var_raw <- fread(f,data.table = F,stringsAsFactors = F)
g2 <- ggplot(df.save.var_raw,aes(x=-log10(thres),y=p)) + geom_line() + geom_point() + theme_bw() +
  labs(title='var_raw') +
  geom_ribbon(aes(ymin=lower,ymax=upper),linetype=2,alpha=0.1) + ylim(0,1) + geom_hline(yintercept=0.5,linetype='dashed',col='red') +
  theme(axis.title=element_blank())

f <- '/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.validation.summary.var_rint.txt'
df.save.var_rint <- fread(f,data.table = F,stringsAsFactors = F)
g3 <- ggplot(df.save.var_rint,aes(x=-log10(thres),y=p)) + geom_line() + geom_point() + theme_bw() +
  labs(title='var_rint') +
  geom_ribbon(aes(ymin=lower,ymax=upper),linetype=2,alpha=0.1) + ylim(0,1) + geom_hline(yintercept=0.5,linetype='dashed',col='red')+
  theme(axis.title=element_blank()) +
  theme(axis.title=element_blank())

f <- '/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.validation.summary.mean.var_raw.txt'
df.save.mean.var_raw <- fread(f,data.table = F,stringsAsFactors = F)
g4 <- ggplot(df.save.mean.var_raw,aes(x=-log10(thres),y=p)) + geom_line() + geom_point() + theme_bw() +
  labs(title='mean.var_raw') +
  geom_ribbon(aes(ymin=lower,ymax=upper),linetype=2,alpha=0.1) + ylim(0,1) + geom_hline(yintercept=0.5,linetype='dashed',col='red')+
  theme(axis.title=element_blank())

f <- '/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.validation.summary.mean.var_rint.txt'
df.save.mean.var_rint <- fread(f,data.table = F,stringsAsFactors = F)
g5 <- ggplot(df.save.mean.var_rint,aes(x=-log10(thres),y=p)) + geom_line() + geom_point() + theme_bw() +
  labs(title='mean.var_rint') +
  geom_ribbon(aes(ymin=lower,ymax=upper),linetype=2,alpha=0.1) + ylim(0,1) + geom_hline(yintercept=0.5,linetype='dashed',col='red')+
  theme(axis.title=element_blank())

plot_grid(g1,g2,g3,g4,g5,nrow = 1)


f <- '/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.validation.summary.only_mean.txt'
df.save.only_mean <- fread(f,data.table = F,stringsAsFactors = F)
g2.1 <- ggplot(df.save.only_mean,aes(x=-log10(thres),y=p)) + geom_line() + geom_point() + theme_bw() +
  labs(title='mean only, no var') +
  geom_ribbon(aes(ymin=lower,ymax=upper),linetype=2,alpha=0.1) + ylim(0,1) + geom_hline(yintercept=0.5,linetype='dashed',col='red') +
  theme(axis.title=element_blank())

f <- '/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.validation.summary.mean.plus_var.txt'
df.save.mean.plus_var <- fread(f,data.table = F,stringsAsFactors = F)
g2.2 <- ggplot(df.save.mean.plus_var,aes(x=-log10(thres),y=p)) + geom_line() + geom_point() + theme_bw() +
  labs(title='mean+var rint/raw') +
  geom_ribbon(aes(ymin=lower,ymax=upper),linetype=2,alpha=0.1) + ylim(0,1) + geom_hline(yintercept=0.5,linetype='dashed',col='red') +
  theme(axis.title=element_blank())

f <- '/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.validation.summary.some_var.txt'
df.save.some_var <- fread(f,data.table = F,stringsAsFactors = F)
g2.3 <- ggplot(df.save.some_var,aes(x=-log10(thres),y=p)) + geom_line() + geom_point() + theme_bw() +
  labs(title='var rint/raw') +
  geom_ribbon(aes(ymin=lower,ymax=upper),linetype=2,alpha=0.1) + ylim(0,1) + geom_hline(yintercept=0.5,linetype='dashed',col='red')+
  theme(axis.title=element_blank()) +
  theme(axis.title=element_blank())

plot_grid(g2.1,g2.2,g2.3,nrow = 1)

df.save.melt <- melt(df.save[,c('thres','winner','winner2')],id.vars = 'thres')
ggplot(df.save.melt,aes(x=-log10(thres),y=value,col=variable)) + geom_line() + geom_point() + theme_bw() +
  labs(title='All') +
  ylim(0,1) + geom_hline(yintercept=0.5,linetype='dashed',col='red') +
  theme(axis.title=element_blank()) + theme(legend.position = 'none')



plot_grid(g1,g2,g3,g4,g5,g2.1,g2.2,g2.3,nrow = 1)
plot_grid(g1,g2.1,g2.2,g4,g5,g2,g3,g2.3,nrow = 1)


f <- '/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.validation.summary.strict_mean.no_var.txt'
df.save.strict_mean.no_var <- fread(f,data.table = F,stringsAsFactors = F); df.save.strict_mean.no_var$snps <- 'no_var'
f <- '/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.validation.summary.strict_mean.plus_var.txt'
df.save.strict_mean.plus_var <- fread(f,data.table = F,stringsAsFactors = F); df.save.strict_mean.plus_var$snps <- 'plus_var'
df.save.strict_mean.var_diff <- rbind(df.save.strict_mean.no_var[,c('thres','p','lower','upper','snps')],df.save.strict_mean.plus_var[,c('thres','p','lower','upper','snps')])

ggplot(df.save.strict_mean.var_diff,aes(x=-log10(thres),y=p,col=snps)) + geom_line() + geom_point() + theme_bw() +
  labs(title='highly significant mean qtls') +
  ylim(0,1) + geom_hline(yintercept=0.5,linetype='dashed',col='red') +
  geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.1) + ylim(0,1) +
  theme(axis.title=element_blank()) + theme(legend.position = 'none')

df.save.only_mean$snps='no_var'
df.save.mean.plus_var$snps='plus_var'
df.save.mean.var_diff <- rbind(df.save.only_mean[,c('thres','p','lower','upper','snps')],df.save.mean.plus_var[,c('thres','p','lower','upper','snps')])
ggplot(df.save.mean.var_diff,aes(x=-log10(thres),y=p,col=snps)) + geom_line() + geom_point() + theme_bw() +
  labs(title='mean qtls') +
  ylim(0,1) + geom_hline(yintercept=0.5,linetype='dashed',col='red') +
  geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.1) + ylim(0,1) +
  theme(axis.title=element_blank()) + theme(legend.position = 'none')

####################################
########################################

f <- '/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxG_2/bmi.GxG.validation.summary.full.txt'
gxg_validation_statistics.df <- fread(f,data.table = F,stringsAsFactors = F)
f <- '/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxG_2/bmi.GxG.validation.summary.mean.txt'
gxg_validation_statistics.df.mean <- fread(f,data.table = F,stringsAsFactors = F)
f <- '/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxG_2/bmi.GxG.validation.summary.var_raw.txt'
gxg_validation_statistics.df.var_raw <- fread(f,data.table = F,stringsAsFactors = F)
f <- '/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxG_2/bmi.GxG.validation.summary.var_rint.txt'
gxg_validation_statistics.df.var_rint <- fread(f,data.table = F,stringsAsFactors = F)

g4.1 <- ggplot(gxg_validation_statistics.df,aes(x=-log10(thres),y=p)) + geom_line() + geom_point() + theme_bw() +
  labs(title='all') +
  geom_ribbon(aes(ymin=lower,ymax=upper),linetype=2,alpha=0.1) + ylim(0,1) + geom_hline(yintercept=0.5,linetype='dashed',col='red') +
  theme(axis.title=element_blank())
g4.2 <- ggplot(gxg_validation_statistics.df.mean,aes(x=-log10(thres),y=p)) + geom_line() + geom_point() + theme_bw() +
  labs(title='mean') +
  geom_ribbon(aes(ymin=lower,ymax=upper),linetype=2,alpha=0.1) + ylim(0,1) + geom_hline(yintercept=0.5,linetype='dashed',col='red') +
  theme(axis.title=element_blank())
g4.3 <- ggplot(gxg_validation_statistics.df.var_raw,aes(x=-log10(thres),y=p)) + geom_line() + geom_point() + theme_bw() +
  labs(title='var raw') +
  geom_ribbon(aes(ymin=lower,ymax=upper),linetype=2,alpha=0.1) + ylim(0,1) + geom_hline(yintercept=0.5,linetype='dashed',col='red') +
  theme(axis.title=element_blank())
g4.4 <- ggplot(gxg_validation_statistics.df.var_rint,aes(x=-log10(thres),y=p)) + geom_line() + geom_point() + theme_bw() +
  labs(title='var rint') +
  geom_ribbon(aes(ymin=lower,ymax=upper),linetype=2,alpha=0.1) + ylim(0,1) + geom_hline(yintercept=0.5,linetype='dashed',col='red') +
  theme(axis.title=element_blank())


plot_grid(g4.1,g4.2,g4.3,g4.4,nrow = 1)

#################################################################

f <- '/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.ALL_RESULTS.trim.txt'
df <- fread(f,data.table = F,stringsAsFactors = F)
library(qqman)
qq(df[,4])
qq(df[,6])



