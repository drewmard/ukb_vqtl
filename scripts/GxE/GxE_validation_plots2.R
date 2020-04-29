library(data.table)
library(ggplot2)
library(cowplot)
thres=0.05; suff <- ifelse(is.null(thres),'','PVAL.')
thres=NULL; suff <- ifelse(is.null(thres),'','PVAL.')

#1
thres=0.05; suff <- ifelse(is.null(thres),'','PVAL.')
f <- paste0('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.validation.mean.raw_matched_snp.full.',suff,'txt')
df.matched <- fread(f,data.table = F,stringsAsFactors = F)
g1<-ggplot(subset(df.matched,thres %in% c(0.05,0.01,0.005,0.001)),aes(x=-log10(thres),y=p)) + geom_line() + geom_point() + theme_bw() +
  labs(title='') +
  geom_ribbon(aes(ymin=lower,ymax=upper),linetype=2,alpha=0.1) + ylim(0,1) +
  theme(axis.title=element_blank(),panel.grid = element_blank()) + scale_x_continuous(breaks = c(1,2,3,4,5))
thres=NULL; suff <- ifelse(is.null(thres),'','PVAL.')

#2
f <- paste0('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.validation.mean.raw_matched_snp.full.',suff,'txt')
df.matched2 <- fread(f,data.table = F,stringsAsFactors = F)
g2<-ggplot(subset(df.matched2,thres %in% c(1,0.1,0.05,0.01,0.005,0.001)),aes(x=-log10(thres),y=p)) + geom_line() + geom_point() + theme_bw() +
  labs(title='') +
  geom_ribbon(aes(ymin=lower,ymax=upper),linetype=2,alpha=0.1) + ylim(0,1) + 
  theme(axis.title=element_blank(),panel.grid = element_blank())

#3
thres=0.05; suff <- ifelse(is.null(thres),'','PVAL.')
f <- paste0('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.validation.summary.full.',suff,'txt')
df.save <- fread(f,data.table = F,stringsAsFactors = F)
df.save$FDR_based <- 1; df.save$FDR_based[df.save$thres %in% c(1,0.1,0.05,0.01,0.005,0.001)] <- 0
g3<-ggplot(subset(df.save,thres <= 0.05),aes(x=-log10(thres),y=p)) + geom_line() + geom_point(aes(col=as.factor(FDR_based))) + theme_bw() +
  labs(title='') + scale_color_manual(values=c('black','red')) + 
  geom_ribbon(aes(ymin=lower,ymax=upper),linetype=2,alpha=0.1) + ylim(0,1) + geom_hline(yintercept=ifelse(is.null(thres),0.5,df.matched$p[df.matched$thres==0.05]),linetype='dashed',col='red') +
  theme(axis.title=element_blank(),panel.grid = element_blank(),legend.position = 'none') + scale_x_continuous(breaks = c(2,3,4,5))

#4
thres=NULL; suff <- ifelse(is.null(thres),'','PVAL.')
f <- paste0('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.validation.summary.full.',suff,'txt')
df.save2 <- fread(f,data.table = F,stringsAsFactors = F)
df.save2$FDR_based <- 1; df.save2$FDR_based[df.save$thres %in% c(1,0.1,0.05,0.01,0.005,0.001)] <- 0
g4<-ggplot(df.save2,aes(x=-log10(thres),y=p)) + geom_line() + geom_point(aes(col=as.factor(FDR_based))) + theme_bw() +
  labs(title='') + scale_color_manual(values=c('black','red')) + 
  geom_ribbon(aes(ymin=lower,ymax=upper),linetype=2,alpha=0.1) + ylim(0,1) + geom_hline(yintercept=ifelse(is.null(thres),0.5,0.025),linetype='dashed',col='red') +
  theme(axis.title=element_blank(),panel.grid = element_blank(),legend.position = 'none')
plot_grid(g2,g4,g1,g3,nrow = 1)

x=1
f <- paste0('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/GxE_validation.QTL_vs_null.png')
png(f,width = 3000,height=1000,res=300)
plot_grid(g2,g4,g1,g3,nrow = 1)
dev.off()

################

thres=0.05; suff <- ifelse(is.null(thres),'','PVAL.')
# thres=NULL; suff <- ifelse(is.null(thres),'','PVAL.')

f <- paste0('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.validation.summary.only_mean.',suff,'txt')
df.only_mean <- fread(f,data.table = F,stringsAsFactors = F); 
df.only_mean$FDR_based <- 1; df.only_mean$FDR_based[df.only_mean$thres %in% c(1,0.1,0.05,0.01,0.005,0.001)] <- 0
g5<-ggplot(subset(df.only_mean,thres <=0.05),aes(x=-log10(thres),y=p)) + geom_line() + geom_point(aes(col=as.factor(FDR_based))) + theme_bw() +
  labs(title='') + scale_color_manual(values=c('black','red')) + 
  geom_ribbon(aes(ymin=lower,ymax=upper),linetype=2,alpha=0.1) + ylim(0,1) + geom_hline(yintercept=ifelse(is.null(thres),0.5,0.025),linetype='dashed',col='red') +
  theme(axis.title=element_blank(),panel.grid = element_blank(),legend.position = 'none')

f <- paste0('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.validation.summary.mean.var_raw_matched_mean.',suff,'txt')
df.mean.var_raw_matched_mean <- fread(f,data.table = F,stringsAsFactors = F); 
df.mean.var_raw_matched_mean$FDR_based <- 1; df.mean.var_raw_matched_mean$FDR_based[df.mean.var_raw_matched_mean$thres %in% c(1,0.1,0.05,0.01,0.005,0.001)] <- 0
g6<-ggplot(subset(df.mean.var_raw_matched_mean,thres <= 0.05),aes(x=-log10(thres),y=p)) + geom_line() + geom_point(aes(col=as.factor(FDR_based))) + theme_bw() +
  labs(title='') + scale_color_manual(values=c('black','red')) + 
  geom_ribbon(aes(ymin=lower,ymax=upper),linetype=2,alpha=0.1) + ylim(0,1) + geom_hline(yintercept=ifelse(is.null(thres),0.5,0.025),linetype='dashed',col='red') +
  theme(axis.title=element_blank(),panel.grid = element_blank(),legend.position = 'none')

f <- paste0('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.validation.summary.var_raw.',suff,'txt')
df.var_raw <- fread(f,data.table = F,stringsAsFactors = F); 
df.var_raw$FDR_based <- 1; df.var_raw$FDR_based[df.var_raw$thres %in% c(1,0.1,0.05,0.01,0.005,0.001)] <- 0
g7<-ggplot(subset(df.var_raw,thres <= 0.05),aes(x=-log10(thres),y=p)) + geom_line() + geom_point(aes(col=as.factor(FDR_based))) + theme_bw() +
  labs(title='') + scale_color_manual(values=c('black','red')) + 
  geom_ribbon(aes(ymin=lower,ymax=upper),linetype=2,alpha=0.1) + ylim(0,1) + geom_hline(yintercept=ifelse(is.null(thres),0.5,0.025),linetype='dashed',col='red') +
  theme(axis.title=element_blank(),panel.grid = element_blank(),legend.position = 'none')

x=1
f <- paste0('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/GxE_validation.vqtl_vs_muqtl.png')
png(f,width = 3000,height=1000,res=300)
plot_grid(g5,g6,g7,nrow=1)
dev.off()

x=1
f <- paste0('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/GxE_validation.final_panel.png')
png(f,width = 3000,height=1000,res=300)
plot_grid(g1,g3,g5,g7,nrow=1)
dev.off()



# 
THRESHOLD=0.05
binom.test(
  prod(df.save[df.save$thres==THRESHOLD,c('p','N')]),
  df.save[df.save$thres==THRESHOLD,c('N')],
  df.matched$p[df.matched$thres==THRESHOLD]
)
df.save$p[df.save$thres==THRESHOLD]/df.matched$p[df.matched$thres==THRESHOLD]

THRESHOLD=0.05
binom.test(
  prod(df.save[6,c('p','N')]),
  df.save[6,c('N')],
  df.matched$p[7]
)
df.save$p[6]/df.matched$p[7]



THRESHOLD=0.05
binom.test(
  prod(df.var_raw[df.var_raw$thres==THRESHOLD,c('p','N')]),
  df.var_raw[df.var_raw$thres==THRESHOLD,c('N')],
  df.only_mean$p[df.only_mean$thres==THRESHOLD]
)
df.var_raw$p[df.var_raw$thres==THRESHOLD]/df.only_mean$p[df.only_mean$thres==THRESHOLD]


THRESHOLD=0.001
binom.test(
  prod(df.var_raw[df.var_raw$thres==THRESHOLD,c('p','N')]),
  df.var_raw[df.var_raw$thres==THRESHOLD,c('N')],
  df.mean.var_raw_matched_mean$p[df.mean.var_raw_matched_mean$thres==THRESHOLD]
)
df.var_raw$p[df.var_raw$thres==THRESHOLD]/df.mean.var_raw_matched_mean$p[df.mean.var_raw_matched_mean$thres==THRESHOLD]

