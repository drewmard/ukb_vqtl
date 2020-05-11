library(data.table)
library(ggplot2)
library(cowplot)
thres=0.05; suff <- ifelse(is.null(thres),'','PVAL.')
# thres=NULL; suff <- ifelse(is.null(thres),'','PVAL.')

#1
f <- paste0('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxG_2/ukbb.bmi.merged_subset2.GxG.FULL.txt')
df <- fread(f,data.table = F,stringsAsFactors = F)
df <- subset(df,!is.na(P.80))
# df <- subset(df,CHR1==CHR2)
# df$FDR <- p.adjust(df$P.80,method = 'fdr')
# subset(df,FDR<0.1)
# min(df$FDR)
df$Validate <- as.numeric(sign(df$BETA_INT.80)==sign(df$BETA_INT.20) & df$P.20 < 0.05)

pvector <- df$P.80
o = -log10(sort(pvector, decreasing = FALSE))
e = -log10(ppoints(length(pvector)))
df.qq <- data.frame(exp=e,obs=o,validate=df$Validate)
g1 <- ggplot(df.qq,aes(x=exp,y=obs)) + geom_point() + geom_abline(slope=1,intercept=0,col='red') + 
  theme_bw() + theme(panel.grid=element_blank(),plot.title = element_text(hjust=0.5)) + 
  scale_y_continuous(breaks=seq(0,max(o),by=2)) +
  scale_x_continuous(breaks=seq(0,max(e),by=2)) +
  labs(x = expression(Expected ~ ~-log[10](italic(p))), y = expression(Observed ~ ~-log[10](italic(p)))) #+

g2 <- ggplot(df,aes(x=(BETA_INT.80),y=(BETA_INT.20),col=as.factor(Validate))) + 
  geom_point(size=rel(0.8)) + 
  geom_abline(slope=1,intercept = 0,col='red',lty='dashed') + geom_smooth(method='lm',col='red',se=F) +
  scale_color_manual(values=c('black','red')) + 
  labs(x='Interaction effects (80% discovery set)',y='Interaction effects (20% validation set)') +
   theme_bw() + theme(legend.position = 'none',panel.grid=element_blank())

df.sub <- subset(df,P.80<0.001)
g4 <- ggplot(df.sub,aes(x=(BETA_INT.80),y=(BETA_INT.20),col=as.factor(Validate))) + 
  geom_point(size=rel(1)) + 
  geom_abline(slope=1,intercept = 0,col='red',lty='dashed') + geom_smooth(method='lm',col='red',se=F) +
  scale_color_manual(values=c('black','red')) + 
  labs(x='Interaction effects (80% discovery set)',y='Interaction effects (20% validation set)') +
  theme_bw() + theme(legend.position = 'none',panel.grid=element_blank())

cor.test((df$BETA_INT.80),(df$BETA_INT.20))
cor.test((df.sub$BETA_INT.80),(df.sub$BETA_INT.20))

thres=0.05; suff <- ifelse(is.null(thres),'','PVAL.')
f <- paste0('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxG_2/ukbb.bmi.merged_subset2.GxG.FULL.valid.txt')
df.valid <- fread(f,data.table = F,stringsAsFactors = F)
g3<-ggplot(subset(df.valid,thres >= 1e-4),aes(x=-log10(thres),y=p)) + geom_line() + geom_point() + theme_bw() +
  labs(x=expression(~-log[10](italic(P)) ~ threshold),y='Validation Rate') +
  geom_hline(yintercept=ifelse(is.null(thres),0.5,0.025),linetype='dashed',col='red') +
  geom_ribbon(aes(ymin=lower,ymax=upper),linetype=2,alpha=0.1) + ylim(0,1) +
  theme(panel.grid = element_blank()) + scale_x_continuous(breaks = c(1,2,3,4,5))
# g3

x=2
f <- '~/Documents/Research/vQTL/ukb_vqtl/output/GxG_2/GxG_panel.png'
png(f,width=4000*x,height=4000*x,res=500*x)
plot_grid(g1,g2,g3,g4,ncol=2,rel_widths = c(1,1.4))
dev.off()


# df.sub <- subset(df,P.80 < 0.001)
# ggplot(df.sub,aes(x=abs(BETA_INT.80),y=abs(BETA_INT.20),col=as.factor(Validate))) + geom_point() + 
#   geom_abline(slope=1,intercept = 0,col='red',lty='dashed') + geom_smooth(method='lm',col='red',se=F) +
#   scale_color_manual(values=c('black','red')) + 
#   labs(x='| Interaction effects | (80% discovery set)',y='| Interaction effects | (20% validation set)') +
#   theme_bw() + theme(legend.position = 'none',panel.grid=element_blank())



