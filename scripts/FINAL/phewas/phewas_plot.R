library(data.table)
library(ggplot2)
library(cowplot)
library(ggplot2)

f <- "/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/phewas/phewas_muqtl_vs_vqtl.full.txt"
df <- fread(f,data.table = F,stringsAsFactors = F)
df$Sig <- as.numeric(df$FDR < 0.1)
g1<-ggplot(df,aes(x=Expected,y=Observed,col=as.factor(Sig))) + geom_point() + geom_abline(slope=1,intercept = 0,col='red',linetype='dashed') + 
  theme_bw() +   
  theme(legend.position = "none",panel.grid = element_blank()) + scale_color_manual(values=c('black','red')) + lims(x=c(0,1),y=c(0,1)) +
  labs(x='Proportion of pure muQTLs associated',y='Proportion of vQTLs associated')

f <- "/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/phewas/phewas_muqtl_vs_muqtl.full.txt"
df <- fread(f,data.table = F,stringsAsFactors = F)
df$Sig <- as.numeric(df$FDR < 0.1)
g2<-ggplot(df,aes(x=Expected,y=Observed,col=as.factor(Sig))) + geom_point() + geom_abline(slope=1,intercept = 0,col='red',linetype='dashed') + 
  theme_bw() +   theme(legend.position = "none",panel.grid = element_blank()) + scale_color_manual(values=c('black','red')) + lims(x=c(0,1),y=c(0,1)) +
  labs(x='Proportion of 21 sampled pure muQTLs associated',y='Proportion of other pure muQTLs associated')

png('~/Documents/Research/vQTL/ukb_vqtl/output/phewas/plots.png',width=4000,height=2000,res=500)
plot_grid(g1,g2,nrow=1)
dev.off()
