library(data.table)
library(ggplot2)
library(cowplot)
scaleFUN <- function(x) sprintf("%.0f", x)

f <- paste0('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.valid.matched.txt')
df.matched <- fread(f,data.table = F,stringsAsFactors = F); 
df.matched$FDR_based <- 1; df.matched$FDR_based[df.matched$thres %in% c(1,0.1,0.05,0.01,0.005,0.001)] <- 0
g1<-ggplot(subset(df.matched,thres <= 0.05 & N >= 5),aes(x=-log10(thres),y=p)) + geom_line() + geom_point(aes(col=as.factor(FDR_based))) + theme_bw() +
  labs(title='') + scale_color_manual(values=c('black','red')) + 
  geom_ribbon(aes(ymin=lower,ymax=upper),linetype=2,alpha=0.1) + ylim(0,1) + geom_hline(yintercept=0.025,linetype='dashed',col='red') +
  theme(axis.title=element_blank(),panel.grid = element_blank(),legend.position = 'none')+
  scale_x_continuous(breaks=c(1:5)) + xlim(-log10(0.05),5)

f <- paste0('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.valid.QTL.txt')
df.QTL <- fread(f,data.table = F,stringsAsFactors = F); 
df.QTL$FDR_based <- 1; df.QTL$FDR_based[df.QTL$thres %in% c(1,0.1,0.05,0.01,0.005,0.001)] <- 0
g2<-ggplot(subset(df.QTL,thres <= 0.05 & N >= 5),aes(x=-log10(thres),y=p)) + geom_line() + geom_point(aes(col=as.factor(FDR_based))) + theme_bw() +
  labs(title='') + scale_color_manual(values=c('black','red')) + 
  geom_ribbon(aes(ymin=lower,ymax=upper),linetype=2,alpha=0.1) + ylim(0,1) + geom_hline(yintercept=0.025,linetype='dashed',col='red') +
  theme(axis.title=element_blank(),panel.grid = element_blank(),legend.position = 'none')+
  scale_x_continuous(breaks=c(1:5)) + xlim(-log10(0.05),5)

f <- paste0('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.valid.only_mean.txt')
df.mean_only <- fread(f,data.table = F,stringsAsFactors = F); 
df.mean_only$FDR_based <- 1; df.mean_only$FDR_based[df.mean_only$thres %in% c(1,0.1,0.05,0.01,0.005,0.001)] <- 0
g3<-ggplot(subset(df.mean_only,thres <= 0.05 & N >= 5),aes(x=-log10(thres),y=p)) + geom_line() + geom_point(aes(col=as.factor(FDR_based))) + theme_bw() +
  labs(title='') + scale_color_manual(values=c('black','red')) + 
  geom_ribbon(aes(ymin=lower,ymax=upper),linetype=2,alpha=0.1) + ylim(0,1) + geom_hline(yintercept=0.025,linetype='dashed',col='red') +
  theme(axis.title=element_blank(),panel.grid = element_blank(),legend.position = 'none') +
  scale_x_continuous(breaks=c(1:5)) + xlim(-log10(0.05),5)

f <- paste0('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.valid.var_raw.txt')
df.var_raw <- fread(f,data.table = F,stringsAsFactors = F); 
df.var_raw$FDR_based <- 1; df.var_raw$FDR_based[df.var_raw$thres %in% c(1,0.1,0.05,0.01,0.005,0.001)] <- 0
g4<-ggplot(subset(df.var_raw,thres <= 0.05 & N >= 5),aes(x=-log10(thres),y=p)) + geom_line() + geom_point(aes(col=as.factor(FDR_based))) + theme_bw() +
  labs(title='') + scale_color_manual(values=c('black','red')) + 
  geom_ribbon(aes(ymin=lower,ymax=upper),linetype=2,alpha=0.1) + ylim(0,1) + geom_hline(yintercept=0.025,linetype='dashed',col='red') +
  theme(axis.title=element_blank(),panel.grid = element_blank(),legend.position = 'none')+
  scale_x_continuous(breaks=c(1:5)) + xlim(-log10(0.05),5)

x=1
f <- paste0('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/GxE_valid.final_panel.png')
png(f,width = 3000,height=1000,res=300)
plot_grid(g1,g2,g3,g4,nrow=1)
dev.off()

#########################

f <- paste0('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.valid.QTL.null.txt')
df.QTL.null <- fread(f,data.table = F,stringsAsFactors = F); 
df.QTL.null$FDR_based <- 1; df.QTL.null$FDR_based[df.var_raw$thres %in% c(1,0.1,0.05,0.01,0.005,0.001)] <- 0
g5<-ggplot(subset(df.QTL.null,N >= 5),aes(x=-log10(thres),y=p)) + geom_line() + geom_point(aes(col=as.factor(FDR_based))) + theme_bw() +
  labs(title='') + scale_color_manual(values=c('black','red')) + 
  geom_ribbon(aes(ymin=lower,ymax=upper),linetype=2,alpha=0.1) + ylim(0,1) + geom_hline(yintercept=0.5,linetype='dashed',col='red') +
  theme(axis.title=element_blank(),panel.grid = element_blank(),legend.position = 'none')+
  scale_x_continuous(breaks=c(0:5)) + xlim(0,5)

f <- paste0('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.valid.matched.null.txt')
df.matched.null <- fread(f,data.table = F,stringsAsFactors = F); 
df.matched.null$FDR_based <- 1; df.matched.null$FDR_based[df.var_raw$thres %in% c(1,0.1,0.05,0.01,0.005,0.001)] <- 0
g6<-ggplot(subset(df.matched.null,N >= 5),aes(x=-log10(thres),y=p)) + geom_line() + geom_point(aes(col=as.factor(FDR_based))) + theme_bw() +
  labs(title='') + scale_color_manual(values=c('black','red')) + 
  geom_ribbon(aes(ymin=lower,ymax=upper),linetype=2,alpha=0.1) + ylim(0,1) + geom_hline(yintercept=0.5,linetype='dashed',col='red') +
  theme(axis.title=element_blank(),panel.grid = element_blank(),legend.position = 'none')+
  scale_x_continuous(breaks=c(0:5)) + xlim(0,5)


x=1
f <- paste0('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/GxE_valid.null_panel.png')
png(f,width = 2250,height=1500,res=400)
plot_grid(g6,g5,nrow=1)
dev.off()


################

library(data.table)
f <- paste0('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.valid.only_mean_strict.txt')
df.mean_only_strict <- fread(f,data.table = F,stringsAsFactors = F); 
df.mean_only_strict$FDR_based <- 1; df.mean_only_strict$FDR_based[df.mean_only_strict$thres %in% c(1,0.1,0.05,0.01,0.005,0.001)] <- 0
g7<-ggplot(subset(df.mean_only_strict,thres <= 0.05 & N >= 5),aes(x=-log10(thres),y=p)) + geom_line() + geom_point(aes(col=as.factor(FDR_based))) + theme_bw() +
  labs(title='') + scale_color_manual(values=c('black','red')) + 
  geom_ribbon(aes(ymin=lower,ymax=upper),linetype=2,alpha=0.1) + ylim(0,1) + geom_hline(yintercept=0.025,linetype='dashed',col='red') +
  theme(axis.title=element_blank(),panel.grid = element_blank(),legend.position = 'none') +
  scale_x_continuous(breaks=c(1:5)) + xlim(-log10(0.05),5)

x=1
f <- paste0('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/GxE_valid.only_mean_strict.png')
png(f,width = 2250,height=2600,res=700)
plot_grid(g7,nrow=1)
dev.off()


###############3

df.QTL.null$QTL <- 1
df.matched.null$QTL <- 0
df.null=rbind(df.QTL.null,df.matched.null)

g8 <- ggplot(subset(df.null,N >= 5),aes(x=-log10(thres),y=winner,col=as.factor(QTL),fill=as.factor(QTL))) + geom_line() + 
  geom_point(aes(col=as.factor(FDR_based))) + 
  # geom_point() + 
  theme_bw() +
  labs(title='') + scale_color_manual(values=c('black','red')) + 
  ylim(0,1) + geom_hline(yintercept=0.5,linetype='dashed',col='red') +
  theme(axis.title=element_blank(),panel.grid = element_blank(),legend.position = 'none')+
  scale_x_continuous(breaks=c(0:5)) + xlim(0,5)

x=1
f <- paste0('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/GxE_valid.winners_curse.png')
png(f,width = 2250,height=2600,res=600)
plot_grid(g8,nrow=1)
dev.off()

#############

qq_ggplot <- function(pvector,maxy=NULL,labs=TRUE) {
  o = -log10(sort(pvector, decreasing = FALSE))
  e = -log10(ppoints(length(pvector)))
  dataf <- data.frame(e,o)
  if (is.null(maxy)) {maxy=max(o)} 
  if (!labs) {
    g <- ggplot(dataf,aes(x=e,y=o)) + geom_point() + geom_abline(slope=1,intercept=0,col='red') + 
      theme_bw() + theme(panel.grid = element_blank(),axis.title = element_blank()) + ylim(0,maxy) +
      scale_x_continuous(breaks=seq(0,max(-log10(e)),by=1))
  } else {
    g <- ggplot(dataf,aes(x=e,y=o)) + geom_point() + geom_abline(slope=1,intercept=0,col='red') + 
      labs(x=expression(Expected ~ ~-log[10](italic(p))),y=expression(Observed ~ ~-log[10](italic(p)))) +
      theme_bw() + theme(panel.grid = element_blank()) + ylim(0,maxy) +
      scale_x_continuous(breaks=seq(0,max(-log10(e)),by=1))
  }
  return(g)
}

f <- paste0('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/sig_results/bmi.sig.GxE.txt')
GxE <- fread(f,data.table = F,stringsAsFactors = F)
g.qtl <- qq_ggplot(GxE$`Pr(>|t|).x`)

f <- paste0('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/sig_results/bmi.matched_QTL.GxE.txt')
matched <- fread(f,data.table = F,stringsAsFactors = F)
g.matched <- qq_ggplot(matched$`Pr(>|t|).x`,maxy=-log10(min(GxE$`Pr(>|t|).x`)))

f <- '/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.discovery.qqplots.png'
png(f,width=4890,height=2240,res=700)
plot_grid(g.qtl,g.matched,ncol=2)
dev.off()

val <- -log10(min(subset(GxE,Raw.vQTL==1)$`Pr(>|t|).x`))
g.vqtl <- qq_ggplot(subset(GxE,Raw.vQTL==1)$`Pr(>|t|).x`,val)
g.muqtl <- qq_ggplot(subset(GxE,Mean.QTL==1 & Raw.vQTL==0)$`Pr(>|t|).x`,val)
g.rintqtl <- qq_ggplot(subset(GxE,Rint.vQTL==1 & Raw.vQTL==0)$`Pr(>|t|).x`,val)
g.dqtl <- qq_ggplot(subset(GxE,dQTL==1 & Raw.vQTL==0)$`Pr(>|t|).x`,val)

f <- '/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.discovery.qqplots.split.png'
png(f,width=9780,height=2240,res=700)
plot_grid(g.vqtl,g.muqtl,g.rintqtl,g.dqtl,ncol=4)
dev.off()

g.matched <- qq_ggplot(matched$`Pr(>|t|).x`,maxy=-log10(min(GxE$`Pr(>|t|).x`)),labs=F)
g.qtl <- qq_ggplot(GxE$`Pr(>|t|).x`,labs=F)
g.muqtl <- qq_ggplot(subset(GxE,Mean.QTL==1 & Raw.vQTL==0)$`Pr(>|t|).x`,val,labs=F)
g.vqtl <- qq_ggplot(subset(GxE,Raw.vQTL==1)$`Pr(>|t|).x`,val,labs=F)

f <- '/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.discovery.qqplots.main_text.png'
png(f,width=6000,height=2240,res=700)
plot_grid(g.matched,g.qtl,g.muqtl,g.vqtl,ncol=4)
dev.off()

x=1
f <- paste0('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.discovery_and_validation.final.png')
png(f,width = 3000,height=2000,res=300)
plot_grid(g.matched,g.qtl,g.muqtl,g.vqtl,
          g1,g2,g3,g4,
          nrow=2,ncol=4)
dev.off()




