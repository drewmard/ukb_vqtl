library(data.table)
library(ggplot2)

f='/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/sig_results/bmi.sig.GxE.txt'
f='~/Documents/Research/vQTL/ukb_vqtl/output/sig_results/bmi.sig.GxE.txt'
results.mg <- fread(f,data.table = F,stringsAsFactors = F)
estimates <- c('BETA.MEAN','BETA.VAR.RAW','BETA.VAR.RINT','dispersion')
envir_factors <- unique(results.mg$E)

cor_res_func <- function(x,k) {
  results.mg.sub <- subset(results.mg,E==envir_factors[k])
  correlation <- cor.test(abs(results.mg.sub$Estimate.x),abs(results.mg.sub[,x]))
  # correlation <- cor.test(abs(results.mg.sub$Estimate.y),abs(results.mg.sub[,x]))
  return(correlation)
}
cor_res_func2 <- function(i,k,cor_res) { 
  data.frame(
    E=envir_factors[k],
    EST=estimates[i],
    COR=as.numeric(cor_res[[i]]$estimate),
    LOW=cor_res[[i]]$conf.int[1],
    HI=cor_res[[i]]$conf.int[2],
    P=cor_res[[i]]$p.value
  )
}
cor_res_func3 <- function(k) {
  cor_res <- lapply(estimates,cor_res_func,k=k)
  cor_res2 <- do.call(rbind,lapply(1:length(estimates),cor_res_func2,k=k,cor_res=cor_res))
  return(cor_res2)
}
df <- do.call(rbind,lapply(1:length(envir_factors),cor_res_func3))

#define colours for dots and bars
dotCOLS = c("#a6d8f0","#f9b282",'steelblue','black')
barCOLS = c("#008fd5","#de6b35",'steelblue4','darkgrey10')
g <- ggplot(data=df, aes(x=E, y=COR, ymin=LOW, ymax=HI,fill=EST,col=EST)) +
  geom_linerange(size=5,position=position_dodge(width = 0.5)) +
  geom_point(size=3, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, lty=2,col='red') +
  # coord_flip() +
  labs(x="Environmental Factor",y="Correlation with GxE effects",col='QTL Type') +
  scale_x_discrete(labels=c('Diet','Alc','PA','SB','Age','Smok','Sex')) +
  theme_bw()  +
  theme(axis.text = element_text(family = 'Helvetica')) +
  scale_fill_manual(values=barCOLS,guide='none')+
  scale_color_manual(values=dotCOLS,labels=c('muQTL','raw vQTL','RINT vQTL','dQTL')) 

x=1.3
# png('~/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/qtl_vs_gxe_effects.png',width=4200*x,height=3100*x,res=400*x)
png('~/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/qtl_vs_gxe_effects.png',width=4200*x,height=1500*x,res=500)
print(g)
dev.off()

g <- ggplot(data=df, aes(x=E, y=COR, ymin=LOW, ymax=HI,fill=EST,col=EST)) +
  geom_linerange(size=5,position=position_dodge(width = 0.5)) +
  geom_point(size=3, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, lty=2,col='red') +
  coord_flip() +
  labs(x="Environmental Factor",y="Correlation with GxE effects",col='QTL Type') +
  scale_x_discrete(labels=c('Diet','Alc','PA','SB','Age','Smok','Sex')) +
  theme_bw()  +
  theme(axis.text = element_text(family = 'Helvetica')) +
  scale_fill_manual(values=barCOLS,guide='none')+
  scale_color_manual(values=dotCOLS,labels=c('muQTL','raw vQTL','RINT vQTL','dQTL')) 

x=1.1
# png('~/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/qtl_vs_gxe_effects.png',width=4200*x,height=3100*x,res=400*x)
png('~/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/qtl_vs_gxe_effects.flip.png',width=3000*x,height=4200*x,res=500)
print(g)
dev.off()
