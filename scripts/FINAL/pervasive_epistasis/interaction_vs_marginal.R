library(ggplot2)
library(data.table)
library(cowplot)
# Initialize
Nindiv=10000
MAF=0.25
# Mcausalsnp=50; Mtotsnp=1000
Mcausalsnp=400; Mtotsnp=1000
Ninteraction <- Mtotsnp * (Mtotsnp-1) / 2

f <- paste0('~/Documents/Research/vQTL/ukb_vqtl/output/sim/epistatic_variance_simulations_dataset.N_',Nindiv,'.MAF_',MAF,'.Mcausal_',Mcausalsnp,'.Mtot_',Mtotsnp,'.txt')
gwas <- fread(f,data.table = F,stringsAsFactors = F)

for (i.to_subset in 1:1000) {
  print(i.to_subset)
  gwas.sub <- subset(gwas,i==i.to_subset | j==i.to_subset)
  # gwas.sub <- subset(gwas,(i==i.to_subset | j==i.to_subset) & (causal.i==1 & causal.j==1))
  df.tmp <- data.frame(i=i.to_subset,true=mean(gwas.sub$true_beta),est=mean(gwas.sub$BETA.GxG))
  if (i.to_subset==1) {
    df.save <- df.tmp
  } else {
    df.save <- rbind(df.save,df.tmp)
  }
}
f <- paste0('~/Documents/Research/vQTL/ukb_vqtl/output/sim/epistatic_variance_simulations_dataset.N_',Nindiv,'.MAF_',MAF,'.Mcausal_',Mcausalsnp,'.Mtot_',Mtotsnp,'.mean_var.txt')
df2 <- fread(f,data.table = F,stringsAsFactors = F)
df.mg <- cbind(df.save,df2)
df.mg.sub <- subset(df.mg,causal==1); 
#cor.test(df.mg.sub$BETA.MEAN,df.mg.sub$BETA.VAR)

# print('Mean-var plots...')
# gA <- ggplot(df.mg.sub,aes(x=BETA.MEAN,y=BETA.VAR)) + geom_point() + geom_smooth(method='lm',col='red',se=F) + theme_bw() + theme(panel.grid = element_blank(),legend.position = 'none') + labs(x='Mean effects',y='Variance effects')
# gB <- ggplot(df.mg.sub,aes(x=-log10(P.MEAN),y=-log10(P.VAR))) + geom_point() + geom_smooth(method='lm',col='red',se=F) + theme_bw() + theme(panel.grid = element_blank(),legend.position = 'none') + labs(x='Mean effects',y='Variance effects')
# x=5
# f <- paste0('~/Documents/Research/vQTL/ukb_vqtl/output/sim/epistatic_variance_simulations_dataset.N_',Nindiv,'.MAF_',MAF,'.Mcausal_',Mcausalsnp,'.Mtot_',Mtotsnp,'.mean_var.png')
# png(f,width = 1700*x,height=2000*x,res = 300*x)
# # plot_grid(gA,gB,ncol = 2)
# gA
# dev.off()

# mean(df.mg$P.MEAN < 0.05/1000)
# mean(df.mg$P.MEAN[1:Mcausalsnp] < 0.05/Mtotsnp)
# mean(df.mg$P.MEAN[-c(1:Mcausalsnp)] < 0.05/Mtotsnp)
# min(df.mg$P.MEAN)
# 
# cor.test(df.mg$true,df.mg$est)
# cor.test(df.mg$true,df.mg$BETA.MEAN)
# cor.test(df.mg$est,df.mg$BETA.MEAN)
# df.mg.sub <- subset(df.mg,causal==1)
# cor.test(df.mg.sub$true,df.mg.sub$est)
# cor.test(df.mg.sub$true,df.mg.sub$BETA.MEAN)
# cor.test(df.mg.sub$est,df.mg.sub$BETA.MEAN)
# cor.test(gwas.sub$BETA.GxG,gwas.sub$true_beta)

# cor.test(df.mg.sub$BETA.MEAN,df.mg.sub$BETA.VAR)

# gD <- ggplot(df.mg,aes(x=true,y=BETA.MEAN))+geom_point(aes(col=as.factor(causal))) + theme_bw() + theme(legend.position = 'none',panel.grid = element_blank()) + labs(x='Mean "true" interaction effect size',y='Observed marginal SNP effects') + geom_smooth(method='lm',col='red',se=F)
# gE <- ggplot(df.mg,aes(x=est,y=BETA.MEAN))+geom_point(aes(col=as.factor(causal))) + theme_bw() + theme(legend.position = 'none',panel.grid = element_blank()) + labs(x='Mean estimated interaction effect size',y='Observed marginal SNP effects') + geom_smooth(method='lm',col='red',se=F)
# gF <- ggplot(df.mg,aes(x=true,y=est))+geom_point(aes(col=as.factor(causal))) + theme_bw() + theme(legend.position = 'none',panel.grid = element_blank()) + labs(x='Mean "true" interaction effect size',y='Mean estimated interaction effect size') + geom_smooth(method='lm',col='red',se=F)
gD <- ggplot(df.mg.sub,aes(x=true,y=BETA.MEAN))+geom_point() + theme_bw() + theme(legend.position = 'none',panel.grid = element_blank()) + labs(x='Mean "true" interaction effect size',y='Observed marginal SNP effects') + geom_smooth(method='lm',col='red',se=F)
gE <- ggplot(df.mg.sub,aes(x=est,y=BETA.MEAN))+geom_point() + theme_bw() + theme(legend.position = 'none',panel.grid = element_blank()) + labs(x='Mean estimated interaction effect size',y='Observed marginal SNP effects') + geom_smooth(method='lm',col='red',se=F)
gF <- ggplot(df.mg.sub,aes(x=true,y=est))+geom_point() + theme_bw() + theme(legend.position = 'none',panel.grid = element_blank()) + labs(x='Mean "true" interaction effect size',y='Mean estimated interaction effect size') + geom_smooth(method='lm',col='red',se=F)
# plot_grid(gD,gF,gE,ncol=3)

gwas.sub <- subset(gwas,causal.i==1 & causal.j==1)
gC <- ggplot(gwas.sub,aes(x=BETA.GxG,y=true_beta)) + geom_point(size=rel(0.8)) + theme_bw() + theme(legend.position = 'none',panel.grid = element_blank()) + labs(x='Estimated interaction effects',y='True interaction effects') + geom_smooth(method='lm',se=F,col='red')
# sum(gwas$P.GxG < .05/nrow(gwas.sub))/nrow(gwas.sub)
# f <- paste0('~/Documents/Research/vQTL/ukb_vqtl/output/sim/epistatic_variance_simulations_dataset.N_',Nindiv,'.MAF_',MAF,'.Mcausal_',Mcausalsnp,'.Mtot_',Mtotsnp,'.est_vs_true_int_eff.png')
# png(f,width = 2000*x,height=2000*x,res = 300*x)
# gC
# dev.off()

x=5
f <- paste0('~/Documents/Research/vQTL/ukb_vqtl/output/sim/epistatic_variance_simulations_dataset.N_',Nindiv,'.MAF_',MAF,'.Mcausal_',Mcausalsnp,'.Mtot_',Mtotsnp,'.interaction_vs_marginal.png')
png(f,width = 3000*x,height=1000*x,res = 300*x)
plot_grid(gD,gE,gC,ncol = 3)
# plot_grid(gD,gF,gE,ncol=3)
dev.off()


# table(sign(gwas.sub$BETA.GxG))
# table(sign(gwas.sub$true_beta))
# mean(gwas.sub$true_beta)
# mean(gwas.sub$BETA.GxG)
