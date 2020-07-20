library(data.table)
library(ggplot2)
library(cowplot)
# df <- fread('~/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/diabetes_bmi_marginal.txt',data.table = F,stringsAsFactors = F)
df <- fread('~/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/diabetes_bmi_marginal.used_bmi.txt',data.table = F,stringsAsFactors = F)

g1 <- ggplot(data=subset(df,pheno=='Diabetes' & Set=='80'), aes(x=as.factor(Val), y=EST, ymin=LOW, ymax=HIGH,col=Val,fill=Val)) +
  geom_linerange(size=5,position=position_dodge(width = 0.5)) +
  geom_point(size=3, shape=21, colour="white", fill='black',stroke = 0.5,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=1, lty=2,col='red') +
  coord_flip() +
  labs(x="PA",y="Estimated rs4743930's T allele effect on diabetes risk (OR)") +
  theme_bw()  +
  theme(axis.text = element_text(family = 'Helvetica'),legend.position = 'none') + 
  scale_color_gradient(low='steelblue1',high='darkgreen') +
  scale_x_discrete(labels=c('Low','Medium','High'))

g2 <- ggplot(data=subset(df,pheno=='bmi.na' & Set=='80'), aes(x=as.factor(Val), y=EST, ymin=LOW, ymax=HIGH,col=Val,fill=Val,width=0.1)) +
  geom_linerange(size=5,position=position_dodge(width = 0.5)) +
  geom_point(size=3, shape=21, colour="white", fill='black',stroke = 0.5,position=position_dodge(width = 0.1)) +
  geom_hline(yintercept=0, lty=2,col='red') +
  coord_flip() +
  labs(x="PA",y=expression("Estimated rs4743930's T allele effect on BMI"~(kg/m^2))) +
  theme_bw()  +
  theme(axis.text = element_text(family = 'Helvetica'),legend.position = 'none')+ 
  scale_color_gradient(low='steelblue1',high='darkgreen') +
  # scale_x_discrete(labels=c('Low','Medium','High'),expand = c(0.9,0))
  scale_x_discrete(labels=c('Low','Medium','High'))


f.out<-'~/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/diabetes_bmi_marginal.used_bmi.png'
png(f.out,width = 5400,height=2500,res=700)
plot_grid(g2,g1,nrow = 2)
dev.off()

