

library(data.table)
library(ggplot2)
# f<-'/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.ALL_RESULTS.trim.txt'
f<-'/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/sig_results/bmi.sig.GxE.txt'
df <- fread(f,data.table = F,stringsAsFactors = F)
colnames(df)[4] <- 'P_GxE'
df$FDR <- p.adjust(df$P_GxE,method = 'fdr')
df$Sig <- 0
df$Sig[df$FDR < 0.1] <- 1
df$Sig[df$FDR < 0.05] <- 2
df$Sig[df$FDR < 0.01] <- 3
df <- subset(df,SNP %in% subset(df,Sig > 0)$SNP)

axis_colors <- c('black','orange')[as.numeric(subset(df,!duplicated(df$SNP))$P.VAR.RAW < 5e-8) + 1]
g <- ggplot(df,aes(x=SNP,y=E,fill=Sig)) +
  geom_tile(col='black') +
  theme_bw()  +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, hjust=1,colour = axis_colors),
        panel.border = element_blank(),
        legend.position = 'none'
  ) +
  scale_fill_gradient2(low = "blue", mid='white',high = "red",name='Power') +
  scale_y_discrete(labels=c('Smoking.E'='Smok','sex'='Sex','SB'='SB','PA'='PA','DIET_SCORE'='Diet','Alcohol_intake_frequency'='Alc','age'='Age')) +
  labs(x='SNP',y='Environmental Factor')

x=1
png('~/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/results_heatmap.png',width=4000*x,height=1100*x,res=320*x)
g
dev.off()

