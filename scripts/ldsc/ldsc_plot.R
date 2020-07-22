library(data.table)
library(ggplot2)

f <- '/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/ldsc/mean.sumSS.Multi_tissue_gene_expr.coef1.cell_type_results.txt'
# f <- '/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/ldsc/mean.sumSS.Multi_tissue_chromatin.coef1.cell_type_results.txt'
mean <- fread(f,data.table = F,stringsAsFactors = F)
colnames(mean)[2:4] <- paste0(colnames(mean)[2:4],'.MEAN')
# mean[order(mean$Coefficient_P_value),][1:5,]
mean$FDR.MEAN <- p.adjust(mean[,4],method = 'fdr')

f <- '/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/ldsc/var.sumSS.Multi_tissue_gene_expr.coef1.cell_type_results.txt'
# f <- '/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/ldsc/var.sumSS.Multi_tissue_chromatin.coef1.cell_type_results.txt'
var <- fread(f,data.table = F,stringsAsFactors = F)
colnames(var)[2:4] <- paste0(colnames(var)[2:4],'.VAR')
var$FDR.VAR <- p.adjust(var[,4],method = 'fdr')

df.mg <- merge(mean,var,by='Name')
df.mg$Stomach <- 0
df.mg$Stomach[df.mg$Name=='A03.556.875.875.Stomach'] <- 1
df.mg$Stomach[df.mg$Name=='Stomach'] <- 2

cor.res<-cor.test(df.mg$Coefficient.MEAN,df.mg$Coefficient.VAR)
cor.res$estimate; cor.res$p.value

t.test(abs(df.mg$Coefficient.MEAN),abs(df.mg$Coefficient.VAR),paired = T)$p.value

g1 <- ggplot(df.mg,aes(x=Coefficient.MEAN,y=Coefficient.VAR)) + 
  geom_point(aes(col=as.factor(Stomach))) + 
  geom_abline(slope=1,intercept=0,col='red',lty='dashed') +
  geom_smooth(method='lm',se=F,col='blue') +
  scale_color_manual(values=c('black','red','orange')) +
  theme_bw() +
  theme(legend.position = 'none',panel.grid = element_blank()) +
  labs(x='Cell-type enrichment of mean heritability',y='Cell-type enrichment of variance heritability')

library(reshape2)
df.mg.melt <- melt(df.mg[,c('Name','Coefficient.MEAN','Coefficient.VAR')],id.vars='Name')
g2 <- ggplot(df.mg.melt,aes(x=variable,y=abs(value),fill=variable)) + 
  geom_boxplot(outlier.shape=NA,col='black') +
  geom_jitter(col='orange',width=0.1,alpha=0.3) +
  labs(x='Heritability',y='| Cell-type enrichment |') + 
  theme_bw() +
  theme(legend.position = 'none',panel.grid = element_blank()) +
  scale_fill_manual(values=c('azure1','azure4')) +
  scale_x_discrete(labels=c('Mean','Variance'))

library(cowplot)
png('~/Documents/Research/vQTL/ukb_vqtl/output/ldsc/mean_vs_var_plot.png',width = 4500,height=2500,res=520)
plot_grid(g1,g2,ncol=2)
dev.off()
subset(df.mg,FDR.VAR<0.1)
df.mg[grep('Stomach',df.mg$Name),]
subset(df.mg,FDR.MEAN < 0.1 | FDR.VAR < 0.1)

library(ggplot2)
library(reshape2)
library(forcats)
library(dplyr)

###############

tmp <- melt(subset(df.mg,FDR.VAR < 0.1)[,c('Name','FDR.MEAN','FDR.VAR')],by='Name')
tmp$log10 <- -log10(tmp$value)
tmp$log10[tmp$variable=='FDR.VAR'] <- -1*tmp$log10[tmp$variable=='FDR.VAR']
# ggplot(tmp,aes(x=Name,y=log10,fill=variable)) + 

tmp$Name <-  as.factor(tmp$Name)
tmp <- tmp %>%
  mutate(name =
           fct_reorder(tmp$Name,tmp$log10,max)
  )

tmp$Col2 <- 'yellow1'
tmp$Col2[tmp$Name!='A03.556.875.875.Stomach' & tmp$variable=='FDR.VAR'] <- 'yellow3'
tmp$Col2[tmp$Name=='A03.556.875.875.Stomach' & tmp$variable=='FDR.MEAN'] <- 'coral1'
tmp$Col2[tmp$Name=='A03.556.875.875.Stomach' & tmp$variable=='FDR.VAR'] <- 'coral3'
tmp$Col2 <- as.factor(tmp$Col2)

g <- ggplot(tmp,
            aes(x=name,
                y=log10,
                fill=Col2)) +
  geom_bar(stat='identity',position='identity') +
  theme_bw() + 
  theme(
    axis.text.x=element_blank(),
    panel.grid=element_blank(),
    axis.title = element_blank(),
    legend.title = element_blank()
  ) +
  geom_abline(slope=0,intercept=-log10(0.1),col='red',lty='dashed') +
  geom_abline(slope=0,intercept=log10(0.1),col='red',lty='dashed') +
  geom_abline(slope=0,intercept=0,col='black') +
  labs(y='-log10 P-values of FDR',x='Tissue') +
  scale_fill_manual(values=c('coral1','coral3','yellow2','yellow3'),labels=c('Stomach','Stomach','CNS','CNS')) +
  scale_x_discrete(labels=ifelse(substring(levels(tmp$name),1,1)=='A','F','G')) +
  scale_y_continuous(labels=c(2,0,2,4,6))
g

library(cowplot)
x=2
png('~/Documents/Research/vQTL/ukb_vqtl/output/ldsc/enriched_annotations2.png',width = 950*x,height=1500*x,res=300*x)
g
dev.off()

##########

# another version of plot

# tmp <- melt(subset(df.mg,FDR.MEAN < 0.1 | FDR.VAR < 0.1)[,c('Name','FDR.MEAN','FDR.VAR')],by='Name')
# # tmp <- melt(subset(df.mg,FDR.VAR < 0.1)[,c('Name','FDR.MEAN','FDR.VAR')],by='Name')
# tmp$log10 <- -log10(tmp$value)
# tmp$log10[tmp$variable=='FDR.VAR'] <- -1*tmp$log10[tmp$variable=='FDR.VAR']
# 
# tmp$Name <-  as.factor(tmp$Name)
# tmp <- tmp %>%
#   mutate(name =
#            fct_reorder(tmp$Name,tmp$log10,max)
#   )
# 
# tmp$Col2 <- 'yellow1'
# tmp$Col2[tmp$Name!='A03.556.875.875.Stomach' & tmp$variable=='FDR.VAR'] <- 'yellow3'
# tmp$Col2[tmp$Name=='A03.556.875.875.Stomach' & tmp$variable=='FDR.MEAN'] <- 'coral1'
# tmp$Col2[tmp$Name=='A03.556.875.875.Stomach' & tmp$variable=='FDR.VAR'] <- 'coral3'
# tmp$Col2 <- as.factor(tmp$Col2)
# 
# g <- ggplot(tmp,
#             aes(x=name,
#                 y=log10,
#                 fill=Col2)) +
#   # fill=variable)) +
#   geom_bar(stat='identity',position='identity') +
#   # ylim(min(tmp$log10[tmp$variable=='FDR.VAR']),max(tmp$log10[tmp$variable=='FDR.MEAN'])) +
#   theme_bw() + 
#   theme(#axis.text.x=element_blank(),
#     # axis.text.x=element_text(hjust=1,angle = 80),
#     axis.text.y=element_text(colour=ifelse(substring(levels(tmp$name),1,1)=='A','yellow4','coral4')),
#     panel.grid=element_blank(),
#     axis.title = element_blank(),
#     legend.title = element_blank()
#     # legend.position = 'none',
#   ) +
#   geom_abline(slope=0,intercept=-log10(0.1),col='red',lty='dashed') +
#   geom_abline(slope=0,intercept=log10(0.1),col='red',lty='dashed') +
#   geom_abline(slope=0,intercept=0,col='black') +
#   labs(y='-log10 P-values of FDR',x='Tissue') +
#   scale_fill_manual(values=c('coral1','coral3','yellow2','yellow3'),labels=c('Stomach','Stomach','CNS','CNS')) +
#   scale_x_discrete(labels=ifelse(substring(levels(tmp$name),1,1)=='A','F','G')) +
#   scale_y_continuous(labels=c(2,0,2,4,6)) +
#   coord_flip(); g
# 
# library(cowplot)
# x=2
# png('~/Documents/Research/vQTL/ukb_vqtl/output/ldsc/enriched_annotations.png',width = 1100*x,height=1500*x,res=300*x)
# g
# dev.off()

