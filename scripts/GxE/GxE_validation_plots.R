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






