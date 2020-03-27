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
