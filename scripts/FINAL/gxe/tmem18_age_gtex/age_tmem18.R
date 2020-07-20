library(data.table)
library(ggplot2)
f <- '/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/Age_GTEx.txt'
df <- fread(f,data.table = F,stringsAsFactors = F)

g <- ggplot(df,aes(x=AGE,y=TMEM18,fill=AGE)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width=0.1,size=rel(0.8)) + scale_fill_brewer(palette=5) +
  theme_bw() + theme(panel.grid = element_blank())

x=1
f <- '/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/adipose_age.png'
png(f,width = 2500*x,height=2500*x,res=500*x)
print(g)
dev.off()
