library(data.table)
library(ggplot2)

f <- '/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxG_2/GxG_correlation_between_disc_and_valid.txt'
df <- fread(f,data.table = F,stringsAsFactors = T)
levels(df$SUB) <- df$SUB
df$i <- 1:nrow(df)
g <- ggplot(data=df, aes(x=SUB, y=cor, ymin=LOW, ymax=HI,col=SUB,fill=SUB)) +
  geom_linerange(size=5,position=position_dodge(width = 0.5)) +
  geom_point(size=3, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, lty=2,col='red') +
  coord_flip() +
  labs(x="Subset",y="Correlation of GxG effects (Discovery vs Replication)",col='Subset') +
  theme_bw()  +
  theme(axis.text = element_text(family = 'Helvetica'),legend.position = 'none')

f.out <- '/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/GxG_2/GxG_correlation_between_disc_and_valid.png'
png(f.out,width = 5000,height=1800,res=500)
g
dev.off()
