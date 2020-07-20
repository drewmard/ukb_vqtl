library(data.table)
library(ggplot2)

gxe_effects <- fread('~/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/FTO_marginal_effects.txt',data.table = F,stringsAsFactors = T)
gxe_effects$Val <- ordered(gxe_effects$Val,levels=gxe_effects$Val)
ggCOL <- c(
  'Daily' = 'lightsalmon',
  '3-4x/wk' = 'lightsalmon1',
  '1-2x/wk' = 'lightsalmon2',
  '1-3x/mo' = 'lightsalmon3',
  'Rarely/never' = 'lightsalmon4',
  '37-49'='lightskyblue1',
  '50-59'='lightskyblue3',
  '60-72'='lightskyblue4',
  'Low' = 'tomato1',
  'Mod' = 'tomato2',
  'High' = 'tomato4',
  'No' = 'midnightblue',
  'Yes' = 'mediumvioletred',
  'Female' = 'cadetblue1',
  'Male' = 'cadetblue4',
  '0-2' = 'bisque1',
  '3-4' = 'bisque2',
  '5-6' = 'bisque3',
  '7+' = 'bisque4',
  'Low BMI' = 'steelblue1',
  'Med BMI' = 'steelblue2',
  'High BMI' = 'steelblue3'
)
g <- ggplot(data=gxe_effects, aes(x=Val, y=Estimate, ymin=LOW, ymax=HI,fill=Val,col=Val)) +
  geom_linerange(size=5,position=position_dodge(width = 0.5)) +
  geom_point(size=3, shape=21, colour="white",fill='black', stroke = 0.5,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, lty=2,col='red') +
  coord_flip() +
  labs(x="Environmental Factor",y=expression('Avg BMI increase per rs56094641 G allele ('~kg/m^2~')'),col='QTL Type') +
  theme_bw()  +
  theme(axis.text = element_text(family = 'Helvetica'),legend.position = 'none') +
  facet_grid(E~.,scales='free',space='free') +
  scale_colour_manual(values=ggCOL)

x <- 1
png('~/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/fto_gxe.png',width=5000,height=2500,res=500)
print(g)
dev.off()
