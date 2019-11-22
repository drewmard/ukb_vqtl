library(data.table)
library(ggplot2)

df.tmp <- fread('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/MHC/stats.txt',data.table = F,stringsAsFactors = F,header=F)
colnames(df.tmp)<- as.character(df.tmp[1,]);df.tmp <- df.tmp[-1,]
# df.tmp <-subset(df.tmp,Haplotype %in% 
         c('000','100','001','010','101','111')
         )

library(ggplot2)
ggplot(df.tmp,aes(x=as.factor(Copies),y=as.numeric(Mean),group=as.factor(Haplotype),col=as.factor(Haplotype))) +
  geom_line() + geom_point(aes(size=log10(as.numeric(N)))) + 
  labs(x='Number of haplotypes observed (counts)',y='RINT covariate-adjusted lymphocyte counts (mean)') +
  theme_bw() + scale_color_brewer(palette="Dark2")


df.tmp <- fread('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/MHC/stats2.txt',data.table = F,stringsAsFactors = F,header=F)
colnames(df.tmp)<- as.character(df.tmp[1,]);df.tmp <- df.tmp[-1,]
# df.tmp <-subset(df.tmp,Haplotype %in% 
c('000','100','001','010','101','111')
)

library(ggplot2)
ggplot(df.tmp,aes(x=as.factor(HaplotypeCount),y=as.numeric(LymphocyteCounts),group=as.factor(Haplotype),col=as.factor(Haplotype))) +
  geom_line() + geom_point() + 
  labs(x='Number of haplotypes observed (counts)',y='RINT covariate-adjusted lymphocyte counts (mean)') +
  theme_bw() + scale_color_brewer(palette="Dark2")
