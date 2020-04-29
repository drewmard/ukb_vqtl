library(data.table)
df <- list()
workdir <- '/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/simulation/'
MAF1=0.4
nindiv=10000
f <- paste0(workdir,'MAF1_',MAF1,'.MAF2_0.4.NSIM_1000.NINDIV_',nindiv,'.TYPE_gxg.NOISE_NORMAL.txt')
df[[1]] <- fread(f,data.table = F,stringsAsFactors = F)
nindiv=20000
workdir <- '/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/simulation/'
f <- paste0(workdir,'MAF1_',MAF1,'.MAF2_0.4.NSIM_1000.NINDIV_',nindiv,'.TYPE_gxg.NOISE_NORMAL.txt')
df[[2]] <- fread(f,data.table = F,stringsAsFactors = F)

df <- do.call(rbind,df)

df.aggre <- aggregate(.~MAF1+MAF2+N+noise+h+type,df,function(x) {mean(x < 0.05)})
mean(subset(df.aggre,N==10000 & h <= 0.05)$DRM) / mean(subset(df.aggre,N==10000 & h <= 0.05)$LT)

# df.aggre <- aggregate(.~MAF1+MAF2+N+noise+h+type,df,mean)

ggplot(subset(df.aggre,h<=0.03),aes(h,N)) + geom_tile(aes(fill=DRM),col='white') +
  scale_fill_gradient(low = "white", high = "steelblue",name='Power') +
  theme_bw()  + theme(panel.grid = element_blank()) + #,legend.title = element_blank()) + 
  labs(x='Variance Explained by GxG',y='Sample Size') + scale_y_continuous(breaks=c(10000,20000))

