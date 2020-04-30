library(data.table)
df <- list()
workdir <- '/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/output/simulation/'
MAF1=0.4

for (i in 1:5) {
  nindiv=i*10000
  f <- paste0(workdir,'MAF1_',MAF1,'.MAF2_0.4.NSIM_1000.NINDIV_',nindiv,'.TYPE_gxg.NOISE_NORMAL.txt')
  df[[i]] <- fread(f,data.table = F,stringsAsFactors = F)
}
df <- do.call(rbind,df)

df.aggre <- aggregate(.~MAF1+MAF2+N+noise+h+type,df,function(x) {mean(x < 0.05)})
mean(subset(df.aggre,N==10000 & h <= 0.05)$DRM) / mean(subset(df.aggre,N==10000 & h <= 0.05)$LT)
mean(subset(df.aggre,N==nindiv & h <= 0.03)$DRM) / mean(subset(df.aggre,N==nindiv & h <= 0.03)$LT)


ggplot(subset(df.aggre,h<=0.03),aes(h,N)) + geom_tile(aes(fill=DRM),col='white') +
  geom_text(aes(label=DRM)) +
  scale_fill_gradient(low = "white", high = "steelblue",name='Power') +
  theme_bw()  + theme(panel.grid = element_blank()) + #,legend.title = element_blank()) + 
  labs(x='Variance Explained by GxG',y='Sample Size') + scale_y_continuous(breaks=seq(min(df.aggre$N),max(df.aggre$N),by=10000))

df.aggre <- aggregate(.~MAF1+MAF2+N+noise+h+type,df,function(x) {mean(x < 0.05)})
ggplot(subset(df.aggre,h<=0.03),aes(N,h)) + geom_tile(aes(fill=DRM),col='white') +
  geom_text(aes(label=DRM)) +
  scale_fill_gradient(low = "white", high = "steelblue",name='Power') +
  theme_bw()  + theme(panel.grid = element_blank()) + #,legend.title = element_blank()) + 
  labs(y='Variance Explained by GxG',x='Sample Size') + scale_x_continuous(breaks=seq(min(df.aggre$N),max(df.aggre$N),by=10000))

df.aggre <- aggregate(.~MAF1+MAF2+N+noise+h+type,df,mean)
ggplot(subset(df,h<=0.03),aes(x=as.factor(h),y=beta,fill=as.factor(N))) + geom_boxplot(alpha=0.5) + 
  theme_bw() + theme(panel.grid = element_blank()) +#,legend.position = 'none') + 
  labs(x='Variance Explained by GxG',y=expression(beta['var']),fill='Sample Size')
  # +
  # scale_color_continuous(low='orange',high='blue')
  scale_fill_continuous(low='yellow',high='red')
