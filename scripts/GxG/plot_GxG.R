library(data.table)
library(ggplot2)

df <- fread('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/GxG/plot/GxG_data.txt',data.table = F,stringsAsFactors = F)
SNP1 <- colnames(df)[1]
SNP2 <- colnames(df)[2]
phenoName <- colnames(df)[3]
df.sub <- df
df.sub <- subset(df,!(is.na(df.sub[,SNP1])) & !(is.na(df.sub[,SNP2])))
# df.sub <- df[1:5000,]

x <-paste0('lymphocyte.count.rint.na','~',SNP2)
x <- formula(x)
mod0 <- lm(x,data=subset(df.sub,df.sub[,SNP1]==0))
mod1 <- lm(x,data=subset(df.sub,df.sub[,SNP1]==1))
mod2 <- lm(x,data=subset(df.sub,df.sub[,SNP1]==2))

ggplot(df.sub,aes(as.numeric(df.sub[,SNP2]),df.sub[,phenoName],col=as.factor(df.sub[,SNP1]))) +
  geom_jitter(position=position_jitter(width=.1, height=0),alpha=0.1) +
  geom_abline(slope=mod0$coefficients[2],intercept = mod0$coefficients[1],col='red') + 
  geom_abline(slope=mod1$coefficients[2],intercept = mod1$coefficients[1],col='green') + 
  geom_abline(slope=mod2$coefficients[2],intercept = mod2$coefficients[1],col='blue') +
  labs(x=SNP2,y=paste0(phenoName)) + scale_color_discrete(name=SNP1)

ggplot(df.sub,aes(as.factor(df.sub[,SNP2]),df.sub[,phenoName],col=as.factor(df.sub[,SNP1]))) +
  geom_boxplot(outlier.shape = NA) +
  labs(x=SNP2,y=paste0(phenoName)) + scale_color_discrete(name=SNP1)


ggplot(df.sub,aes(as.factor(df.sub[,SNP1]),df.sub[,phenoName],col=as.factor(df.sub[,SNP2]))) +
  geom_boxplot(outlier.shape = NA) +
  labs(x=SNP1,y=paste0(phenoName)) + scale_color_discrete(name=SNP2)

aggregate(df.sub$lymphocyte.count.rint.na,by=list(df.sub$rs887468),var,na.rm=T)
aggregate(df.sub$lymphocyte.count.rint.na,by=list(df.sub$rs2516491),var,na.rm=T)
