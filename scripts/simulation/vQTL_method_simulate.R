library(car)
library(lmtest)
library(dglm)

source('/Users/andrewmarderstein/Documents/Research/vQTL/ukb_vqtl/bin/rntransform.R')

nindiv <- 100000
MAF1 <- 0.4
MAF2 <- 0.4
nsim <- 50; 
simMethod=2
v.epi.exp.vec <- seq(0.005,0.015,by=0.005); j = 1
# v.epi.exp.vec <- seq(0.005,0.05,by=0.005); j = 1

transform = 'RINT' # 'NONE'
phenotype.vec <- c('NORMAL','CHISQ4'); phenotype <- phenotype.vec[1]
CHISQ_df <- 15

DeviationRegressionModel <- function(Y,x) {
  Y.i <- tapply(Y, as.factor(x), median)
  Z.ij <- abs(Y - Y.i[as.factor(x)])
  p <- summary(lm(Z.ij~x))$coef[2,]
  return(p)
}

p.mat <- rep(list(as.data.frame(matrix(NA,nrow=nsim*length(v.epi.exp.vec),ncol=9))),2)
p.mat[[1]][,9] <- rep(v.epi.exp.vec,each=nsim)
p.mat[[2]][,9] <- rep(v.epi.exp.vec,each=nsim)

p.mat2 <- rep(list(as.data.frame(matrix(NA,nrow=nsim*length(v.epi.exp.vec),ncol=9))),2)
p.mat2[[1]][,9] <- rep(v.epi.exp.vec,each=nsim)
p.mat2[[2]][,9] <- rep(v.epi.exp.vec,each=nsim)

for (phenotype in 1:length(phenotype.vec)) {
  for (simMethod in 1:3) {
    p.mat[[1]]$phenotype <- phenotype
    p.mat[[1]]$simMethod <- simMethod
    p.mat[[2]]$phenotype <- phenotype
    p.mat[[2]]$simMethod <- simMethod
    p.mat2[[1]]$phenotype <- phenotype
    p.mat2[[1]]$simMethod <- simMethod
    p.mat2[[2]]$phenotype <- phenotype
    p.mat2[[2]]$simMethod <- simMethod
    
  }
}
for (j in 1:length(v.epi.exp.vec)) { 
  print(paste0(j,'/',length(v.epi.exp.vec),': ',v.epi.exp.vec[j]))
  for (i in 1:nsim) {
    
    print(i)
    G <- cbind(rbinom(nindiv,2,MAF1),rbinom(nindiv,2,MAF2))
    Y.epi <- G[,1]*G[,2]
    v.epi.exp <- v.epi.exp.vec[j]
    scalar <- (1-(v.epi.exp))/(v.epi.exp); var.env <- scalar*var(Y.epi)

    if (simMethod==1) {
      
      if (phenotype=='NORMAL') {
        Y2 <- Y.epi + rnorm(nindiv,0,sqrt(var.env))
      } else if (phenotype=='CHISQ4') {
        # Y.env <- as.numeric(sqrt(var.env)*scale(rchisq(nindiv,df=4)))
        # Y <- Y.epi + Y.env
        Y2 <- rchisq(nindiv,(Y.epi+var.env)/4)
      } 
    } else if (simMethod==2) {
      Y <- Y.epi + rnorm(nindiv,0,sqrt(var.env))
      if (phenotype=='NORMAL') {
        Y2.tmp <- qnorm((1:length(Y)) / (length(Y)+1),0,1)
      } else if (phenotype=='CHISQ4') {
        Y2.tmp <- qchisq((1:length(Y)) / (length(Y)+1),4)
      }
      Y2 <- Y2.tmp[match(rank(Y),rank(Y2.tmp))]
    } else if (simMethod==3) {
      if (phenotype=='NORMAL') {
        Y2 <- Y.epi + rnorm(nindiv,0,sqrt(var.env))
      } else if (phenotype=='CHISQ4') {
        Y.env <- as.numeric(sqrt(var.env)*scale(rchisq(nindiv,df=4)))
        Y2 <- Y.epi + Y.env
      }
    }
    Y2 <- scale(Y2)[,1]
    if (transform=='RINT') {
      Y2.RINT <- rntransform(Y2)
      # Y2.RINT <- log10(Y2-min(Y2)+1)
    }
    
    # G <- cbind(rbinom(nindiv,2,MAF1),rbinom(nindiv,2,MAF2))
    
    res <- as.numeric(DeviationRegressionModel(Y2,G[,1]))
    beta.VAR.GxG <- res[1]
    p.VAR.GxG <- res[4]
    
    res <- summary(lm(Y2~G[,1]))$coef[2,]
    beta.MEAN.GxG <- res[1]
    p.MEAN.GxG <- res[4]
    
    res <- as.numeric(DeviationRegressionModel(Y2.RINT,G[,1]))
    beta.VAR.GxG.RINT <- res[1]
    p.VAR.GxG.RINT <- res[4]
    
    res <- summary(lm(Y2.RINT~G[,1]))$coef[2,]
    beta.MEAN.GxG.RINT <- res[1]
    p.MEAN.GxG.RINT <- res[4]
    
    p.mat[[1]][(nsim*(j-1))+i,1:8] <- c(beta.MEAN.GxG,beta.MEAN.GxG.RINT,beta.VAR.GxG,beta.VAR.GxG.RINT,
                           p.MEAN.GxG,p.MEAN.GxG.RINT,p.VAR.GxG,p.VAR.GxG.RINT)
    
    x <- G[,1]
    p.mat2[[1]][(nsim*(j-1))+i,1] <- as.numeric(DeviationRegressionModel(Y2,G[,1]))[4]
    p.mat2[[1]][(nsim*(j-1))+i,2] <- leveneTest(Y2~as.factor(x),center=median)$"Pr(>F)"[1]
    p.mat2[[1]][(nsim*(j-1))+i,3] <- as.numeric(bptest(Y2~x)$p.value)
    p.mat2[[1]][(nsim*(j-1))+i,4] <- coef(summary(dglm(Y2~x))$dispersion.summary)[1,4]
    
    p.mat2[[1]][(nsim*(j-1))+i,5] <- as.numeric(DeviationRegressionModel(Y2.RINT,G[,1]))[4]
    p.mat2[[1]][(nsim*(j-1))+i,6] <- leveneTest(Y2.RINT~as.factor(x),center=median)$"Pr(>F)"[1]
    p.mat2[[1]][(nsim*(j-1))+i,7] <- as.numeric(bptest(Y2.RINT~x)$p.value)
    p.mat2[[1]][(nsim*(j-1))+i,8] <- coef(summary(dglm(Y2.RINT~x))$dispersion.summary)[1,4]
    
    
    ####### Additive
    
    Y.epi <- G[,1]
    v.epi.exp <- v.epi.exp.vec[j]
    scalar <- (1-(v.epi.exp))/(v.epi.exp); var.env <- scalar*var(Y.epi)
    
    if (phenotype=='NORMAL') {
      Y2 <- Y.epi + rnorm(nindiv,0,sqrt(var.env))
    } else if (phenotype=='CHISQ4') {
      # Y.env <- as.numeric(sqrt(var.env)*scale(rchisq(nindiv,df=4)))
      # Y <- Y.epi + Y.env
      Y2 <- rchisq(nindiv,(Y.epi+var.env)/4)
    }
    
    Y2 <- scale(Y2)[,1]
    if (transform=='RINT') {
      Y2.RINT <- rntransform(Y2)
      # Y2.RINT <- log10(Y2-min(Y2)+1)
    }
    
    res <- as.numeric(DeviationRegressionModel(Y2,G[,1]))
    beta.VAR.G <- res[1]
    p.VAR.G <- res[4]
    
    res <- summary(lm(Y2~G[,1]))$coef[2,]
    beta.MEAN.G <- res[1]
    p.MEAN.G <- res[4]
    
    res <- as.numeric(DeviationRegressionModel(Y2.RINT,G[,1]))
    beta.VAR.G.RINT <- res[1]
    p.VAR.G.RINT <- res[4]
    
    res <- summary(lm(Y2.RINT~G[,1]))$coef[2,]
    beta.MEAN.G.RINT <- res[1]
    p.MEAN.G.RINT <- res[4]
    
    p.mat[[2]][(nsim*(j-1))+i,1:8] <- c(beta.MEAN.G,beta.MEAN.G.RINT,beta.VAR.G,beta.VAR.G.RINT,
                           p.MEAN.G,p.MEAN.G.RINT,p.VAR.G,p.VAR.G.RINT)
    
    x <- G[,1]
    p.mat2[[2]][(nsim*(j-1))+i,1] <- as.numeric(DeviationRegressionModel(Y2,G[,1]))[4]
    p.mat2[[2]][(nsim*(j-1))+i,2] <- leveneTest(Y2~as.factor(x),center=median)$"Pr(>F)"[1]
    p.mat2[[2]][(nsim*(j-1))+i,3] <- as.numeric(bptest(Y2~x)$p.value)
    p.mat2[[2]][(nsim*(j-1))+i,4] <- coef(summary(dglm(Y2~x))$dispersion.summary)[1,4]
    
    p.mat2[[2]][(nsim*(j-1))+i,5] <- as.numeric(DeviationRegressionModel(Y2.RINT,G[,1]))[4]
    p.mat2[[2]][(nsim*(j-1))+i,6] <- leveneTest(Y2.RINT~as.factor(x),center=median)$"Pr(>F)"[1]
    p.mat2[[2]][(nsim*(j-1))+i,7] <- as.numeric(bptest(Y2.RINT~x)$p.value)
    p.mat2[[2]][(nsim*(j-1))+i,8] <- coef(summary(dglm(Y2.RINT~x))$dispersion.summary)[1,4]
    
  }
  
}
colnames(p.mat[[1]]) <- c('beta.MEAN.GxG','beta.MEAN.GxG.RINT','beta.VAR.GxG','beta.VAR.GxG.RINT',
                          'p.MEAN.GxG','p.MEAN.GxG.RINT','p.VAR.GxG','p.VAR.GxG.RINT','genetic.var.exp')
colnames(p.mat[[2]]) <- c('beta.MEAN.G','beta.MEAN.G.RINT','beta.VAR.G','beta.VAR.G.RINT',
                          'p.MEAN.G','p.MEAN.G.RINT','p.VAR.G','p.VAR.G.RINT','genetic.var.exp')

colnames(p.mat2[[1]]) <- c('DRM','LT','BP','DGLM','DRM.RINT','LT.RINT','BP.RINT','DGLM.RINT',
                           'genetic.var.exp')
colnames(p.mat2[[2]]) <- c('DRM','LT','BP','DGLM','DRM.RINT','LT.RINT','BP.RINT','DGLM.RINT',
                           'genetic.var.exp')
aggregate(p.mat2[[1]],by=list(p.mat2[[1]]$genetic.var.exp),function(x) {mean(x < 0.05,na.rm=T)})
aggregate(p.mat2[[2]],by=list(p.mat2[[2]]$genetic.var.exp),function(x) {mean(x < 0.05,na.rm=T)})

library(ggplot2)
g1 <- ggplot(p.mat[[2]],aes(x=beta.MEAN.G,y=beta.VAR.G)) + geom_point() + geom_smooth(method='lm',col='red'); g1
g2 <- ggplot(p.mat[[2]],aes(x=beta.MEAN.G.RINT,y=beta.VAR.G.RINT)) + geom_point() + geom_smooth(method='lm',col='red'); g2

# cor(p.mat[[2]]$beta.MEAN.G,p.mat[[2]]$beta.VAR.G)
# cor(p.mat[[2]]$beta.MEAN.G.RINT,p.mat[[2]]$beta.VAR.G.RINT)

g3 <- ggplot(p.mat[[1]],aes(x=beta.MEAN.GxG,y=beta.VAR.GxG)) + geom_point() + geom_smooth(method='lm',col='red'); g3
g4 <- ggplot(p.mat[[1]],aes(x=beta.MEAN.GxG.RINT,y=beta.VAR.GxG.RINT)) + geom_point() + geom_smooth(method='lm',col='red'); g4
# cor(p.mat[[1]]$beta.MEAN.GxG,p.mat[[1]]$beta.VAR.GxG)
# cor(p.mat[[1]]$beta.MEAN.GxG.RINT,p.mat[[1]]$beta.VAR.GxG.RINT)

library(cowplot)
plot_grid(g1,g2,g3,g4,ncol=2)


head(p.mat[[2]])


aggregate(Y,by=list(G[,1]),var)
aggregate(Y.RINT,by=list(G[,1]),var)
aggregate(Y.epi,by=list(G[,1]),var)
aggregate(rntransform(Y.epi),by=list(G[,1]),var)

df <- data.frame(G=G[,1],P=Y2)
ggplot(df,aes(x=P,fill=as.factor(G))) + geom_density(alpha=0.2)
df <- data.frame(G=G[,1],P=Y2.RINT)
ggplot(df,aes(x=P,fill=as.factor(G))) + geom_density(alpha=0.2)

quantile(subset(df,G==0)$P)
quantile(subset(df,G==2)$P)

summary(subset(df,G==2)$P)

as.numeric(DeviationRegressionModel(log10(Y-min(Y)+1),G[,1]))
as.numeric(DeviationRegressionModel(Y,G[,1]))

as.numeric(DeviationRegressionModel(Y,G[,1]))



x1 <- rchisq(100000,4)
x2 <- rchisq(100000,50)

quantile(x1,0.99)
quantile(x2,0.99)
median(x1)
median(x2)
