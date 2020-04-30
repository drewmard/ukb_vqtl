N=1000
G <- rbinom(N,2,0.4)
Y <- rep(NA,N)
Y[G==0] <- rnorm(sum(G==0),0,1)
Y[G==1] <- rnorm(sum(G==1),0.5,2)
Y[G==2] <- rnorm(sum(G==2),1,3)
G2 <- G
G2[G==0] <- 'AA'
G2[G==1] <- 'AT'
G2[G==2] <- 'TT'

library(ggplot2)
g <- ggplot(data.frame(G=G2,Y), aes(x=as.factor(G),y=Y)) + geom_jitter(width=0.1,col='orange',alpha=0.4) + geom_boxplot(outlier.shape = NA) + theme_bw() + 
  theme(panel.grid = element_blank()) + labs(y='Phenotype') + scale_y_continuous(breaks=c(min(Y),max(Y)),labels=c('Low','High')) + theme(axis.title = element_blank())
g
x <- 1.5
png(filename='~/Documents/Research/vQTL/vQTL.png',height = 4.95/x,width=7.66/x,units = 'in',res = 500)
g
dev.off()




N=1000
G <- rbinom(N,2,0.4)
Y <- rep(NA,N)
Y[G==0] <- rnorm(sum(G==0),0,1)
Y[G==1] <- rnorm(sum(G==1),0,1)
Y[G==2] <- rnorm(sum(G==2),0,1)
G2 <- G
G2[G==0] <- 'AA'
G2[G==1] <- 'AT'
G2[G==2] <- 'TT'
G2.noeff <- G2
Y.noeff <- Y

library(ggplot2)
g <- ggplot(data.frame(G=G2,Y), aes(x=as.factor(G),y=Y)) + geom_jitter(width=0.1,col='orange',alpha=0.4) + geom_boxplot(outlier.shape = NA) + theme_bw() + 
  theme(panel.grid = element_blank()) + labs(y='Phenotype') + scale_y_continuous(breaks=c(min(Y),max(Y)),labels=c('Low','High')) + theme(axis.title = element_blank())
g

N=1000
G <- rbinom(N,2,0.4)
Y <- rep(NA,N)
Y[G==0] <- rnorm(sum(G==0),0,1)
Y[G==1] <- rnorm(sum(G==1),1,1)
Y[G==2] <- rnorm(sum(G==2),2,1)
G2 <- G
G2[G==0] <- 'AA'
G2[G==1] <- 'AT'
G2[G==2] <- 'TT'

library(ggplot2)
g <- ggplot(data.frame(G=G2,Y), aes(x=as.factor(G),y=Y)) + geom_jitter(width=0.1,col='orange',alpha=0.4) + geom_boxplot(outlier.shape = NA) + theme_bw() + 
  theme(panel.grid = element_blank()) + labs(y='Phenotype') + scale_y_continuous(breaks=c(min(Y),max(Y)),labels=c('Low','High')) + theme(axis.title = element_blank())
g

df=data.frame(G=c(G2.noeff,G2),Y=c(Y.noeff,Y),E=rep(c(0,1),each=N))
library(ggplot2)
g <- ggplot(df, aes(x=as.factor(G),y=Y,col=as.factor(E),fill=as.factor(E))) + geom_jitter(width=0.1,alpha=0.4,aes(group=E)) + geom_boxplot(outlier.shape = NA) + theme_bw() + 
  theme(panel.grid = element_blank()) + labs(y='Phenotype') + scale_y_continuous(breaks=c(min(Y),max(Y)),labels=c('Low','High')) + theme(axis.title = element_blank())
g




