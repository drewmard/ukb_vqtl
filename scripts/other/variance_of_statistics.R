library(cowplot)
library(ggplot2)
set.seed(03191995)
var.var.vec <- c()
var.mean.vec <- c()
var.var.vec2 <- c()
var.mean.vec2 <- c()
N.range <- seq(20,1000,by=20)
for (N in N.range) {
  print(N)
  nsim=10000
  mn=0
  vr=1
  x <- matrix(rnorm(N*nsim,mn,sqrt(vr)),nrow=N,ncol=nsim)
  var.var <- sd(var.vec <- apply(x,2,var))
  var.mean <- sd(mean.vec <- apply(x,2,mean))
  var.var.vec <- c(var.var.vec,var.var)
  var.mean.vec <- c(var.mean.vec,var.mean)
  
  df=15
  x <- matrix(rchisq(N*nsim,df),nrow=N,ncol=nsim)
  var.var <- sd(var.vec <- apply(x,2,var))
  var.mean <- sd(mean.vec <- apply(x,2,mean))
  var.var.vec2 <- c(var.var.vec2,var.var)
  var.mean.vec2 <- c(var.mean.vec2,var.mean)
}

var.var.vec4 <- c()
var.mean.vec4 <- c()
var.var.vec5 <- c()
var.mean.vec5 <- c()
df.range <- seq(1,20,by=1)
N=100
for (df in df.range) {
  
  mn=mn
  vr=1
  x <- matrix(rnorm(N*nsim,mn,1),nrow=N,ncol=nsim)
  var.var <- sd(var.vec <- apply(x,2,var))
  var.mean <- sd(mean.vec <- apply(x,2,mean))
  var.var.vec5 <- c(var.var.vec5,var.var)
  var.mean.vec5 <- c(var.mean.vec5,var.mean)
  
  
  x <- matrix(rchisq(N*nsim,df),nrow=N,ncol=nsim)
  var.var <- sd(var.vec <- apply(x,2,var))
  var.mean <- sd(mean.vec <- apply(x,2,mean))
  var.var.vec4 <- c(var.var.vec4,var.var)
  var.mean.vec4 <- c(var.mean.vec4,var.mean)

}

var.var.vec3 <- c()
var.mean.vec3 <- c()
df.range2.orig <- (c(0.01,0.1,0.25,0.5,0.75,1,4,9,20))
df.range2 <- sqrt(df.range2.orig)
(df.range2)
for (df in df.range2) {
  print(df)
  nsim=10000
  mn=0
  vr=df
  x <- matrix(rnorm(N*nsim,mn,sqrt(vr)),nrow=N,ncol=nsim)
  var.var <- sd(var.vec <- apply(x,2,var))
  var.mean <- sd(mean.vec <- apply(x,2,mean))
  var.var.vec3 <- c(var.var.vec3,var.var)
  var.mean.vec3 <- c(var.mean.vec3,var.mean)
  
}

df1 <- data.frame(Dist='Normal',
                  N=N.range,
                  Param=1,
                  Var.Mean=var.mean.vec,
                  Var.Var=var.var.vec,
                  Ratio=var.var.vec/var.mean.vec)
df2 <- data.frame(Dist='Chi-Square',
                  N=N.range,
                  Param=15,
                  Var.Mean=var.mean.vec2,
                  Var.Var=var.var.vec2,
                  Ratio=var.var.vec2/var.mean.vec2)
df3 <- data.frame(Dist='Normal',
           N=N,
           Param=df.range2.orig,
           Var.Mean=var.mean.vec3,
           Var.Var=var.var.vec3,
           Ratio=var.var.vec3/var.mean.vec3)
df4 <- data.frame(Dist='Chi-Square',
                  N=N,
                  Param=df.range,
                  Var.Mean=var.mean.vec4,
                  Var.Var=var.var.vec4,
                  Ratio=var.var.vec4/var.mean.vec4)
df5 <- data.frame(Dist='Normal',
                  N=N,
                  Param=df.range,
                  Var.Mean=var.mean.vec5,
                  Var.Var=var.var.vec5,
                  Ratio=var.var.vec5/var.mean.vec5)

#############

g1.1 <- ggplot(df1,aes(x=N)) +
  geom_line(aes(y=Var.Mean),col='black') +
  geom_line(aes(y=Var.Var),col='red') +
  theme_bw() + theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5)) +
  labs(x='Sample size',y='Standard error',title=df1$Dist[1])

g1.2 <- ggplot(df1,aes(x=N.range,y=Ratio)) +
  geom_line(col='blue') +
  theme_bw() + theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5)) +
  labs(x='Sample Size',y='Ratio',title=df1$Dist[1])

g2.1 <- ggplot(df2,aes(x=N)) +
  geom_line(aes(y=Var.Mean),col='black') +
  geom_line(aes(y=Var.Var),col='red') +
  theme_bw() + theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5)) +
  labs(x='Sample size',y='Standard error',title=df2$Dist[1])

g2.2 <- ggplot(df2,aes(x=N.range,y=Ratio)) +
  geom_line(col='blue') +
  theme_bw() + theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5)) +
  labs(x='Sample Size',y='Ratio',title=df2$Dist[1])

g3.1 <- ggplot(df3,aes(x=Param)) +
  geom_line(aes(y=Var.Mean),col='black') +
  geom_line(aes(y=Var.Var),col='red') +
  theme_bw() + theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5)) +
  labs(x='Normal distribution variance',y='Standard error',title=df3$Dist[1])

g3.2 <- ggplot(df3,aes(x=Param,y=Ratio)) +
  geom_line(col='blue') +
  theme_bw() + theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5)) +
  labs(x='Normal distribution variance',y='Ratio',title=df3$Dist[1])

g4.1 <- ggplot(df4,aes(x=Param)) +
  geom_line(aes(y=Var.Mean),col='black') +
  geom_line(aes(y=Var.Var),col='red') +
  theme_bw() + theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5)) +
  labs(x='Degrees of freedom',y='Standard error',title=df4$Dist[1])

g4.2 <- ggplot(df4,aes(x=Param,y=Ratio)) +
  geom_line(col='blue') +
  theme_bw() + theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5)) +
  labs(x='Degrees of freedom',y='Ratio',title=df4$Dist[1])

g5.1 <- ggplot(df5,aes(x=Param)) +
  geom_line(aes(y=Var.Mean),col='black') +
  geom_line(aes(y=Var.Var),col='red') +
  theme_bw() + theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5)) +
  labs(x='Normal distribution mean',y='Standard error',title=df3$Dist[1])

g5.2 <- ggplot(df5,aes(x=Param,y=Ratio)) +
  geom_line(col='blue') +
  theme_bw() + theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5)) +
  labs(x='Normal distribution mean',y='Ratio',title=df3$Dist[1])

x <- 2
png('~/Documents/Research/vQTL/ukb_vqtl/output/variance_of_parameter/variance_of_parameter_panel.png',width = 2000*x,height = 3800*x,res=350*x)
plot_grid(g1.1,g1.2,g2.1,g2.2,g5.1,g5.2,g3.1,g3.2,g4.1,g4.2,nrow = 5)
dev.off()

N <- 10
vr=0.01
x <- matrix(rnorm(N*nsim,mn,sqrt(vr)),nrow=N,ncol=nsim)
var.var <- sd(var.vec <- apply(x,2,var))
var.mean <- sd(mean.vec <- apply(x,2,mean))
var.var/var.mean
