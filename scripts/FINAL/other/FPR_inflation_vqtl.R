set.seed(031995)

N = 10000
MAF = 0.3
nsim = 10000

deviation_regression_model <-	function(SNP) 
{
  r <- tryCatch(
    {
      X <- as.factor(SNP)
      Y.i <- tapply(PHENO, X, median,na.rm=T)
      Z.ij <- abs(PHENO - Y.i[X])
      summary(lm(Z.ij~SNP))$coef[2,]
    },
    error=function(cond) {
      return(rep(NA,4))
    })
  
  return(c( length(r) , r ))
}

results_run <- function(i) {
  G1 <- rbinom(N,2,MAF)
  G2 <- rbinom(N,2,0.3)
  PHENO <- rnorm(N)
  
  DRM <- deviation_regression_model(G1)
  GxG <- summary(lm(PHENO~G1*G2))$coef
  return( c(
    DRM[2],
    DRM[5],
    GxG[1,4],
    GxG[4,4]
  ) )
}

library(parallel)
res <- mclapply(1:nsim,results_run,mc.cores = 8)
res.df <- as.data.frame(do.call(rbind,res),stringsAsFactors = F)
colnames(res.df) <- c('BETA.VAR','P.VAR','BETA.GxG','P.GxG')

res.df.sub <- subset(res.df,P.VAR < 0.05)
library(ggplot2)
g1=ggplot(res.df.sub,aes(x=BETA.VAR,y=BETA.GxG)) + geom_point() + geom_smooth(method = 'lm',se=F) + theme_bw() + theme(panel.grid=element_blank()) + labs(x=expression(beta['var']),y=expression(beta['interaction']))
g2=ggplot(res.df,aes(x=BETA.VAR,y=BETA.GxG)) + geom_point() + geom_smooth(method = 'lm',se=F) + theme_bw() + theme(panel.grid=element_blank()) + labs(x=expression(beta['var']),y=expression(beta['interaction']))
cor.test(res.df.sub$BETA.VAR,res.df.sub$BETA.GxG)
cor.test(res.df$BETA.VAR,res.df$BETA.GxG)

g3=ggplot(res.df.sub,aes(x=-log10(P.VAR),y=-log10(P.GxG))) + geom_point() + geom_smooth(method = 'lm',se=F) + theme_bw() + theme(panel.grid=element_blank()) + 
  labs(x=expression(-log[10](italic(P)['var'])),y=expression(-log[10](italic(P)['interaction'])))
g4=ggplot(res.df,aes(x=-log10(P.VAR),y=-log10(P.GxG))) + geom_point() + geom_smooth(method = 'lm',se=F) + theme_bw() + theme(panel.grid=element_blank()) + 
  labs(x=expression(-log[10](italic(P)['var'])),y=expression(-log[10](italic(P)['interaction'])))
cor.test(-log10(res.df.sub$P.VAR),-log10(res.df.sub$P.GxG))
cor.test(-log10(res.df$P.VAR),-log10(res.df$P.GxG))


library(data.table)
f <- '~/Documents/Research/vQTL/ukb_vqtl/output/other/simulated_var_gxg_fpr.txt'
res.df <- fread(f,data.table = F,stringsAsFactors = F)
fwrite(res.df,file = f,row.names = F,col.names = T,sep = '\t',quote = F,na = 'NA')

library(cowplot)
plot_grid(g2,g1,g4,g3,ncol=2)

tmp <- subset(res.df,P.VAR < 0.05); binom.test(sum(as.numeric(tmp$P.GxG<0.05)),nrow(tmp),0.05)

