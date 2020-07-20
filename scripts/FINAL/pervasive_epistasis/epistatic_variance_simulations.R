library(parallel)
library(ggplot2)
library(cowplot)
library(stringr)
library(data.table)
set.seed(03191995)

# Initialize
Nindiv=10000
MAF=0.25
Mcausalsnp=400; Mtotsnp=1000
Ninteraction <- Mtotsnp * (Mtotsnp-1) / 2

print('Generating genotypes...')
G <- matrix(rbinom(Nindiv*Mtotsnp,2,MAF),nrow = Nindiv,ncol = Mtotsnp)

print('Generating interactions file...')
tmp <- model.matrix( ~.^2, data=as.data.frame(G[,1:Mcausalsnp]))[,-1]

print('Generating phenotypes...')
X <- tmp[,(Mcausalsnp+1):ncol(tmp)]
beta <- rnorm(ncol(X))
# beta <- sort(beta,decreasing = FALSE) # sort1
beta <- sort(beta,decreasing = TRUE) # sort2

# pheno.additive <- as.numeric(tmp[,1:Mcausalsnp] %*% beta[1:Mcausalsnp])
# pheno.interaction <- as.numeric(tmp[,(Mcausalsnp+1):ncol(tmp)] %*% beta[(Mcausalsnp+1):ncol(tmp)])

pheno.interaction <- as.numeric(X %*% beta)

# pheno <- scale(pheno.additive) + scale(pheno.interaction)
pheno <- scale(pheno.interaction)
# pheno <- pheno.interaction
pheno <- as.numeric(pheno)

# truth dataframe
# tmp.col_names <- colnames(tmp)[(Mcausalsnp+1):ncol(tmp)]
# i.j <- as.numeric(unlist(strsplit(str_replace_all(tmp.col_names,'V',''),':')))
# i <- i.j[seq(1,length(i.j),by=2)]
# j <- i.j[seq(2,length(i.j),by=2)]
# df.truth <- data.frame(i,j,true_beta=beta[(Mcausalsnp+1):ncol(tmp)])

X.col_names <- colnames(X)
i.j <- as.numeric(unlist(strsplit(str_replace_all(X.col_names,'V',''),':')))
i <- i.j[seq(1,length(i.j),by=2)]
j <- i.j[seq(2,length(i.j),by=2)]
df.truth <- data.frame(i,j,true_beta=beta)

genetic_association <- function(i) {
  mod <- lm(pheno~G[,i])
  res <- summary(mod)$coef[2,c(1,4)]
  return(res)
}
deviation_regression_model <-	function(i) {
  SNP <- G[,i]
  X <- as.factor(SNP)
  Y.i <- tapply(pheno, X, median,na.rm=T)
  Z.ij <- abs(pheno - Y.i[X])
  res <- summary(lm(Z.ij~SNP))$coef[2,c(1,4)]
  return(res)
}
gxg_association <- function(i,j) {
  mod <- lm(pheno~G[,i]*G[,j])
  p <- summary(mod)$coef[4,c(1,4)]
  p <- c(i,j,p)
  return(p)
}
screen_gxg <- function(i) {
  print(i)
  res <- lapply((i+1):Mtotsnp,function(j) {gxg_association(i,j)})
  res <- do.call(rbind,res)
  return(res)
}
qq_ggplot <- function(pvector) {
  o = -log10(sort(pvector, decreasing = FALSE))
  e = -log10(ppoints(length(pvector)))
  dataf <- data.frame(e,o)
  g <- ggplot(dataf,aes(x=e,y=o)) + geom_point() + geom_abline(slope=1,intercept=0,col='red') + 
    labs(x=expression(Expected ~ ~-log[10](italic(p))),y=expression(Observed ~ ~-log[10](italic(p)))) +
    theme_bw() + theme(panel.grid = element_blank())
  return(g)
}
# generate_ij <- function(max_n) {
#   for (i in 1:(max_n-1)) {
#     tmp <- data.frame(i=i,j=seq(i+1,max_n))
#     if (i==1) {
#       tmp.sv <- tmp
#     } else {
#       tmp.sv <- rbind(tmp.sv,tmp)
#     }
#   }
#   return(tmp.sv)
# }

print('muGWAS...')
res <- mclapply(1:Mtotsnp,genetic_association,mc.cores = 8)
res.add <- do.call(rbind,res)
colnames(res.add) <- c('BETA.MEAN','P.MEAN')
# g1 <- qq_ggplot(res.add);#g1

print('vGWAS...')
res <- mclapply(1:Mtotsnp,deviation_regression_model,mc.cores = 8)
res.var <- do.call(rbind,res)
colnames(res.var) <- c('BETA.VAR','P.VAR')

res.mg <- as.data.frame(cbind(res.add,res.var),stringsAsFactors = F)
res.mg$causal <- 0; res.mg$causal[1:Mcausalsnp] <- 1

# print('Mean-var plots...')
# gA <- ggplot(res.mg,aes(x=BETA.MEAN,y=BETA.VAR)) + geom_point(aes(col=as.factor(causal))) + geom_smooth(method='lm',col='red') + theme_bw() + theme(panel.grid = element_blank(),legend.position = 'none')
# gB <- ggplot(res.mg,aes(x=-log10(P.MEAN),y=-log10(P.VAR))) + geom_point(aes(col=as.factor(causal))) + geom_smooth(method='lm',col='red') + theme_bw() + theme(panel.grid = element_blank(),legend.position = 'none')
# x=5
# f <- paste0('~/Documents/Research/vQTL/ukb_vqtl/output/sim/epistatic_variance_simulations_dataset.N_',Nindiv,'.MAF_',MAF,'.Mcausal_',Mcausalsnp,'.Mtot_',Mtotsnp,'.mean_var.png')
# png(f,width = 3000*x,height=1000*x,res = 200*x)
# plot_grid(gA,gB,ncol = 2)
# dev.off()

# f <- paste0('~/Documents/Research/vQTL/ukb_vqtl/output/sim/epistatic_variance_simulations_dataset.N_',Nindiv,'.MAF_',MAF,'.Mcausal_',Mcausalsnp,'.Mtot_',Mtotsnp,'.mean_var.txt')
f <- paste0('~/Documents/Research/vQTL/ukb_vqtl/output/sim/epistatic_variance_simulations_dataset.N_',Nindiv,'.MAF_',MAF,'.Mcausal_',Mcausalsnp,'.Mtot_',Mtotsnp,'.mean_var.sort2.txt')
fwrite(res.mg,f,sep = '\t',quote = F,na='NA',col.names = T,row.names = F)

print('GxG testing...')
res <- mclapply(1:(Mtotsnp-1),screen_gxg,mc.cores = 8)
print('GxG testing completed.')
gxg.df <- as.data.frame(do.call(rbind,res),stringsAsFactors = F)
colnames(gxg.df) <- c('i','j','BETA.GxG','P.GxG')

tmp2 <- res.mg; colnames(tmp2) <- paste0(colnames(tmp2),'.j'); tmp2$j <- 1:Mtotsnp
gwas <- merge(gxg.df,tmp2,by='j')
tmp2 <- res.mg; colnames(tmp2) <- paste0(colnames(tmp2),'.i'); tmp2$i <- 1:Mtotsnp
gwas <- merge(gwas,tmp2,by='i')
gwas <- merge(gwas,df.truth,by=c('i','j'),all = TRUE)
gwas$true_beta[is.na(gwas$true_beta)] <- 0

# f <- paste0('~/Documents/Research/vQTL/ukb_vqtl/output/sim/epistatic_variance_simulations_dataset.N_',Nindiv,'.MAF_',MAF,'.Mcausal_',Mcausalsnp,'.Mtot_',Mtotsnp,'.txt')
f <- paste0('~/Documents/Research/vQTL/ukb_vqtl/output/sim/epistatic_variance_simulations_dataset.N_',Nindiv,'.MAF_',MAF,'.Mcausal_',Mcausalsnp,'.Mtot_',Mtotsnp,'.sort2.txt')
fwrite(gwas,f,sep = '\t',quote = F,na='NA',col.names = T,row.names = F)

g1 <- ggplot(data.frame(i=1:nrow(res.mg),P=-log10(res.mg$P.MEAN),causal=c(rep(1,Mcausalsnp),rep(0,Mtotsnp-Mcausalsnp))),aes(i,P,col=causal)) + geom_point() + labs(x='index',y='-log10 p-value (marginal mean)') + theme_bw() +  theme(legend.position = 'none')
g2 <- ggplot(data.frame(i=1:nrow(res.mg),BETA=(res.mg$BETA.MEAN),causal=c(rep(1,Mcausalsnp),rep(0,Mtotsnp-Mcausalsnp))),aes(i,BETA,col=causal)) + geom_point() + labs(x='index',y='marginal mean effect') + theme_bw() +  theme(legend.position = 'none')
plot_grid(g1,g2,ncol=2)
# f <- paste0('~/Documents/Research/vQTL/ukb_vqtl/output/sim/epistatic_variance_simulations_dataset.N_',Nindiv,'.MAF_',MAF,'.Mcausal_',Mcausalsnp,'.Mtot_',Mtotsnp,'.est_vs_true_int_eff.png')
# png(f,width = 2000*x,height=2000*x,res = 300*x)
# gC <- ggplot(gwas,aes(x=BETA.GxG,y=true_beta)) + geom_point(aes(col=as.factor(causal.i==1 & causal.j==1))) + theme_bw() + theme(legend.position = 'none',panel.grid = element_blank()) + labs(x='Estimated interaction effects',y='True interaction effects')
# gC
# dev.off()

# g1 <- qq_ggplot(res.mg$P.MEAN);#g1
# g2 <- qq_ggplot(gwas$P.GxG);#g1
# gwas.sub <- subset(gwas,P.MEAN.i < 5e-8 & P.MEAN.j < 5e-8)
# g3 <- qq_ggplot(gwas.sub$P.GxG);#g2

# x=5
# f <- paste0('~/Documents/Research/vQTL/ukb_vqtl/output/sim/epistatic_variance_simulations_dataset.N_',Nindiv,'.MAF_',MAF,'.Mcausal_',Mcausalsnp,'.Mtot_',Mtotsnp,'.QQgrid.png')
# png(f,width = 3000*x,height=1000*x,res = 300*x)
# plot_grid(g1,g2,g3,ncol = 3)
# dev.off()


