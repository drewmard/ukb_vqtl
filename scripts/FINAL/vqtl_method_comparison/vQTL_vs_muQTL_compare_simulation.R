library(data.table)
library(parallel)

# Function to perform DRM test for vqtls
DeviationRegressionModel <- function(PHENOTYPE,SNP) {
  X <- as.factor(SNP)
  Y.i <- tapply(PHENOTYPE, X, median)
  Z.ij <- abs(PHENOTYPE - Y.i[X])
  p <- as.numeric(summary(lm(Z.ij~SNP))$coef[2,c(1,4)])
  return(p)
}

testing <- function(j,i=1,type='gxg') {
  
  genetic_variance_explained <- genetic_variance_explained.vec[i]
  if (j %% 50 == 0) {print(paste0('Iteration ',j,', ',i,'/',length(genetic_variance_explained.vec)))}

  # generate SNPs & causal effects on phenotype
  if (type=='gxg') {
    G <- cbind(rbinom(nindiv,2,MAF1),rbinom(nindiv,2,MAF2))
    phenotype.g <- (G[,1])*(G[,2]-1)
    # phenotype.g <- (G[,1])*(G[,2])
    SNP <- G[,1]
  } else if (type=='mean') {
    SNP <- rbinom(nindiv,2,MAF1)
    phenotype.g <- SNP
  }

  # generate phenotype using environmental noise & point genetic effects
  scalar <- (1-(genetic_variance_explained))/(genetic_variance_explained);
  environmental_variance <- scalar*var(phenotype.g)
  if (phenotype_noise=='NORMAL') {
    phenotype.e <- rnorm(nindiv,0,sqrt(environmental_variance))
  } else if (phenotype_noise=='CHISQ4') {
    phenotype.e <- as.numeric(sqrt(environmental_variance)*scale(rchisq(nindiv,df=4)))
  }
  PHENOTYPE <- scale(phenotype.g + phenotype.e)[,1]
  
  # perform QTL tests
  res <- as.numeric(DeviationRegressionModel(PHENOTYPE,SNP))
  beta <- res[1]
  DRM <- res[2]
  LINEAR_REGRESSION <- as.numeric(summary(lm(PHENOTYPE~SNP))$coef[2,c(1,4)])
  
  # save results
  res <- data.frame(type=type,MAF1=MAF1,MAF2=MAF2,N=nindiv,noise=phenotype_noise,
                    h=genetic_variance_explained,
                    LR_beta=LINEAR_REGRESSION[1],LR_p=LINEAR_REGRESSION[2],
                    DRM_beta=beta,DRM_p=DRM)
  
  # library(ggplot2)
  # tmp <- data.frame(SNP1=G[,1],SNP2=G[,2],PHENO=PHENOTYPE)
  # ggplot(tmp,aes(SNP1,PHENOTYPE,col=as.factor(SNP2))) + geom_jitter(width=0.1) + geom_smooth(method='lm',se=F) + theme_bw()
  
  return(list(res,SNP,PHENOTYPE))
}

runSimulation <- function(i,type='gxg') {
  tests <- mclapply(1:nsim,testing,i=i,type=type,mc.cores = 8)
  results.tmp <- do.call(rbind,lapply(tests,function(x){x[[1]]}))
  genotypes.tmp <- do.call(cbind,lapply(tests,function(x){x[[2]]}))
  phenotypes.tmp <- do.call(cbind,lapply(tests,function(x){x[[3]]}))
  return(list(results.tmp,genotypes.tmp,phenotypes.tmp))
}


# parameters
MAF1 <- 0.4
MAF2 <- 0.4
nsim <- 1000; 
nindiv <- 10000
simulation_type='gxg'
genetic_variance_explained.vec <- seq(0.01,0.1,by=0.01);
phenotype_noise.vec <- c('NORMAL','CHISQ4')
results <- list(); genotypes <- list(); phenotypes <- list()

genetic_variance_explained.vec <- genetic_variance_explained.vec[1:3]
for (k in 1:1) {
  set.seed(03191995)
  phenotype_noise <- phenotype_noise.vec[k]
  simulation_results <- lapply(1:length(genetic_variance_explained.vec),runSimulation,type=simulation_type)
  results <- do.call(rbind,lapply(simulation_results,function(x){x[[1]]}))
  # genotypes <- do.call(cbind,lapply(simulation_results,function(x){x[[2]]}))
  # phenotypes <- do.call(cbind,lapply(simulation_results,function(x){x[[3]]}))
  results <- as.data.frame(results,stringsAsFactors = FALSE)

  #1
  results$LR_p.reject <- (results$LR_p<0.05)
  results$DRM_p.reject <- (results$DRM_p<0.05)
  res <- aggregate(results[,c('LR_p.reject','DRM_p.reject')],by = list(h=results$h),mean)
  colnames(res)[2:3] <- c('LR','DRM')
  png('~/Documents/Research/vQTL/ukb_vqtl/output/simulation/muqtl_vs_vqtl.png',width = 2700,height = 2700,res=650)
  ggplot(melt(res,id.vars='h'),aes(x=h,y=value,col=variable)) + geom_line() + geom_point() + theme_bw() + theme(panel.grid=element_blank()) + labs(col='Method',x='Variance explained by GxG',y='Power')
  dev.off()
  #2
  val <- 0.01
  results.sub <- subset(results,h==val)
  tab <- table(results.sub[,c('LR_p.reject','DRM_p.reject')])
  fisher.test(tab)$p.value
  tab[1,2]/(tab[1,1]+tab[1,2])
  tab[2,2]/(tab[2,1]+tab[2,2])
  
  #3
  val <- 0.02
  results.sub <- subset(results,h==val)
  tab <- table(results.sub[,c('LR_p.reject','DRM_p.reject')])
  fisher.test(tab)$p.value
  tab[1,2]/(tab[1,1]+tab[1,2])
  tab[2,2]/(tab[2,1]+tab[2,2])
  
  #4
  val <- 0.03
  results.sub <- subset(results,h==val)
  tab <- table(results.sub[,c('LR_p.reject','DRM_p.reject')])
  fisher.test(tab)$p.value
  tab[1,2]/(tab[1,1]+tab[1,2])
  tab[2,2]/(tab[2,1]+tab[2,2])
  
  
  dir <- '~/Documents/Research/vQTL'
  dir <- '/athena/elementolab/scratch/anm2868/vQTL'
  f <- paste0(dir,'/ukb_vqtl/output/simulation/MAF1_',MAF1,'.MAF2_',MAF2,'.NSIM_',nsim,'.NINDIV_',nindiv,'.TYPE_',simulation_type,'.NOISE_',phenotype_noise,'.txt')
  fwrite(results,f,quote = F,na = 'NA',sep = '\t',row.names = F,col.names = T)
  
  f <- paste0(dir,'/ukb_vqtl/output/simulation/MAF1_',MAF1,'.MAF2_',MAF2,'.NSIM_',nsim,'.NINDIV_',nindiv,'.TYPE_',simulation_type,'.NOISE_',phenotype_noise,'.geno.txt')
  fwrite(genotypes,f,quote = F,na = 'NA',sep = '\t',row.names = F,col.names = T)
  
  f <- paste0(dir,'/ukb_vqtl/output/simulation/MAF1_',MAF1,'.MAF2_',MAF2,'.NSIM_',nsim,'.NINDIV_',nindiv,'.TYPE_',simulation_type,'.NOISE_',phenotype_noise,'.pheno.txt')
  fwrite(phenotypes,f,quote = F,na = 'NA',sep = '\t',row.names = F,col.names = T)
  
  rm(results);rm(genotypes);rm(phenotypes)
}


# yep
phenotype_noise <- phenotype_noise.vec[k]
simulation_results <- lapply(1:length(genetic_variance_explained.vec),runSimulation,type=simulation_type)
results <- do.call(rbind,lapply(simulation_results,function(x){x[[1]]}))
genotypes <- do.call(cbind,lapply(simulation_results,function(x){x[[2]]}))
phenotypes <- do.call(cbind,lapply(simulation_results,function(x){x[[3]]}))
results <- as.data.frame(results,stringsAsFactors = FALSE)
genotypes <- as.data.frame(genotypes,stringsAsFactors=FALSE)
phenotypes <- as.data.frame(phenotypes,stringsAsFactors=FALSE)

