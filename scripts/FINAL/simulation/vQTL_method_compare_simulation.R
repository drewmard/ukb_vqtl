# initialize
library(car)
# library(lmtest)
library(dglm)
library(parallel)
library(data.table)
DeviationRegressionModel <- function(PHENOTYPE,SNP) {
  # Function to perform DRM test for vqtls
  X <- as.factor(SNP)
  Y.i <- tapply(PHENOTYPE, X, median)
  Z.ij <- abs(PHENOTYPE - Y.i[X])
  p <- as.numeric(summary(lm(Z.ij~SNP))$coef[2,c(1,4)])
  return(p)
}

testing <- function(j,i=1,type='gxg') {
  
  genetic_variance_explained <- genetic_variance_explained.vec[i]
  if (j %% 50 == 0) {print(paste0('Iteration ',j,', ',i,'/',length(genetic_variance_explained.vec)))}
  # if (j %% 1 == 0) {print(paste0('Iteration ',j,', ',i,'/',length(genetic_variance_explained.vec)))}
  
  # generate SNPs & causal effects on phenotype
  if (type=='gxg') {
    G <- cbind(rbinom(nindiv,2,MAF1),rbinom(nindiv,2,MAF2))
    phenotype.g <- G[,1]*G[,2]
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
  
  # perform vqtl tests
  res <- as.numeric(DeviationRegressionModel(PHENOTYPE,SNP))
  beta <- res[1]
  DRM <- res[2]
  LT <- leveneTest(PHENOTYPE,as.factor(SNP))$"Pr(>F)"[1]
  BT <- bartlett.test(PHENOTYPE,as.factor(SNP))$p.value
  FK <- fligner.test(PHENOTYPE,as.factor(SNP))$p.value
  DGLM <- summary(dglm(PHENOTYPE~SNP))$dispersion.summary$coef[1,4]
  TSSR <- summary(lm(PHENOTYPE^2 ~ SNP))$coef[2,4]
  
  # save results
  res <- data.frame(type=type,MAF1=MAF1,MAF2=MAF2,N=nindiv,noise=phenotype_noise,
                        h=genetic_variance_explained,
                        beta,DRM,LT,BT,FK,DGLM,TSSR)
  return(list(res,SNP,PHENOTYPE))
}

runSimulation <- function(i,type='gxg') {
  tests <- mclapply(1:nsim,testing,i=i,type=type,mc.cores = 8)
  results.tmp <- do.call(rbind,lapply(tests,function(x){x[[1]]}))
  genotypes.tmp <- do.call(cbind,lapply(tests,function(x){x[[2]]}))
  phenotypes.tmp <- do.call(cbind,lapply(tests,function(x){x[[3]]}))
  return(list(results.tmp,genotypes.tmp,phenotypes.tmp))
}

#######################

# parameters
MAF1 <- 0.4
MAF2 <- 0.4
nsim <- 1000; 
nindiv <- 10000
simulation_type='gxg'
# genetic_variance_explained.vec <- seq(0.005,0.02,by=0.005);
# genetic_variance_explained.vec <- seq(0.01,0.1,by=0.01);
genetic_variance_explained.vec <- seq(0.01,0.03,by=0.01);
phenotype_noise <- 'NORMAL' # CHISQ4
phenotype_noise.vec <- c('NORMAL','CHISQ4')
results <- list(); genotypes <- list(); phenotypes <- list()
# for (k in 1:2) {
#   phenotype_noise <- phenotype_noise.vec[k]
#   simulation_results <- lapply(1:length(genetic_variance_explained.vec),runSimulation,type=simulation_type)
#   results[[k]] <- do.call(rbind,lapply(simulation_results,function(x){x[[1]]}))
#   genotypes[[k]] <- do.call(cbind,lapply(simulation_results,function(x){x[[2]]}))
#   phenotypes[[k]] <- do.call(cbind,lapply(simulation_results,function(x){x[[3]]}))
# }
# results <- as.data.frame(do.call(rbind,results),stringsAsFactors = FALSE)
# genotypes <- as.data.frame(do.call(cbind,genotypes),stringsAsFactors=FALSE)
# phenotypes <- as.data.frame(do.call(cbind,phenotypes),stringsAsFactors=FALSE)

for (nindiv in c(30000,40000,50000)) { 

for (k in 1:1) {
  phenotype_noise <- phenotype_noise.vec[k]
  simulation_results <- lapply(1:length(genetic_variance_explained.vec),runSimulation,type=simulation_type)
  results <- do.call(rbind,lapply(simulation_results,function(x){x[[1]]}))
  genotypes <- do.call(cbind,lapply(simulation_results,function(x){x[[2]]}))
  phenotypes <- do.call(cbind,lapply(simulation_results,function(x){x[[3]]}))
  results <- as.data.frame(results,stringsAsFactors = FALSE)
  genotypes <- as.data.frame(genotypes,stringsAsFactors=FALSE)
  phenotypes <- as.data.frame(phenotypes,stringsAsFactors=FALSE)
  
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

}

# aggregate(.~MAF1+MAF2+N+noise+h+type,results[,-c(7)],mean)
# aggregate(.~MAF1+MAF2+N+noise+h+type,results[,-c(7)],function(x) {mean(x < 0.05)})


