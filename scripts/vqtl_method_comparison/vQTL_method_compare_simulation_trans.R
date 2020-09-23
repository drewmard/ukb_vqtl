# initialize
library(car)
library(dglm)
library(parallel)
library(data.table)
library(quantreg)
# library(gJLS)

source('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/bin/rntransform.R')

DeviationRegressionModel <- function(PHENOTYPE,SNP) {
  # Function to perform DRM test for vqtls
  X <- as.factor(SNP)
  Y.i <- tapply(PHENOTYPE, X, median)
  Z.ij <- abs(PHENOTYPE - Y.i[X])
  p <- as.numeric(summary(lm(Z.ij~SNP))$coef[2,c(1,4)])
  return(p)
}

run_SVLM <- function(PHENOTYPE,SNP) {
  resid <- residuals(summary(lm(PHENOTYPE ~ SNP))); SVLM <- summary(lm((resid^2)~SNP))$coef[2,4]
  return(SVLM)
}


gS_test <- function (model, data, correlation = NULL) {
  if (sum(is.na(data)) > 0) 
    stop("missing value(s) not allowed")
  model <- terms(model)
  lm1 <- rq(model, tau = 0.5, data = data)
  data$d1 <- abs(resid(lm1))
  model2 = reformulate(attr(model, "term.labels"), response = "d1")
  if (is.null(correlation)) {
    model2 = terms(model2)
    model0 = model2[-(1:2)]
    fit0 <- lm(model0, data = data)
    fit <- lm(model2, data = data)
    aovfit <- anova(fit0, fit)
    gS_F <- aovfit[2, 5]
    numDF <- aovfit[2, 3]
    denDF <- aovfit[2, 1]
    gS_p <- aovfit[2, 6]
  }
  else {
    fit <- gls(model2, data = data, correlation = correlation, 
               method = "ML", control = lmeControl(opt = "optim"))
    aovfit <- anova(fit, Terms = 2:dim(anova(fit))[1])
    gS_F <- aovfit[1, 2]
    numDF <- aovfit[1, 1]
    denDF <- fit$dims$N - fit$dims$p
    gS_p <- aovfit[1, 3]
  }
  return(cbind(gS_F, numDF, denDF, gS_p))
}

gL_test <- function (model, data, correlation = NULL) {
  if (sum(is.na(data)) > 0) 
    stop("missing value(s) not allowed")
  if (is.null(correlation)) {
    model = terms(model)
    model0 = model[-(1:2)]
    fit0 <- lm(model0, data = data)
    fit <- lm(model, data = data)
    aovfit <- anova(fit0, fit)
    gL_F <- aovfit[2, 5]
    numDF <- aovfit[2, 3]
    denDF <- aovfit[2, 1]
    gL_p <- aovfit[2, 6]
  }
  else {
    fit <- gls(model, data = data, correlation = correlation, 
               method = "ML", control = lmeControl(opt = "optim"))
    aovfit <- anova(fit, Terms = 2:dim(anova(fit))[1])
    gL_F <- aovfit[1, 2]
    numDF <- aovfit[1, 1]
    denDF <- fit$dims$N - fit$dims$p
    gL_p <- aovfit[1, 3]
  }
  return(cbind(gL_F, numDF, denDF, gL_p))
}

gJLS_test <- function (model.loc, model.scale, data, correlation = NULL) {
  if (sum(is.na(data)) > 0) 
    stop("missing value(s) not allowed")
  gL <- gL_test(model = model.loc, data = data, correlation = correlation)
  gS <- gS_test(model = model.scale, data = data, correlation = correlation)
  t_JLS <- -2 * log(gL[4]) - 2 * log(gS[4])
  p_JLS <- 1 - pchisq(t_JLS, 4)
  dataf = data.frame(rbind(gL, gS, c(t_JLS, 4, "NA", p_JLS)))
  names(dataf) = c("Statistic", "df1", "df2", "p-value")
  rownames(dataf) = c("gL", "gS", "gJLS")
  return(dataf)
}

testing <- function(j,i=1,type='gxg') {
  
  genetic_variance_explained <- genetic_variance_explained.vec[i]
  if (j %% 50 == 0) {print(paste0('Iteration ',j,', ',i,'/',length(genetic_variance_explained.vec)))}

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
  
  # transform phenotypes
  min_phenotype <- abs(min(PHENOTYPE))+1
  tested_phenotype <- list()
  tested_phenotype[['orig']]  <- PHENOTYPE
  tested_phenotype[['rint']]  <- rntransform(PHENOTYPE)
  tested_phenotype[['log']] <- log(PHENOTYPE+min_phenotype)
  tested_phenotype[['recip']]  <- 1/(PHENOTYPE+min_phenotype)
  
  # perform vqtl tests
  res <- list()
  for (transformation in c('orig','rint','log','recip')) {
  # for (transformation in c('orig')) {
    print(transformation)
    PHENOTYPE <- tested_phenotype[[transformation]]
    DRM.time <- system.time(DRM_res <- as.numeric(DeviationRegressionModel(PHENOTYPE,SNP)))['elapsed']
    beta <- DRM_res[1]
    DRM <- DRM_res[2]
    LT.time <- system.time(LT <- leveneTest(PHENOTYPE,as.factor(SNP),center= 'mean')$"Pr(>F)"[1])['elapsed']
    BF.time <- system.time(BF <- leveneTest(PHENOTYPE,as.factor(SNP),center = 'median')$"Pr(>F)"[1])['elapsed']
    BT.time <- system.time(BT <- bartlett.test(PHENOTYPE,as.factor(SNP))$p.value)['elapsed']
    FK.time <- system.time(FK <- fligner.test(PHENOTYPE,as.factor(SNP))$p.value)['elapsed']
    DGLM.time <- system.time(DGLM <- summary(dglm(PHENOTYPE~SNP,~SNP))$dispersion.summary$coef[2,4])['elapsed']
    TSSR.time <- system.time(TSSR <- summary(lm(PHENOTYPE^2 ~ SNP))$coef[2,4])['elapsed']
    SVLM.time <- system.time(SVLM <- run_SVLM(PHENOTYPE,SNP))['elapsed']
    gJLS.time <- system.time(gJLS_res <- gJLS_test(lm(PHENOTYPE ~ SNP),lm(PHENOTYPE ~ SNP),data=data.frame(SNP,PHENOTYPE)))['elapsed']
    gS <- (gJLS_res['gS','p-value'])
    gJLS <- (gJLS_res['gJLS','p-value'])
    # save results
    res[[transformation]] <- data.frame(type=type,transformation=transformation,
                                        MAF1=MAF1,MAF2=MAF2,N=nindiv,noise=phenotype_noise,
                                        h=genetic_variance_explained,
                                        beta,DRM,LT,BF,BT,FK,DGLM,TSSR,SVLM,gS,gJLS,
                                        DRM.time,LT.time,BF.time,BT.time,FK.time,DGLM.time,TSSR.time,SVLM.time,gJLS.time)
  }
  res <- do.call(rbind,res)
  return(res)
}

runSimulation <- function(i,type='gxg') {
  tests <- mclapply(1:nsim,testing,i=i,type=type,mc.cores = 16)
  results.tmp <- do.call(rbind,tests)
  # genotypes.tmp <- do.call(cbind,lapply(tests,function(x){x[[2]]}))
  # phenotypes.tmp <- do.call(cbind,lapply(tests,function(x){x[[3]]}))
  return(results.tmp)
}

#######################

# parameters
MAF1 <- 0.4
MAF2 <- 0.4
nsim <- 1000; 
# nindiv <- 1000
nindiv <- 250000
# genetic_variance_explained.vec <- seq(0.002,0.02,by=0.002)
genetic_variance_explained.vec <- seq(0.002,0.02,by=0.002)[3:5]
# phenotype_noise <- 'NORMAL'
# phenotype_noise <- 'CHISQ4'
# simulation_type='gxg'
for (phenotype_noise in c('CHISQ4','NORMAL')) {
  for (simulation_type in c('gxg','mean')) {
    simulation_results <- lapply(1:length(genetic_variance_explained.vec),runSimulation,type=simulation_type)
    results <- do.call(rbind,simulation_results)
    results$gS <- as.numeric(as.character(results$gS))
    results$gJLS <- as.numeric(as.character(results$gJLS))
    # dir <- '~/Documents/Research/vQTL'
    dir <- '/athena/elementolab/scratch/anm2868/vQTL'
    f <- paste0(dir,'/ukb_vqtl/output/simulation/MAF1_',MAF1,'.MAF2_',MAF2,'.NSIM_',nsim,'.NINDIV_',nindiv,'.TYPE_',simulation_type,'.NOISE_',phenotype_noise,'.txt')
    fwrite(results,f,quote = F,na = 'NA',sep = '\t',row.names = F,col.names = T)
  }
}
#############

# aggregate(results[,c('DRM','LT','BF','BT','FK','DGLM','TSSR','SVLM','gJLS')],by=list(results$h,results$transformation),function(x) mean(x<0.05))
# aggregate(results[,c('DRM.time','LT.time','BF.time','BT.time','FK.time','DGLM.time','TSSR.time','SVLM.time','gJLS.time')],by=list(results$h),mean)



