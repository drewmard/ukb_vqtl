GxE_acrossSNPs <- function(k) {
  
  vQTL <- index[k]
  print(paste0('SNP ',k,'/',length(index),': ',vQTL)) # for debugging
  
  GxE_acrossEnvirFactor <- function(j) {
    ENVIR_FACTOR <- ENVIR_NAMES[j]
    print(paste0('Environmental factor ',j,'/',length(ENVIR_NAMES),': ',ENVIR_FACTOR)) # for debugging
    
    if (ENVIR_FACTOR=='Smoking.E') {
      mod.formula <- paste(paste0('resid2'),' ~ ',
                           vQTL,'*',ENVIR_FACTOR)
    } else if (ENVIR_FACTOR=='alcohol.freq.E') { 
      mod.formula <- paste(paste0('resid3'),' ~ ',
                           vQTL,'*',ENVIR_FACTOR)
    } else if (ENVIR_FACTOR=='age') {
      mod.formula <- paste(paste0('resid4'),' ~ ',
                           vQTL,'*',ENVIR_FACTOR)
    #   # need to do
    # } else if (ENVIR_FACTOR=='sex') {
    #   pass
    } else {
      mod.formula <- paste(paste0('resid1'),' ~ ',
                           vQTL,'*',ENVIR_FACTOR)
    }
    mod.formula <- formula(mod.formula)
      
    mod1 <- lm(mod.formula,
               data=df2,na.action=na.exclude)

    mod1.coef <- coeftest(mod1, vcov = vcovHC(mod1))
    
    # }
    res1 <- mod1.coef[paste0(vQTL),]
    res2 <- mod1.coef[paste0(ENVIR_FACTOR),]
    res3 <- tryCatch({mod1.coef[paste0(vQTL,':',ENVIR_FACTOR),]},
                     error=function(e) { mod1.coef[paste0(ENVIR_FACTOR,':',vQTL),]})
    
    return(c(phenoName,vQTL,ENVIR_FACTOR,res1,res2,res3))
  }
  
  y <- lapply(1:length(ENVIR_NAMES),GxE)
  y.df <- as.data.frame(do.call(rbind,y))
  colnames(y.df) <- c('Phenotype','vQTL','E','BETA_vQTL','SE_vQTL','T_vQTL','P_vQTL',
                      'BETA_E','SE_E','T_E','P_E',
                      'BETA_GxE','SE_GxE','T_GxE','P_GxE')
  return(y.df)
}
