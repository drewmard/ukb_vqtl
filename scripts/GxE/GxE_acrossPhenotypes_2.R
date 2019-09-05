GxE_acrossPhenotypes <- function(i) {
  phenoName <- PHENOTYPE_NAMES[i]
  print(paste0('Phenotype ',i,'/',length(PHENOTYPE_NAMES),': ',phenoName )) # for debugging
  
  GxE <- function(j) {
    ENVIR_FACTOR <- ENVIR_NAMES[j]
    print(paste0('Environmental factor ',j,'/',length(ENVIR_NAMES),': ',ENVIR_FACTOR)) # for debugging
    
      
      mod.formula <- paste(paste0('resid'),' ~ ',
                           vQTL,'*',ENVIR_FACTOR)
      
      mod1 <- lm(mod.formula,
                 data=df2,na.action=na.exclude)
      mod1
    },error=function(e) {
      mod.formula <- paste(paste0(phenoName),' ~ age+age2+sex+
                           PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
                           PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+',
                           vQTL,'*',ENVIR_FACTOR)
      
      # time.since.period2+time.since.period2.dummy+menopause2+ # women only
      # bmi2.dummy+bmi2+bmi2*age+',
      if (ENVIR_FACTOR!='Smoking.E') {
        mod.formula <- paste0(mod.formula,'+Smoking+Smoking.dummy')
      } 
      if (ENVIR_FACTOR!='alcohol.freq.E') { 
        mod.formula <- paste0(mod.formula,'+alcohol.freq2+alcohol.freq2.dummy')
      }
      mod.formula <- formula(mod.formula)
      
      mod1 <- lm(mod.formula,
                 data=df2,na.action=na.exclude)
      mod1
    })
    
    # mod1.coef <- summary(mod1)$coef
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
