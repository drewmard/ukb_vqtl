Rplink <- function(PHENO,GENO,CLUSTER,COVAR) {
  
  PHENO[which(PHENO==-9)] <- NA

  association_test <- function(SNP)
  {
    r <- tryCatch(
      {
      summary(lm(PHENO~SNP))$coef[2,]
      },
      error=function(cond) {
        return(rep(NA,4))
      })
    return(c( length(r) , r ))
  }
  deviation_regression_model <-	function(SNP) 
  {
    r <- tryCatch(
      {
      X <- as.factor(SNP)
      Y.i <- tapply(PHENO, X, median,na.rm=T)
#      Y.i <- tapply(PHENO, X, function(x) {return(as.numeric(quantile(x,probs=0.9,na.rm=T)))})
      Z.ij <- abs(PHENO - Y.i[X])
      summary(lm(Z.ij~SNP))$coef[2,]
      },
      error=function(cond) {
        return(rep(NA,4))
      })

    return(c( length(r) , r ))
  }
#  return(as.numeric(apply(GENO, 2 , association_test)))
  return(as.numeric(apply(GENO, 2 , deviation_regression_model)))
}


