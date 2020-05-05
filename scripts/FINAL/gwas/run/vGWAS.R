# Written by Andrew Marderstein (2020). Contact: anm2868@med.cornell.edu

# Script: run a DRM across the SNPs present in GENO. Used in combination with PLINK.

Rplink <- function(PHENO,GENO,CLUSTER,COVAR) {
  
  PHENO[which(PHENO==-9)] <- NA

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
  
  return(as.numeric(apply(GENO, 2 , deviation_regression_model)))
}


