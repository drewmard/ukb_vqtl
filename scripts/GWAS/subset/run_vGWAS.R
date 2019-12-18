# load
library('BEDMatrix')
library(data.table)
library(parallel)

# system arguments
args = commandArgs(trailingOnly=TRUE)
ID.df <- fread('/athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/andrew_copies/subset/ID.impute.txt',data.table = F,stringsAsFactors = F)
x1=ID.df[args[1],1]
x2=ID.df[args[1],2]
phenotype=args[2]
num_cores=args[3]


# manual input
# i=123
# ID.df <- fread('/athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/andrew_copies/subset/ID.impute.txt',data.table = F,stringsAsFactors = F)
# x1=ID.df[i,1]
# x2=ID.df[i,2]
# phenotype='lymphocyte.count.rint.ALL'
# num_cores=4

path <- paste0('/athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/andrew_copies/subset/','ukbb.',x1,'.',x2,'.impute.bed')
geno <- BEDMatrix(path)
f.pheno='/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/phenotypes_processed.80.txt'
pheno <- fread(f.pheno,data.table = F,stringsAsFactors = F); pheno <- subset(pheno,!duplicated(pheno$IID)); pheno[pheno[,phenotype]==-9,phenotype] <- NA 

geno_names <- unlist(lapply(strsplit(rownames(geno),'_'),function(x) {return(x[2])}))
ind <- which(geno_names %in% pheno$IID)
pheno <- pheno[match(geno_names[ind],pheno$IID),]

DeviationRegressionModel <- function(i) {
  if (i %% 100 == 0) {print(i)}
  # print(i)
  SNP <- geno[ind,i]
  PHENO <- pheno[,phenotype]
  X <- as.factor(SNP)
  Y.i <- tapply(PHENO, X, median,na.rm=T)
  Z.ij <- abs(PHENO - Y.i[X])
  res <- summary(lm(Z.ij~SNP))$coef[2,]
}

Fit_Model <- function(start=1,p=5000) {
  # df.save <- mclapply(1:p,function(idx) DeviationRegressionModel(as(geno$genotypes[,idx],"numeric")),mc.cores=8)
  # start=4990; p = 5000
  # nsnp <- length(colnames(geno))
  # p <- min(p,nsnp)
  df.save <- mclapply(start:p,DeviationRegressionModel,mc.cores=num_cores)
  df.save <- do.call(rbind,df.save)
  df.save <- as.data.frame(df.save)
  df.save$SNP <- colnames(geno)[start:p]
  colnames(df.save)[1:4] <- c('Estimate','Std. Error','t value','Pr(>|t|)')
  return(df.save)
}

nsnp <- length(colnames(geno))
df.results <- Fit_Model(start=1,p=nsnp)

print('Saving...')
# f.out <- paste0(prefix,'.save.txt')
f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_subset/','ukbb.',x1,'.',x2,'.',phenotype,'.txt')
print(paste0('Writing: ',f.out))
fwrite(df.results,f.out,sep='\t',quote=F,col.names = T,row.names = F,na="NA")
