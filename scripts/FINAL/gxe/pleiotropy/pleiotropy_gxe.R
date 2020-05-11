library(data.table)
library(BEDMatrix)
# read in CAD
f<-'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/CAD.txt'
CAD <- fread(f,data.table = F,stringsAsFactors = F)
f<-'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/Diabetes.txt'
Diabetes <- fread(f,data.table = F,stringsAsFactors = F)
f<-'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/HBP.txt'
HBP <- fread(f,data.table = F,stringsAsFactors = F)

# align to covariate file
for (s in c('80','20')) {
  pheno='bmi'; 
  # s='80'; 
  f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/full_data.',pheno,'.GxE.',s,'.txt')
  full_dataset <- fread(f,data.table = F,stringsAsFactors = F)
  f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/diet_data.',pheno,'.GxE.',s,'.txt')
  diet_dataset <- fread(f,data.table = F,stringsAsFactors = F)
  full_dataset$DIET_SCORE <- diet_dataset[match(full_dataset$IID,diet_dataset$IID),'DIET_SCORE']
  full_dataset$CAD <- CAD[match(full_dataset$IID,CAD$eid),'CAD']
  full_dataset$CAD[is.na(full_dataset$bmi.na)] <- NA
  full_dataset$Diabetes <- Diabetes[match(full_dataset$IID,Diabetes$eid),'Diabetes']
  full_dataset$Diabetes[is.na(full_dataset$bmi.na)] <- NA
  full_dataset$HBP <- HBP[match(full_dataset$IID,HBP$eid),'HBP']
  full_dataset$HBP[is.na(full_dataset$bmi.na)] <- NA

  # read in geno information
  f.geno <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.',pheno,'.merged_subset2')
  geno <- BEDMatrix(f.geno)
  geno_names <- unlist(lapply(strsplit(rownames(geno),'_'),function(x) {return(x[2])}))
  SNP_names <- colnames(geno); g <- regexpr("_[^_]*$", SNP_names)-1; SNP_names <- substring(SNP_names,1,g)
  ind <- which(geno_names %in% full_dataset$IID)
  
  # read gxe results
  library(data.table)
  library(BEDMatrix)
  library(parallel)
  
  # identify all significant interactions of interest
  # try P < 0.05 to start. the GxE replicate well
  library(data.table)
  pheno <- 'bmi'
  
  #################################################################
  # GxE
  
  s2='20';results.20 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.',s2,'.diet_score.more_snp.txt'),data.table = F,stringsAsFactors = F)
  s2='80';results.80 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.',s2,'.diet_score.more_snp.txt'),data.table = F,stringsAsFactors = F)
  results.mg.diet <- merge(results.80,results.20,by=c('SNP','E'))
  
  s2='20';results.20 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.',s2,'.ext.more_snp.txt'),data.table = F,stringsAsFactors = F)
  s2='80';results.80 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.',s2,'.ext.more_snp.txt'),data.table = F,stringsAsFactors = F)
  results.mg.all <- merge(results.80,results.20,by=c('SNP','E'))
  results.mg <- rbind(results.mg.diet,results.mg.all)
  g <- regexpr("_[^_]*$", results.mg$SNP)-1
  results.mg$SNP <- substring(results.mg$SNP,1,g)
  results <- results.mg
  # results <- subset(results.mg,E!='DIET_SCORE')
  results$FDR <- p.adjust(results[,4],method = 'fdr')
  results <- subset(results,(FDR<0.1) & (sign(results[,3])==sign(results[,5])))
  
  # iterator through environments
  GxE <- function(i) {
    print(i)
    SNP <- results$SNP[i]
    envir_name <- results$E[i]
    index <- which(SNP_names == SNP)
    full_dataset$SNP <- geno[ind,index]
    mod.formula <- formula(paste(phenoName,' ~ age+age2+genotyping.array+sex+age*sex+age2*sex+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
               PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+',
                               envir_name,'*SNP'))
    mod <- glm(mod.formula,data=full_dataset,family=binomial(link="logit"))
    res <- summary(mod)$coef[paste0(envir_name,':SNP'),c("Estimate","Pr(>|z|)")]
    return(res)
  }
  
  Fit_Model_Int <- function(start=1,p=5000) {
    df.save <- mclapply(start:p,GxE,mc.cores=8)
    df.save <- do.call(rbind,df.save)
    df.save <- as.data.frame(df.save)
    return(df.save)
  }
  
  phenoName.vec <- c('Diabetes','HBP','CAD')
  for (i in 1:length(phenoName.vec)) {
    phenoName <- phenoName.vec[i]
    df.results <- Fit_Model_Int(1,nrow(results))
    # df.results <- Fit_Model_Int(1,14)
    df.results$SNP <- results$SNP
    df.results$E <- results$E
    df.results$pheno <- phenoName
    if (i==1) {
      df.results.save <- df.results
    } else {
      df.results.save <- rbind(df.results.save,df.results)
    }
  }
  
  f<-paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.',s,'.pleiotropy.txt')
  fwrite(df.results.save,f,quote = F,sep = '\t',na = 'NA',row.names = F,col.names = T)
  print(f)
}

s='80'; results.80 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.',s,'.pleiotropy.txt'),data.table = F,stringsAsFactors = F)
s='20'; results.20 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/',pheno,'.GxE.',s,'.pleiotropy.txt'),data.table = F,stringsAsFactors = F)
df.results.save <- merge(results.80,results.20,by=c('SNP','E','pheno'))
df.results.save$FDR <- p.adjust(df.results.save[,5],method = 'fdr')
# df.results.save$FDR <- p.adjust(df.results.save[,2],method = 'fdr')
subset(df.results.save,FDR<0.1)

