library(data.table)
library(BEDMatrix)
# read in CAD
f<-'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/CAD.txt'
CAD <- fread(f,data.table = F,stringsAsFactors = F)
f<-'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/Diabetes.txt'
Diabetes <- fread(f,data.table = F,stringsAsFactors = F)
f<-'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/HBP.txt'
HBP <- fread(f,data.table = F,stringsAsFactors = F)

pheno='bmi'; 
results.df.lst <- list()
for (s in c('80','20')) { 
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
  
  
  i=38
  phenoName='Diabetes'
  SNP <- results$SNP[i]
  envir_name <- results$E[i]
  index <- which(SNP_names == SNP)
  full_dataset$SNP <- geno[ind,index]
  # mod.formula <- formula(paste(phenoName,' ~ age+age2+genotyping.array+sex+age*sex+age2*sex+
  #                PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
  #                PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+',
  #                              envir_name,'*SNP'))
  # mod <- glm(mod.formula,data=full_dataset,family=binomial(link="logit"))
  # res <- summary(mod)$coef[paste0(envir_name,':SNP'),c("Estimate",'Std. Error',"Pr(>|z|)")]
  # mod.formula <- formula(paste(phenoName,' ~ age+age2+genotyping.array+sex+age*sex+age2*sex+
  #                PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
  #                PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+',
  #                              envir_name,'*SNP+bmi.na'))
  # mod <- glm(mod.formula,data=full_dataset,family=binomial(link="logit"))
  # res <- summary(mod)$coef[paste0(envir_name,':SNP'),c("Estimate",'Std. Error',"Pr(>|z|)")]
  
  results.df <- list()
  for (Val in 1:3) {
    phenoName='Diabetes'
    # mod.formula <- formula(paste(phenoName,' ~ age+age2+genotyping.array+sex+age*sex+age2*sex+
    #            PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
    #            PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+SNP'))
    mod.formula <- formula(paste(phenoName,' ~ age+age2+genotyping.array+sex+age*sex+age2*sex+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
               PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+SNP+bmi.na'))
    
    mod <- glm(mod.formula,data=subset(full_dataset,PA==Val),family=binomial(link="logit"))
    res <- summary(mod)$coef['SNP',c("Estimate",'Std. Error',"Pr(>|z|)")]
    
    results.df[[Val]] <- data.frame(pheno=phenoName,
                                    E=envir_name,
                                    Val=Val,
                                    Set=s,
                                    EST=as.numeric(exp(res[1])),
                                    LOW=as.numeric(exp(res[1]-1.96*res[2])),
                                    HIGH=as.numeric(exp(res[1]+1.96*res[2])),
                                    P=as.numeric(res[3]))
  }
  for (Val in 1:3) {
    phenoName <- 'bmi.na'
    mod.formula <- formula(paste(phenoName,' ~ age+age2+genotyping.array+sex+age*sex+age2*sex+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
               PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+SNP'))
    mod <- lm(mod.formula,data=subset(full_dataset,PA==Val))
    res <- summary(mod)$coef['SNP',c("Estimate",'Std. Error',"Pr(>|t|)")]
    results.df[[Val+3]] <- data.frame(pheno=phenoName,
                                      E=envir_name,
                                      Val=Val,
                                      Set=s,
                                      EST=as.numeric((res[1])),
                                      LOW=as.numeric((res[1]-1.96*res[2])),
                                      HIGH=as.numeric((res[1]+1.96*res[2])),
                                      P=as.numeric(res[3]))
    
  }
  results.df.lst[[s]] <- do.call(rbind,results.df)
}

results.df.save <- as.data.frame(do.call(rbind,results.df.lst))
fwrite(results.df.save,'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/diabetes_bmi_marginal.used_bmi.txt',quote = F,na='NA',sep = '\t',row.names = F,col.names = T)
# fwrite(results.df.save,'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/diabetes_bmi_marginal.txt',quote = F,na='NA',sep = '\t',row.names = F,col.names = T)

# mod.formula <- formula(paste(phenoName,' ~ age+age2+genotyping.array+sex+age*sex+age2*sex+
#                PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
#                PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+SNP'))
# mod <- glm(mod.formula,data=full_dataset,family=binomial(link="logit"))
# res <- summary(mod)$coef['SNP',c("Estimate",'Std. Error',"Pr(>|t|)")]

# cor.test(full_dataset$SNP,full_dataset$PA)

# do the same for the 20% set.


