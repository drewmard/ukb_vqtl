library(data.table)

LD_statistics <- function(data,SNP1,SNP2) {
  
  df2 <- data[,c(3,10:ncol(df))]
  rownames(df2) <- df2$ID; df2 <- df2[,-1]; df2 <- as.data.frame(t(df2)); df2$eid <- rownames(df2); rownames(df2) <- NULL
  
  f <- '/home/kulmsc/athena/ukbiobank/haplotypes/ukb_hap_chr6_v2.sample'
  samp_id <- fread(f,data.table=F,stringsAsFactors = F)
  samp_id <- samp_id[-1,]
  
  df2$eid <- samp_id$ID_2
  
  f <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/phenotypes_processed.80.txt'
  pheno <- fread(f,data.table = F,stringsAsFactors = F)
  pheno2 <- pheno[,c('IID','lymphocyte.count.rint.ALL')]
  df.mg <- merge(df2,pheno2,by.x='eid',by.y='IID')
  df.mg[which(df.mg[,'lymphocyte.count.rint.ALL']==-9),'lymphocyte.count.rint.ALL'] <- NA
  df.mg <- subset(df.mg,(!is.na(lymphocyte.count.rint.ALL)))
  
  i1 <- grep(SNP1,colnames(df.mg))
  i2 <- grep(SNP2,colnames(df.mg))
  
  df.mg[,i1] <- as.character(df.mg[,i1])
  x <- strsplit(unlist(lapply(strsplit(df.mg[,i1],':'),function(x) x[1])),'|')
  df.mg[,paste0(SNP1,'_1')] <- unlist(lapply(x,function(x) x[1]))
  df.mg[,paste0(SNP1,'_2')] <- unlist(lapply(x,function(x) x[3]))
  
  df.mg[,i2] <- as.character(df.mg[,i2])
  x <- strsplit(unlist(lapply(strsplit(df.mg[,i2],':'),function(x) x[1])),'|')
  df.mg[,paste0(SNP2,'_1')] <- unlist(lapply(x,function(x) x[1]))
  df.mg[,paste0(SNP2,'_2')] <- unlist(lapply(x,function(x) x[3]))

  df.mg[,'Haplotype1'] <- paste0(df.mg[,paste0(SNP1,'_1')],df.mg[,paste0(SNP2,'_1')])
  df.mg[,'Haplotype2'] <- paste0(df.mg[,paste0(SNP1,'_2')],df.mg[,paste0(SNP2,'_2')])

  pA = (
    sum(as.numeric(df.mg[,paste0(SNP1,'_1')])) + 
      sum(as.numeric(df.mg[,paste0(SNP1,'_2')]))
    ) / 
    (2*nrow(df.mg))
  pB = (
    sum(as.numeric(df.mg[,paste0(SNP2,'_1')])) + 
      sum(as.numeric(df.mg[,paste0(SNP2,'_2')]))
  ) / 
    (2*nrow(df.mg))
  
  pAB = (sum(df.mg$Haplotype1=='11') + sum(df.mg$Haplotype2=='11'))/(2*nrow(df.mg))
  
  D=pAB - pA*pB
  if (D < 0) { Dmax=max(-1*(1-pA)*(1-pB),-1*pA*pB) } else { Dmax=min(pA*(1-pB),(1-pA)*pB) }
  Dprime=D/Dmax;
  # Dprime
  
  r = D / sqrt((pA*(1-pA)*pB*(1-pB)));# r
  r2 = D^2 / (pA*(1-pA)*pB*(1-pB)); r2;# r^2
  
  print(paste0("D' = ",Dprime));
  print(paste0("r2 = ", r2));
}

# rm(data);rm(SNP1);rm(SNP2)
# 
# f <- '/athena/elementolab/scratch/anm2868/vQTL/UKB/subset/phased/ukb_hap_chr6_v2.subset3.vcf'
# df <- fread(f,data.table = F,stringsAsFactors = F)
# SNP1='rs1265158'
# SNP2='rs2516491'
# SNP3='rs3131643'
# LD_statistics(data=df,SNP1=SNP1,SNP2=SNP2)
# LD_statistics(data=df,SNP1=SNP1,SNP2=SNP3)
# LD_statistics(data=df,SNP1=SNP2,SNP2=SNP3)

