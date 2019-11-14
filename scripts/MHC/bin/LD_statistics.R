library(data.table)

LD_statistics <- function(df.mg,SNP1,SNP2) {
  
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

