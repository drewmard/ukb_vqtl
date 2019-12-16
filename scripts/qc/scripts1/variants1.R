library(data.table)
# chrNum <- 18
args = commandArgs(trailingOnly=TRUE)
chrNum <- args[1]
# for (chrNum in 3:22) {
  print(chrNum)
  UKB_info_MAF.sub <- fread(paste0('/home/kulmsc/athena/ukbiobank/qc/imputed/ukb_mfi_chr',chrNum,'_v3.txt'),data.table = F,stringsAsFactors = F)
  
  UKB_info_MAF.sub <- UKB_info_MAF.sub[,c(2,6,8)]
  colnames(UKB_info_MAF.sub) <- c('SNP','MAF','INFO')
  UKB_info_MAF.sub <- subset(UKB_info_MAF.sub,MAF > 0.05 & INFO > 0.3)
  
  fwrite(data.frame(UKB_info_MAF.sub$SNP),
         paste0('/athena/elementolab/scratch/anm2868/vQTL/UKB/imputed/variants1/ukb.',chrNum,'.imputed.variants1.txt'),
         sep = '\t',quote = F,na = 'NA',row.names = F,col.names = F)
  rm(UKB_info_MAF.sub)
# }
