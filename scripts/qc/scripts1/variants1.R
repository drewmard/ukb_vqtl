library(data.table)
split_num=20
# chrNum <- 18
args = commandArgs(trailingOnly=TRUE)
chrNum <- args[1]
# for (chrNum in 3:22) {
  print(chrNum)
  UKB_info_MAF.sub <- fread(paste0('/home/kulmsc/athena/ukbiobank/qc/imputed/ukb_mfi_chr',chrNum,'_v3.txt'),data.table = F,stringsAsFactors = F)
  
  UKB_info_MAF.sub <- UKB_info_MAF.sub[,c(2,6,8)]
  colnames(UKB_info_MAF.sub) <- c('SNP','MAF','INFO')
  UKB_info_MAF.sub <- subset(UKB_info_MAF.sub,MAF > 0.05 & INFO > 0.3)
  
  df.tmp <- data.frame(UKB_info_MAF.sub$SNP)
  x=1:nrow(df.tmp)
  max=length(x)/split_num
  lst <- split(x,ceiling(x/max))
  
  for (j in 1:split_num) {
    df.tmp2 <- df.tmp[lst[[j]],]
    fwrite(as.data.frame(df.tmp2),
           paste0('/athena/elementolab/scratch/anm2868/vQTL/UKB/imputed/variants1/ukb.',chrNum,'.imputed.variants1.',j,'.txt'),
           sep = '\t',quote = F,na = 'NA',row.names = F,col.names = F)
    
  }
  
  # fwrite(data.frame(UKB_info_MAF.sub$SNP),
  #        paste0('/athena/elementolab/scratch/anm2868/vQTL/UKB/imputed/variants1/ukb.',chrNum,'.imputed.variants1.',iter,'.txt'),
  #        sep = '\t',quote = F,na = 'NA',row.names = F,col.names = F)
  rm(UKB_info_MAF.sub)
# }
