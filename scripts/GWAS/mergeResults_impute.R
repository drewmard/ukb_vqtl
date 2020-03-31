library(data.table)
phenoName='bmi.ALL'


for (CHR in 1:22) {
  print(paste('CHR:',CHR))

  print('GWAS...')
  df.mean <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/imputed/results/ukbb.',CHR,'.impute.',phenoName,'.muGWAS.qassoc'),data.table = F,stringsAsFactors = F)
  
  # df.SNP_data <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/andrew_copies/ukbb.',CHR,'.impute.bim'),data.table = F,stringsAsFactors = F)
  # df.SNP_data <- df.SNP_data[,c(2,4)]
  # colnames(df.SNP_data) <- c('SNP','POS')
  
  # freq
  print('Freq...')
  df.freq <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/imputed/MAF/ukbb.',CHR,'.impute.frq'),data.table = F,stringsAsFactors = F)
  df.freq <- df.freq[,c('SNP','A1','A2','MAF')]
  
  print('Merge...')
  df.mg2 <- merge(df.mean,df.freq,by='SNP')
    
  # print('Saving...')
  if (CHR==1) {
    df.save <- df.mg2
  } else {
    df.save <- rbind(df.save,df.mg2)
  }
  
}

f.save <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/imputed/results/ukbb.',phenoName,'.results.txt')
fwrite(df.save,f.save,col.names = T,row.names = F,na='NA',quote=F,sep='\t')

