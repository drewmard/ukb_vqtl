library(data.table)
phenoName='lymphocyte.count.log'

for (CHR in 1:22) {
  print(paste('CHR:',CHR))
  RELEVANT_COLUMNS=c('SNP','CHR','BP','NMISS','BETA.x','P.x','BETA.y','P.y','REF','ALT')
  #mean
  #print('Mean...')
  df.mean <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/UKB/results/ukbb.',CHR,'.',phenoName,'.muGWAS.qassoc'),data.table = F,stringsAsFactors = F)
  
  #var
  #print('Var...')
  tryCatch({
    df.var <- read.table(paste0('/athena/elementolab/scratch/anm2868/vQTL/UKB/results/ukbb.',CHR,'.',phenoName,'.vGWAS.auto.R'),stringsAsFactors = F)
    colnames(df.var) <- c('CHR','SNP','BP','REF',
                          'BETA','SE','T','P')
    
    # freq
    #print('Freq...')
    df.freq <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/UKB/results/ukbb.',CHR,'.frq'),data.table = F,stringsAsFactors = F)
    colnames(df.freq)[4] <- 'ALT'
    
    #print('Merge...')
    df.mg <- merge(df.mean,df.var,by=c('SNP','CHR','BP'))
    df.mg2 <- merge(df.mg,df.freq[,c('SNP','ALT','MAF')])
    
  },error=function(e) { print(paste0('ERROR: ',CHR)) })
  
  # print('Saving...')
  if (CHR==1) {
    df.save <- df.mg2
  } else {
    df.save <- rbind(df.save,df.mg2)
  }
  
}

f.save <- paste0('/athena/elementolab/scratch/anm2868/vQTL/UKB/results/ukbb.',phenoName,'.results.txt')
fwrite(df.save,f.save,col.names = T,row.names = F,na='NA',quote=F,sep='\t')

