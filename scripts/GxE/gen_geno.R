library(data.table)

# initialize
user_direc <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl'
phenoName <- 'lymphocyte.count.rint'

# var hits
f <- paste0(user_direc,'/output/GWAS/results2/',phenoName,'/ukbb.',phenoName,'.results.var.clumped.cut.txt')
index <- fread(f,data.table = F,stringsAsFactors = F)

# mean hits
# f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/UKB/results/ukbb.',phenoName,'.results.mean.clumped.cut.txt')
# index <- fread(f,data.table = F,stringsAsFactors = F)

# index <- fread('/athena/elementolab/scratch/anm2868/vQTL/UKB/results/ukbb.bmi.results.custom.txt',data.table = F,stringsAsFactors = F)
print('Creating genotype file...')
for (i in 1:nrow(index)) {
  print(paste0(i,'/',nrow(index),' SNPs...'))
  vQTL=index[i,2]
  CHR_vQTL=index[i,1]
  tryCatch(
    {
      f.vQTL <- paste0(user_direc,'/output/GWAS/subset/',phenoName,'/ukbb.',CHR_vQTL,'.',vQTL,'.raw')
      df.vQTL <- fread(f.vQTL,data.table = F,stringsAsFactors = F)
      
      df.vQTL <- df.vQTL[,c(2,grep(vQTL,colnames(df.vQTL)))]
      colnames(df.vQTL)[2] <- vQTL
      
      if (i != 1) {
        df.geno <- merge(df.geno,df.vQTL,by='IID')
      } else {
        df.geno <- df.vQTL
      }
    },error=function(e) {
      print(paste0('ERROR: ',i))
      break
    }
  )
  
}

f.out <- paste0(user_direc,'/output/GWAS/subset/',phenoName,'/ukbb.ALL_vQTL.raw')
fwrite(df.geno,f.out,sep='\t',quote=F,row.names = F,col.names = T,na='NA')
