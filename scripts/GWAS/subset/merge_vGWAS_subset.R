library(data.table)
ID.df <- fread('/athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/andrew_copies/subset/ID.impute.txt',data.table = F,stringsAsFactors = F)
phenotype='lymphocyte.count.rint.ALL'

for (i in 1:nrow(ID.df)) {
# for (i in 1:5) {
  if (i %% 100 == 0) {print(i)}
    
  x1=ID.df[i,1]
  x2=ID.df[i,2]
  
  tryCatch({
    f.res.tmp <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_subset/','ukbb.',x1,'.',x2,'.',phenotype,'.txt')
    df.res.tmp <- fread(f.res.tmp,data.table = F,stringsAsFactors = F)

    if (i==1) {
      df.res.save <- df.res.tmp
    } else {
      df.res.save <- rbind(df.res.save,df.res.tmp)
    }
    
  },error=function(e) { print(paste0('ERROR: ',i,'(',x1,',',x2,')')) })
}


colnames(df.res.save)[4] <- 'P'
f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_subset/','ukbb.',phenotype,'.vGWAS.txt')
print(paste0('Writing: ',f.out))
fwrite(df.res.save,f.out,sep='\t',quote=F,col.names = T,row.names = F,na="NA")


