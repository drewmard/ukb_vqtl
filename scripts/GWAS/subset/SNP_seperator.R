library(data.table)

for (i in 1:22) {
  df <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/andrew_copies/ukbb.',i,'.impute.bim'),data.table = F,stringsAsFactors = F)
  rng <- seq(1,nrow(df),by = 5000)
  for (j in 1:length(rng)) {
    if (j==length(rng)) {
      x <- df[rng[j]:nrow(df),2]
    } else {
      x <- df[rng[j]:(rng[j+1]-1),2]
    }
    fwrite(data.frame(x=x),paste0('/athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/andrew_copies/subset/SNP.',i,'.',j,'.impute.txt'),row.names = F,col.names = F,quote = F,na='NA',sep = '\t')
  }
}

flist <- list.files('/athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/andrew_copies/subset/')
x <- strsplit(flist,'\\.')
x1 <- unlist(lapply(x,function(x) x[2]))
x2 <- unlist(lapply(x,function(x) x[3]))
fwrite(data.frame(x1,x2),paste0('/athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/andrew_copies/subset/ID.impute.txt'),row.names = F,col.names = F,quote = F,na='NA',sep = '\t')




