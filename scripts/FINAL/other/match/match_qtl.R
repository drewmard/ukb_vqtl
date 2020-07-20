library(data.table)
f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/andrew_copies/maf/ukbb.','ALL','.impute.frqx')
df.frqx <- fread(f,data.table = F,stringsAsFactors = F)

f <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/sig_results/bmi.SNP_sig_list2.txt'
df <- fread(f,data.table = F,stringsAsFactors = F,header = F)

df.sub <- subset(df.frqx,SNP%in%df[,1])

# match
ind_lst <- list()
nperm <- 10
for (j in 1:nrow(df.sub)) {
  print(j)
  theta <- 0.01
  ind <- c()
  while (length(ind) < nperm) {
    hom_a1.theta <- theta
    MAF.theta <- theta
    MISS.theta <- df.sub$MISS[j] * theta
    ind <- which(df.frqx$hom_a1 >= (df.sub$hom_a1[j] - hom_a1.theta) &
                   df.frqx$hom_a1 <= (df.sub$hom_a1[j] + hom_a1.theta) &
                   df.frqx$MAF >= (df.sub$MAF[j] - MAF.theta) &
                   df.frqx$MAF <= (df.sub$MAF[j] + MAF.theta) &
                   df.frqx$MISS >= (df.sub$MISS[j] - MISS.theta) &
                   df.frqx$MISS <= (df.sub$MISS[j] + MISS.theta) &
                   df.frqx$CHR != df.sub$CHR[j]
                 )
    theta <- theta + 0.01
  }
  df.frqx.sub <- df.frqx[ind,]
  m=nrow(df.frqx.sub); print(m)
  perm_snp <- df.frqx.sub[sample(1:m,nperm,replace=FALSE),'SNP']
  ind_lst[[j]] <- perm_snp
}
ind_lst.2 <- unlist(ind_lst)

f.out <- paste0('/athena/elementolab/scratch/anm2868/open_targets/query/bmi.QTL.matched_snp.txt')
fwrite(data.frame(ind_lst.2),f.out,row.names = F,col.names = F,sep = '\t',na = 'NA',quote = F)

##################
input=/athena/elementolab/scratch/anm2868/open_targets/query/bmi.QTL.matched_snp.txt
output=/athena/elementolab/scratch/anm2868/open_targets/query/bmi.QTL.matched_snp.query.txt
rm $output
while read -r SNP; do
echo $SNP
if [[ "$SNP" == *":"* ]]; then
SNP=$(echo $SNP | sed 's/:/_/g')
SNP=${SNP}_b37
else
  SNP=${SNP},
fi
query=$(grep -m 1 $SNP /athena/elementolab/scratch/anm2868/open_targets/csv_files/variant_index.csv)
echo $query >> $output
done < $input
###################

# spack load -r r@3.5.0
library(data.table)
source('/athena/elementolab/scratch/anm2868/open_targets/scripts/open_targets_query_func.R')
f <- paste0('/athena/elementolab/scratch/anm2868/open_targets/query/bmi.QTL.matched_snp.query.txt')
query <- fread(f,data.table = F,stringsAsFactors = F,header = F,fill=T)
query <- subset(query,query[,2]!='')

pheWAS_search_df <- function(i) {
  snp <- query[i,2]
  tryCatch({
    res <- pheWAS(snp)
    res.df <- do.call(rbind,
                      lapply(res,function(x) {
                        data.frame(rs=query[i,1],
                                   SNP=query[i,2],
                                   trait=x$study$traitReported,
                                   category=x$study$traitCategory,
                                   beta=x$beta,pval=x$pval,
                                   ntotal=x$nTotal,
                                   nCases=ifelse(is.null(x$nCases),
                                                 'NA',
                                                 x$nCases),
                                   stringsAsFactors=FALSE)
                      }))
  }, error=function(e) {
    res.df <- data.frame(rs=query[i,1],
                         SNP=query[i,2],
                         trait=NA,
                         category=NA,
                         beta=NA,
                         pval=NA,
                         ntotal=NA,
                         nCases=NA,
                         stringsAsFactors=FALSE)
  })
  return(res.df)
}
library(parallel)
# res.df <- mclapply(1:500,pheWAS_search_df,mc.cores = 8)
res.df <- mclapply(1:nrow(query),pheWAS_search_df,mc.cores = 8)
# res.df <- lapply(1:2,pheWAS_search_df)
res.df.save <- do.call(rbind,res.df)

f.out <- paste0("/athena/elementolab/scratch/anm2868/open_targets/query/bmi.QTL.matched_snp.query.pheWAS.txt")
fwrite(res.df.save,file=f.out,col.names = T,row.names = F,sep = '\t',na = 'NA',quote = F)

