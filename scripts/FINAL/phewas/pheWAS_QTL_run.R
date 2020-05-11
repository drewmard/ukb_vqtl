library(data.table)
source('/athena/elementolab/scratch/anm2868/open_targets/scripts/open_targets_query_func.R')
f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/sig_results/bmi.SNP_sig_list.query.txt')
query <- fread(f,data.table = F,stringsAsFactors = F,header = F,fill=T)
query <- subset(query,query[,2]!='')
# for (i in 1:nrow(query)) {
# for (i in 1:50) {
#   print(i)
#   snp <- query[i,2]
#   tryCatch({
#     res <- pheWAS(snp)
#     res.df <- do.call(rbind,lapply(res,function(x) {data.frame(rs=query[i,1],SNP=query[i,2],trait=x$study$traitReported,category=x$study$traitCategory,beta=x$beta,pval=x$pval,ntotal=x$nTotal,nCases=ifelse(is.null(x$nCases),'NA',x$nCases),stringsAsFactors=FALSE)}))
#     if (i==1) {
#       res.df.save <- res.df
#     } else {
#       res.df.save <- rbind(res.df.save,res.df)
#     }
#   }, error=function(e) {NULL}
#   )
# }

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
res.df <- mclapply(1:nrow(query),pheWAS_search_df,mc.cores = 8)
res.df.save <- do.call(rbind,res.df)

f.out <- paste0("/athena/elementolab/scratch/anm2868/open_targets/query/bmi.QTL.SNP_sig_list.query.pheWAS.txt")
fwrite(res.df.save,file=f.out,col.names = T,row.names = F,sep = '\t',na = 'NA',quote = F)
