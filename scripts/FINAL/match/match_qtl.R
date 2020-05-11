library(data.table)
f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/andrew_copies/maf/ukbb.','ALL','.impute.frqx')
df.frqx <- fread(f,data.table = F,stringsAsFactors = F)

f <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/sig_results/bmi.SNP_sig_list2.txt'
df <- fread(f,data.table = F,stringsAsFactors = F,header = F)

df.sub <- subset(df.frqx,SNP%in%df[,1])
# 
# f <- '/athena/elementolab/scratch/anm2868/open_targets/query/bmi.mean.txt'
# df <- fread(f,data.table = F,stringsAsFactors = F,header = F)
# 
# mean.raw <- subset(df.frqx,SNP%in%df[,1])
# mean.raw <- subset(mean.raw,!(SNP%in%var.raw[,1]))
# t.test(mean.raw[,2],var.raw[,2])
# t.test(mean.raw[,3],var.raw[,3])
# t.test(mean.raw[,4],var.raw[,4])

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
# res.df <- mclapply(1:500,pheWAS_search_df,mc.cores = 8)
res.df <- mclapply(1:nrow(query),pheWAS_search_df,mc.cores = 8)
# res.df <- lapply(1:2,pheWAS_search_df)
res.df.save <- do.call(rbind,res.df)

f.out <- paste0("/athena/elementolab/scratch/anm2868/open_targets/query/bmi.QTL.matched_snp.query.pheWAS.txt")
fwrite(res.df.save,file=f.out,col.names = T,row.names = F,sep = '\t',na = 'NA',quote = F)

#####################################################

library(data.table)
# field <- 'category'
field <- 'trait'
PVAL_THRES <- 0.05
SAMP_SIZE_MIN <- 10000
pheno='bmi'
QTLtype='var.raw'
# QTLtype='mean.raw'
# background_set <- 'matched_snp'
background_set <- 'mean'
category_keep <- fread('/athena/elementolab/scratch/anm2868/open_targets/category_keep.csv',data.table = F,stringsAsFactors = F,sep=',',header = F)
# QTL data
# f <- paste0('/athena/elementolab/scratch/anm2868/open_targets/query/bmi.var.raw.query.txt')
# n=nrow(fread(f,data.table = F,stringsAsFactors = F))
f <- paste0("/athena/elementolab/scratch/anm2868/open_targets/query/",pheno,".",QTLtype,".matched_snp.query.pheWAS.txt")
qtl.full <- fread(f,data.table = F,stringsAsFactors = F)
qtl.full <- subset(qtl.full,SNP %in% unique(qtl.full$SNP)[seq(1,length(unique(qtl.full$SNP)),by=9)])
# f <- paste0("/athena/elementolab/scratch/anm2868/open_targets/query/",pheno,".",QTLtype,".query.pheWAS.txt")
# qtl.full <- fread(f,data.table = F,stringsAsFactors = F)
qtl <- subset(qtl.full,pval < PVAL_THRES)
qtl <- subset(qtl,!duplicated(qtl[,c('rs',field)]))
qtl <- subset(qtl,category %in% category_keep[,1])
qtl <- subset(qtl,ntotal > SAMP_SIZE_MIN)
n <- length(unique(qtl$rs))
qtl <- table(qtl[,field])/n
qtl <- data.frame(qtl)

# background snps, eg mean or matched
# f <- paste0("/athena/elementolab/scratch/anm2868/open_targets/query/bmi.var.raw.matched_snp.query.txt")
# m=nrow(subset(fread(f,data.table = F,stringsAsFactors = F,fill=T,header=F),V2!=''))
# f <- paste0("/athena/elementolab/scratch/anm2868/open_targets/query/bmi.var.raw.matched_snp.query.pheWAS.txt")
if (background_set=='matched_snp') {
  f <- paste0("/athena/elementolab/scratch/anm2868/open_targets/query/",pheno,".",QTLtype,".matched_snp.query.pheWAS.txt")
  df.full <- fread(f,data.table = F,stringsAsFactors = F)
  suff='bg_matched_snp'
} else if (background_set=='mean') {
  f <- paste0("/athena/elementolab/scratch/anm2868/open_targets/query/bmi.mean.raw.query.pheWAS.txt")
  df.full <- fread(f,data.table = F,stringsAsFactors = F)
  suff='bg_mean'
}
df <- subset(df.full,!(rs %in% qtl.full$rs))
df <- subset(df,pval < PVAL_THRES)
df <- subset(df,!duplicated(df[,c('rs',field)]))
df <- subset(df,category %in% category_keep[,1])
# x <- unique(df[,c('category','trait')])
# x <- x[order(x$category,x$trait),]
# fwrite(x,"/athena/elementolab/scratch/anm2868/open_targets/query/trait_category.txt",quote = F,row.names = F,col.names = T,na='NA',sep = '\t')
m <- length(unique(df$rs))
df <- table(df[,field])/m
df <- data.frame(df)


colnames(df)[2] <- 'Expected'
colnames(qtl)[2] <- 'Observed'

df.mg <- merge(df,qtl,by='Var1',all=T)
df.mg[is.na(df.mg)] <- 0
df.mg$Diff <- df.mg$Observed - df.mg$Expected
df.mg$Diff.Ratio <- df.mg$Observed/df.mg$Expected
# df.mg[order(df.mg$Diff,decreasing = T)[1:10],]
# df.mg[order(df.mg$Diff.Ratio,decreasing = T)[1:10],]
# df.mg[order(df.mg$Observed,decreasing = T)[1:10],]
# df.mg[order(df.mg$Expected,decreasing = F)[1:10],]

df.mg2 <- subset(df.mg,Observed > 0)
# df.mg2 <- subset(df.mg,Expected > 0)
p.vec <- c()
for (i in 1:nrow(df.mg2)) {
  x <- df.mg2[i,]
  p <- binom.test(x$Observed*n,n,x$Expected,alternative = 'greater')$p.value
  # p <- binom.test(x$Observed*n,n,x$Expected,alternative = 'less')$p.value
  p.vec <- c(p.vec,p)
}
df.mg2$PVAL <- p.vec
# df.mg2[order(df.mg2$PVAL,-df.mg2$Observed,decreasing = F)[1:5],]
# df.mg2[order(df.mg2$Observed,decreasing = T)[1:5],]
df.mg2$FDR <- p.adjust(df.mg2$PVAL,method = 'fdr')
# df.mg2.sub[order(df.mg2.sub$FDR,decreasing = F)[1:5],]
f <- paste0("/athena/elementolab/scratch/anm2868/open_targets/query/",pheno,".",QTLtype,".query.pheWAS.txt")
qtl.full <- fread(f,data.table = F,stringsAsFactors = F)
qtl.full <- unique(qtl.full[,c('trait','category')])
df.mg2 <- merge(qtl.full,df.mg2,by.y='Var1',by.x='trait')
df.mg2.sub <- subset(df.mg2,FDR < 0.1 & Observed > 0)#.5)

# f.out=paste0("/athena/elementolab/scratch/anm2868/open_targets/query/",pheno,".",QTLtype,".query.pheWAS.",suff,".sig.txt")
# fwrite(df.mg2.sub,f.out,quote = F,na = 'NA',sep = '\t',row.names = F,col.names = T)
# f.out=paste0("/athena/elementolab/scratch/anm2868/open_targets/query/",pheno,".",QTLtype,".query.pheWAS.",suff,".txt")
# fwrite(df.mg2,f.out,quote = F,na = 'NA',sep = '\t',row.names = F,col.names = T)
# f.out=paste0("/athena/elementolab/scratch/anm2868/open_targets/query/",pheno,".",QTLtype,".matched_snp.query.pheWAS.",suff,".sig.txt")
# fwrite(df.mg2.sub,f.out,quote = F,na = 'NA',sep = '\t',row.names = F,col.names = T)
# f.out=paste0("/athena/elementolab/scratch/anm2868/open_targets/query/",pheno,".",QTLtype,".matched_snp.query.pheWAS.",suff,".txt")
# fwrite(df.mg2,f.out,quote = F,na = 'NA',sep = '\t',row.names = F,col.names = T)
f.out=paste0("/athena/elementolab/scratch/anm2868/open_targets/query/",pheno,".",QTLtype,".matched_snp_1.query.pheWAS.",suff,".sig.txt")
fwrite(df.mg2.sub,f.out,quote = F,na = 'NA',sep = '\t',row.names = F,col.names = T)
f.out=paste0("/athena/elementolab/scratch/anm2868/open_targets/query/",pheno,".",QTLtype,".matched_snp_1.query.pheWAS.",suff,".txt")
fwrite(df.mg2,f.out,quote = F,na = 'NA',sep = '\t',row.names = F,col.names = T)

# 
# library(data.table)
# # field <- 'category'
# field <- 'trait'
# PVAL_THRES <- 0.05
# pheno='bmi'
# QTLtype='var.raw'
# # QTLtype='mean.raw'
# background_set <- 'matched_snp'
# suff='bg_mean'
# f=paste0("/athena/elementolab/scratch/anm2868/open_targets/query/",pheno,".",QTLtype,".query.pheWAS.",suff,".txt")
# df <- fread(f,data.table = F,stringsAsFactors = F)
# 
# f <- paste0("/athena/elementolab/scratch/anm2868/open_targets/query/",pheno,".",QTLtype,".query.pheWAS.txt")
# qtl.full <- fread(f,data.table = F,stringsAsFactors = F)
# qtl.full <- unique(qtl.full[,c('trait','category')])
# df.mg <- merge(df,qtl.full,by.x='Var1',by.y='trait')
# subset(df.mg,category=="Immune system")
# subset(df.mg,category=="Hematological measurement")
# subset(df.mg,category=="Endocrine system")
# 
# subset(df.mg2.sub,Var1=="Diabetes | non-cancer illness code, self-reported")

