# Run pheWAS enrichment analysis
pheWAS_enrichment_analysis <- function(qtl.full,df.full) { 
  qtl <- subset(qtl.full,pval <= PVAL_THRES)
  qtl <- subset(qtl,!duplicated(qtl[,c('rs',field)]))
  qtl <- subset(qtl,category %in% category_keep[,1])
  qtl <- subset(qtl,ntotal > SAMP_SIZE_MIN)
  n <- length(unique(qtl$rs))
  qtl <- table(qtl[,field])/n
  qtl <- data.frame(qtl)
  
  # df <- subset(df.full,!(rs %in% qtl.full$rs))
  df <- subset(df.full,pval <= PVAL_THRES)
  df <- subset(df,!duplicated(df[,c('rs',field)]))
  df <- subset(df,category %in% category_keep[,1])
  m <- length(unique(df$rs))
  df <- table(df[,field])/m
  df <- data.frame(df)
  
  colnames(df)[2] <- 'Expected'
  colnames(qtl)[2] <- 'Observed'
  
  df.mg <- merge(df,qtl,by='Var1',all=T)
  df.mg[is.na(df.mg)] <- 0
  df.mg$Diff <- df.mg$Observed - df.mg$Expected
  df.mg$Diff.Ratio <- df.mg$Observed/df.mg$Expected
  
  df.mg2 <- subset(df.mg,Observed > 0)
  df.mg2 <- subset(df.mg,Observed > 0 & Expected > 0)
  p.vec <- c()
  for (i in 1:nrow(df.mg2)) {
    x <- df.mg2[i,]
    p <- binom.test(x$Observed*n,n,x$Expected,alternative = 'greater')$p.value
    p.vec <- c(p.vec,p)
  }
  df.mg2$PVAL <- p.vec
  df.mg2$FDR <- p.adjust(df.mg2$PVAL,method = 'fdr')
  return(df.mg2)
}

library(data.table)

category_keep <- fread('/athena/elementolab/scratch/anm2868/open_targets/category_keep.csv',data.table = F,stringsAsFactors = F,sep=',',header = F)
# field <- 'category'
field <- 'trait'
PVAL_THRES <- 0.05
SAMP_SIZE_MIN <- 10000

# Query data
f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/sig_results/bmi.SNP_sig_list.query.txt')
query <- fread(f,data.table = F,stringsAsFactors = F,header = F,fill=T)

# sig results data
sig_results <- fread('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/sig_results/bmi.sig.txt',data.table = F,stringsAsFactors = F)
# vqtls vs muqtls. rint vqtls vs vqtls. dqtls vs rint vqtls. rint vqtls vs muqtls. dqtls vs muqtls.

# QTL data
f <- paste0("/athena/elementolab/scratch/anm2868/open_targets/query/bmi.QTL.SNP_sig_list.query.pheWAS.txt")
qtl.full <- fread(f,data.table = F,stringsAsFactors = F)
ind <- which(sig_results$Raw.vQTL==1)
# ind <- which(sig_results$Rint.vQTL==1)
# ind <- which(sig_results$dQTL==1)

qtl.full <- subset(qtl.full,rs %in% query[ind,1])

# Background data
# f <- paste0("/athena/elementolab/scratch/anm2868/open_targets/query/bmi.QTL.matched_snp.query.pheWAS.txt")
# df.full <- fread(f,data.table = F,stringsAsFactors = F)

f <- paste0("/athena/elementolab/scratch/anm2868/open_targets/query/bmi.QTL.SNP_sig_list.query.pheWAS.txt")
df.full <- fread(f,data.table = F,stringsAsFactors = F)
ind <- which(sig_results$Raw.vQTL==0 & sig_results$Mean.QTL==1)
# ind <- which(sig_results$Raw.vQTL==1)
df.full <- subset(df.full,rs %in% query[ind,1])


df.mg2 <- pheWAS_enrichment_analysis(qtl.full,df.full)
df.mg2[order(df.mg2$FDR)[1:10],]


# raw vQTLs vs pure muQTLs
f <- paste0("/athena/elementolab/scratch/anm2868/open_targets/query/bmi.QTL.SNP_sig_list.query.pheWAS.txt")
qtl.full <- fread(f,data.table = F,stringsAsFactors = F)
ind1 <- which(sig_results$Raw.vQTL==1)
f <- paste0("/athena/elementolab/scratch/anm2868/open_targets/query/bmi.QTL.SNP_sig_list.query.pheWAS.txt")
df.full <- fread(f,data.table = F,stringsAsFactors = F)
ind2 <- which(sig_results$Raw.vQTL==0 & sig_results$Mean.QTL==1)
res <- pheWAS_enrichment_analysis(subset(qtl.full,rs %in% query[ind1,1]),
                           subset(df.full,rs %in% query[ind2,1]))
res[order(res$FDR)[1:10],]

# RINT vQTLs vs pure muQTLs
f <- paste0("/athena/elementolab/scratch/anm2868/open_targets/query/bmi.QTL.SNP_sig_list.query.pheWAS.txt")
qtl.full <- fread(f,data.table = F,stringsAsFactors = F)
ind1 <- which(sig_results$Rint.vQTL==1)
f <- paste0("/athena/elementolab/scratch/anm2868/open_targets/query/bmi.QTL.SNP_sig_list.query.pheWAS.txt")
df.full <- fread(f,data.table = F,stringsAsFactors = F)
ind2 <- which(sig_results$Rint.vQTL==0 & sig_results$Mean.QTL==1)
res <- pheWAS_enrichment_analysis(subset(qtl.full,rs %in% query[ind1,1]),
                                  subset(df.full,rs %in% query[ind2,1]))
res[order(res$FDR)[1:10],]

# dQTLs vs pure muQTLs
f <- paste0("/athena/elementolab/scratch/anm2868/open_targets/query/bmi.QTL.SNP_sig_list.query.pheWAS.txt")
qtl.full <- fread(f,data.table = F,stringsAsFactors = F)
ind1 <- which(sig_results$dQTL==1)
f <- paste0("/athena/elementolab/scratch/anm2868/open_targets/query/bmi.QTL.SNP_sig_list.query.pheWAS.txt")
df.full <- fread(f,data.table = F,stringsAsFactors = F)
ind2 <- which(sig_results$dQTL==0 & sig_results$Mean.QTL==1)
res <- pheWAS_enrichment_analysis(subset(qtl.full,rs %in% query[ind1,1]),
                                  subset(df.full,rs %in% query[ind2,1]))
res[order(res$FDR)[1:10],]

# QTLs vs matched SNPs
f <- paste0("/athena/elementolab/scratch/anm2868/open_targets/query/bmi.QTL.SNP_sig_list.query.pheWAS.txt")
qtl.full <- fread(f,data.table = F,stringsAsFactors = F)
ind1 <- which(sig_results$dQTL==1)
f <- paste0("/athena/elementolab/scratch/anm2868/open_targets/query/bmi.QTL.matched_snp.query.pheWAS.txt")
df.full <- fread(f,data.table = F,stringsAsFactors = F)
res <- pheWAS_enrichment_analysis(qtl.full,
                                  df.full)
res[order(res$FDR)[1:10],]


# used to browse a unique category/traits file
# x <- unique(df[,c('category','trait')])
# x <- x[order(x$category,x$trait),]
# fwrite(x,"/athena/elementolab/scratch/anm2868/open_targets/query/trait_category.txt",quote = F,row.names = F,col.names = T,na='NA',sep = '\t')
















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
