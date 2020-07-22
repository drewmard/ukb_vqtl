# Run pheWAS enrichment analysis
pheWAS_enrichment_analysis <- function(qtl.full,df.full) { 
  qtl <- subset(qtl.full,pval <= PVAL_THRES)
  qtl <- subset(qtl,!duplicated(qtl[,c('rs',field)]))
  qtl <- subset(qtl,category %in% category_keep[,1])
  qtl <- subset(qtl,ntotal > SAMP_SIZE_MIN)
  n <- length(unique(qtl$rs))
  qtl <- table(qtl[,field])/n
  qtl <- data.frame(qtl)
  
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

# QTL data
f <- paste0("/athena/elementolab/scratch/anm2868/open_targets/query/bmi.QTL.SNP_sig_list.query.pheWAS.txt")
qtl.full <- fread(f,data.table = F,stringsAsFactors = F)
ind <- which(sig_results$Raw.vQTL==1)

qtl.full <- subset(qtl.full,rs %in% query[ind,1])

# Background data:
# matched
# f <- paste0("/athena/elementolab/scratch/anm2868/open_targets/query/bmi.QTL.matched_snp.query.pheWAS.txt")
# df.full <- fread(f,data.table = F,stringsAsFactors = F)

# muqtls
f <- paste0("/athena/elementolab/scratch/anm2868/open_targets/query/bmi.QTL.SNP_sig_list.query.pheWAS.txt")
df.full <- fread(f,data.table = F,stringsAsFactors = F)
ind <- which(sig_results$Raw.vQTL==0 & sig_results$Mean.QTL==1)
df.full <- subset(df.full,rs %in% query[ind,1])

# vqtl vs background
df.mg2 <- pheWAS_enrichment_analysis(qtl.full,df.full)
df.mg2[order(df.mg2$Observed,decreasing = T)[1:10],]
df.mg2[order(df.mg2$FDR,decreasing = F)[1:10],]
df.mg2.sub <- subset(df.mg2,FDR<0.1)
fwrite(df.mg2.sub,'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/phewas/phewas_muqtl_vs_vqtl.txt',quote = F,na = 'NA',sep = '\t',row.names = F,col.names = T)
fwrite(df.mg2,'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/phewas/phewas_muqtl_vs_vqtl.full.txt',quote = F,na = 'NA',sep = '\t',row.names = F,col.names = T)

# sampling 21 from background
rs.uniq=unique(df.full$rs)
i <- sample(1:length(rs.uniq),21,replace = F)
df.mg2 <- pheWAS_enrichment_analysis(subset(df.full,rs%in%rs.uniq[i]),
                                     subset(df.full,!(rs%in%rs.uniq[i])))
df.mg2.sub <- subset(df.mg2,FDR<0.1)
df.mg2.sub
fwrite(df.mg2.sub,'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/phewas/phewas_muqtl_vs_muqtl.txt',quote = F,na = 'NA',sep = '\t',row.names = F,col.names = T)
fwrite(df.mg2,'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/phewas/phewas_muqtl_vs_muqtl.full.txt',quote = F,na = 'NA',sep = '\t',row.names = F,col.names = T)


