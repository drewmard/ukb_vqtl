library(data.table)

s <- '80'
f.res <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/results/ukbb.lymphocyte.count.rint.ALL.results.txt')
results <- fread(f.res,data.table = F,stringsAsFactors = F)
res.sub <- subset(results,MAF > 0.01)
res.sub <- res.sub[,c('SNP','CHR','BP','REF','ALT','BETA.x','P.x')]
colnames(res.sub)[(ncol(res.sub)-1):ncol(res.sub)] <- c('BETA','P')
fwrite(res.sub,'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/downstream/ukbb.lymphocyte.count.rint.ALL.results.PGS.txt',col.names = T,row.names = F,sep = '\t',na = 'NA',quote = F)


