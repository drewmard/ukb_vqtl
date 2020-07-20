library(data.table)

pheno = 'bmi'

# Matched snps
s='20';results.20 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.',s,'.diet_score.more_snp.QTL_matched_snp.txt'),data.table = F,stringsAsFactors = F)
s='80';results.80 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.',s,'.diet_score.more_snp.QTL_matched_snp.txt'),data.table = F,stringsAsFactors = F)
results.mg.diet <- merge(results.80,results.20,by=c('SNP','E'))

s='20';results.20 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.',s,'.ext.more_snp.QTL_matched_snp.txt'),data.table = F,stringsAsFactors = F)
s='80';results.80 <- fread(paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/bmi.GxE.',s,'.ext.more_snp.QTL_matched_snp.txt'),data.table = F,stringsAsFactors = F)
results.mg.all <- merge(results.80,results.20,by=c('SNP','E'))

results.mg <- rbind(results.mg.diet,results.mg.all)
# results.mg <- results.mg.all
g <- regexpr("_[^_]*$", results.mg$SNP)-1
results.mg$SNP <- substring(results.mg$SNP,1,g)

fwrite(results.mg,'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/sig_results/bmi.matched_QTL.GxE.txt',quote = F,na='NA',sep = '\t',row.names = F,col.names = T)

# qtl snps
f <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/sig_results/bmi.sig.GxE.txt'
sig_results <- fread(f,data.table = F,stringsAsFactors = F,header = T)
sig_results$FDR <- p.adjust(sig_results[,4],method='fdr')

# remove any snps (if any) that are also qtls
results.mg <- subset(results.mg,!(SNP %in% sig_results$SNP))

set.seed(031995)
THRES=0.1
nperm=10000
p.matched.vec <- rep(NA,nperm)
for (i in 1:nperm) {
  SNP_sampled <- sample(unique(results.mg$SNP),length(unique(sig_results$SNP)),replace = F)
  results.mg.sampled <- subset(results.mg,SNP %in% SNP_sampled)
  results.mg.sampled$FDR <- p.adjust(results.mg.sampled[,4],method = 'fdr')
  results.mg.sampled[order(results.mg.sampled$FDR)[1:5],]
  
  results.mg.sampled.sub <- subset(results.mg.sampled,FDR < THRES)
  p.matched <- nrow(results.mg.sampled.sub)
  p.matched.vec[i] <- p.matched
}
p.qtl <- nrow(subset(sig_results,FDR<THRES))/nrow(sig_results)
p.qtl/mean(p.matched.vec/nrow(results.mg.sampled))
test_res <- binom.test(nrow(subset(sig_results,FDR<=THRES)),nrow(sig_results),mean(p.matched.vec/nrow(results.mg.sampled)),alternative = 'greater')
test_res$p.value



