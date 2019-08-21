library(data.table)

set.seed(99)
fam2 <- fread('/home/kulmsc/athena/ukbiobank/calls/ukbb.1.fam',data.table = F,stringsAsFactors = F)
i.20 <- sample(1:nrow(fam2),floor(0.2*nrow(fam2)),replace = F)
fwrite(data.frame(row=i.20,IID=fam2[,2][i.20]),'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/phenotypes_blood.indiv_id.txt',sep='\t',quote = F,col.names = T,row.names = F)