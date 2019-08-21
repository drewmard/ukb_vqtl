library(data.table)

pheno <- fread('/home/kulmsc/athena/ukbiobank/phenotypes/ukb26867.csv.gz',data.table=F,stringsAsFactors = F)
df.disease <- fread('/athena/elementolab/scratch/anm2868/vQTL/UKB/blood_disease.indiv_id.txt',data.table = F,stringsAsFactors = F,header = T)

# phenotypes
pheno3 <- pheno[,c('eid','30120-0.0','30130-0.0','30140-0.0','30200-0.0','30000-0.0')]
pheno3 <- subset(pheno3, !(eid %in% df.disease$eid))  # remove indiv w/ disease
PHENOTYPE_NAMES <- c('lymphocyte.count','monocyte.count','neutrophil.count','neutrophil.percentage','wbc.leukocyte.count')
colnames(pheno3)[2:ncol(pheno3)] <- PHENOTYPE_NAMES

fwrite(pheno3,'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/phenotypes.txt',na='NA',row.names = F,col.names = T,quote = F,sep = '\t')


