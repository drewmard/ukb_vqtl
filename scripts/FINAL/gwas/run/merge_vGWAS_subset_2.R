library(data.table)
args = commandArgs(trailingOnly=TRUE)
phenotype=args[1]
# phenotype='lymphocyte.count.rint.ALL'

f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_subset/','ukbb.',phenotype,'.vGWAS.old.txt')
df.res.save <- fread(f.out,data.table = F,stringsAsFactors = F)

g <- regexpr("_[^_]*$", df.res.save$SNP)-1
df.res.save$rs <- substring(df.res.save$SNP,1,g)

colnames(df.res.save) <- c('BETA','SE','T','P','SNP','rs')

f.out <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_subset/','ukbb.',phenotype,'.vGWAS.txt')
print(paste0('Writing: ',f.out))
fwrite(df.res.save,f.out,sep='\t',quote=F,col.names = T,row.names = F,na="NA")

