tis <- list('____')

library(data.table)
df.exp <- fread('/athena/elementolab/scratch/anm2868/GTEx/EXPRESSION/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct',data.table = F,stringsAsFactors = F)
GENE <- list('TMEM18');
EXPR.df <- df.exp[which(df.exp$Description %in% GENE),-1]
rownames(EXPR.df) <- EXPR.df$Description; EXPR.df <- EXPR.df[,-1]
EXPR.df <- as.data.frame(t(EXPR.df))
EXPR.df$SAMP <- rownames(EXPR.df); rownames(EXPR.df) <- NULL
paste.s <- function(x) {return(paste(x[1:2],collapse='-'))}
EXPR.df$SAMP2 <- sapply(strsplit(EXPR.df$SAMP,"-"),paste.s)

df.infil <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-T.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)

df.mg <- merge(EXPR.df[,c('SAMP',GENE[[1]])],df.infil,by.x='SAMP',by.y='Input Sample')
df.mg.sub <- subset(df.mg,SMTSD %in% c('Adipose - Subcutaneous','Adipose - Visceral (Omentum)'))
df.mg.sub <- subset(df.mg,SMTSD %in% c('Adipose - Subcutaneous'))
df.mg.sub <- subset(df.mg,SMTSD %in% c('Adipose - Visceral (Omentum)'))

t.test(subset(df.mg.sub,SMTSD=='Adipose - Subcutaneous')[,GENE[[1]]],
       subset(df.mg.sub,SMTSD=='Adipose - Visceral (Omentum)')[,GENE[[1]]])

cor.test(df.mg.sub[,GENE[[1]]],as.numeric(as.factor(df.mg.sub$AGE)))
cor.test(df.mg.sub[,GENE[[1]]],as.numeric(as.factor(df.mg.sub$AGE)))

tmp <- df.mg.sub[,c(GENE[[1]],'AGE')]
fwrite(tmp,'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/Age_GTEx.txt',quote = F,na = 'NA',sep='\t',row.names = F,col.names = T)

