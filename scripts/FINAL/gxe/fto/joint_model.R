library(data.table)
library(BEDMatrix)

# phenotype data
s <- '80'; pheno <- 'bmi'
f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/full_data.',pheno,'.GxE.',s,'.txt')
df <- fread(f,data.table = F,stringsAsFactors = F)
f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/diet_data.',pheno,'.GxE.',s,'.txt')
diet_score <- fread(f,data.table = F,stringsAsFactors = F)
df <- merge(df,diet_score,by='IID',all=TRUE)

# genetic data
f.geno <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.',pheno,'.merged_subset2')
geno <- BEDMatrix(f.geno)
geno_names <- unlist(lapply(strsplit(rownames(geno),'_'),function(x) {return(x[2])}))

# input
environmental_factors <- c(
  'DIET_SCORE',
  'age','Alcohol_intake_frequency',
  'PA','SB','sex','Smoking.E')
SNP_name='rs56094641'

# pre process
i <- grep(SNP_name,colnames(geno))
ind <- which(geno_names %in% df$IID)
df$SNP <- geno[ind,i]
phenoName <- paste0(pheno,'.na')
mod.formula <- formula(paste(phenoName,' ~ age+age2+genotyping.array+sex+age*sex+age2*sex+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
               PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+SNP'))
mod <- lm(mod.formula,data = df)
summary(mod)$coef['SNP',c(1,4)]

mod.formula <- formula(paste(phenoName,' ~ age+age2+genotyping.array+sex+age*sex+age2*sex+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
               PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+
                             SNP*Smoking.E+
                             SNP*age+
                             SNP*DIET_SCORE+
                             SNP*SB+
                             SNP*PA+
                             SNP*Alcohol_intake_frequency'))
mod <- lm(mod.formula,data = df)
summary(mod)

tmp <- as.data.frame(
  summary(mod)$coef[c('SNP:Smoking.E',
                      'age:SNP',
                      'SNP:DIET_SCORE',
                      'SNP:SB',
                      'SNP:PA',
                      'SNP:Alcohol_intake_frequency'),]
)
f <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/FTO_joint_model_interaction_effects.txt'
fwrite(tmp,f,sep = '\t',quote = F,na = 'NA',row.names = T,col.names = T)
