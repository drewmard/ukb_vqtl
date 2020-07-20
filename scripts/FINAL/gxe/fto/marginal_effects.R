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

# env factors
val=1
mod <- list()
for (val in 1:4) {
  envir_name='Alcohol_intake_frequency'
  mod[[val]]<-lm(mod.formula,data = subset(df,Alcohol_intake_frequency==val))
}
mod[[5]]<-lm(mod.formula,data = subset(df,Alcohol_intake_frequency%in%c(5,6)))
gxe_effects <- as.data.frame(do.call(rbind,lapply(1:5,function(i){summary(mod[[i]])$coef['SNP',c(1,2)]})))
gxe_effects$E <- 'Alc'
gxe_effects$LOW <- gxe_effects$Estimate-1.96*gxe_effects$`Std. Error`
gxe_effects$HI <- gxe_effects$Estimate+1.96*gxe_effects$`Std. Error`
x <- c('Daily','3-4x/wk','1-2x/wk','1-3x/mo','Rarely/never')
gxe_effects$Val <- factor(x,levels = x)

mod <- list()
mod[[1]] <- lm(mod.formula,data = subset(df,age<=49))
mod[[2]] <- lm(mod.formula,data = subset(df,age>=50 & age<=59))
mod[[3]] <- lm(mod.formula,data = subset(df,age>=60))
gxe_effects.tmp <- as.data.frame(do.call(rbind,lapply(1:3,function(i){summary(mod[[i]])$coef['SNP',c(1,2)]})))
gxe_effects.tmp$E <- 'Age'
gxe_effects.tmp$LOW <- gxe_effects.tmp$Estimate-1.96*gxe_effects.tmp$`Std. Error`
gxe_effects.tmp$HI <- gxe_effects.tmp$Estimate+1.96*gxe_effects.tmp$`Std. Error`
gxe_effects.tmp$Val <- c('37-49','50-59','60-72')
gxe_effects <- rbind(gxe_effects,gxe_effects.tmp)

mod <- list()
mod[[1]] <- lm(mod.formula,data = subset(df,PA==1))
mod[[2]] <- lm(mod.formula,data = subset(df,PA==2))
mod[[3]] <- lm(mod.formula,data = subset(df,PA==3))
gxe_effects.tmp <- as.data.frame(do.call(rbind,lapply(1:3,function(i){summary(mod[[i]])$coef['SNP',c(1,2)]})))
gxe_effects.tmp$LOW <- gxe_effects.tmp$Estimate-1.96*gxe_effects.tmp$`Std. Error`
gxe_effects.tmp$HI <- gxe_effects.tmp$Estimate+1.96*gxe_effects.tmp$`Std. Error`
gxe_effects.tmp$E <- 'PA'
gxe_effects.tmp$Val <- c('Low','Mod','High')
gxe_effects <- rbind(gxe_effects,gxe_effects.tmp)

mod <- list()
mod[[1]] <- lm(mod.formula,data = subset(df,Smoking.E==0))
mod[[2]] <- lm(mod.formula,data = subset(df,Smoking.E==1))
gxe_effects.tmp <- as.data.frame(do.call(rbind,lapply(1:2,function(i){summary(mod[[i]])$coef['SNP',c(1,2)]})))
gxe_effects.tmp$LOW <- gxe_effects.tmp$Estimate-1.96*gxe_effects.tmp$`Std. Error`
gxe_effects.tmp$HI <- gxe_effects.tmp$Estimate+1.96*gxe_effects.tmp$`Std. Error`
gxe_effects.tmp$E <- 'Smok'
gxe_effects.tmp$Val <- c('No','Yes')
gxe_effects <- rbind(gxe_effects,gxe_effects.tmp)

mod <- list()
mod[[1]] <- lm(mod.formula,data = subset(df,sex==0))
mod[[2]] <- lm(mod.formula,data = subset(df,sex==1))
gxe_effects.tmp <- as.data.frame(do.call(rbind,lapply(1:2,function(i){summary(mod[[i]])$coef['SNP',c(1,2)]})))
gxe_effects.tmp$LOW <- gxe_effects.tmp$Estimate-1.96*gxe_effects.tmp$`Std. Error`
gxe_effects.tmp$HI <- gxe_effects.tmp$Estimate+1.96*gxe_effects.tmp$`Std. Error`
gxe_effects.tmp$E <- 'Sex'
gxe_effects.tmp$Val <- c('Female','Male')
gxe_effects <- rbind(gxe_effects,gxe_effects.tmp)

mod <- list()
mod[[1]] <- lm(mod.formula,data = subset(df,SB>=0 & SB<=2))
mod[[2]] <- lm(mod.formula,data = subset(df,SB>=3 & SB<=4))
mod[[3]] <- lm(mod.formula,data = subset(df,SB>=5 & SB<=6))
mod[[4]] <- lm(mod.formula,data = subset(df,SB>=7))
gxe_effects.tmp <- as.data.frame(do.call(rbind,lapply(1:4,function(i){summary(mod[[i]])$coef['SNP',c(1,2)]})))
gxe_effects.tmp$LOW <- gxe_effects.tmp$Estimate-1.96*gxe_effects.tmp$`Std. Error`
gxe_effects.tmp$HI <- gxe_effects.tmp$Estimate+1.96*gxe_effects.tmp$`Std. Error`
gxe_effects.tmp$E <- 'SB'
gxe_effects.tmp$Val <- c('0-2','3-4','5-6','7+')
gxe_effects <- rbind(gxe_effects,gxe_effects.tmp)

mod <- list()
x <- quantile(df$DIET_SCORE,probs=c(0,0.2,0.8,1),na.rm=T)
mod[[1]] <- lm(mod.formula,data = subset(df,DIET_SCORE <= x[2]))
mod[[2]] <- lm(mod.formula,data = subset(df,DIET_SCORE > x[2] & DIET_SCORE <= x[3]))
mod[[3]] <- lm(mod.formula,data = subset(df,DIET_SCORE > x[3]))
gxe_effects.tmp <- as.data.frame(do.call(rbind,lapply(1:3,function(i){summary(mod[[i]])$coef['SNP',c(1,2)]})))
gxe_effects.tmp$LOW <- gxe_effects.tmp$Estimate-1.96*gxe_effects.tmp$`Std. Error`
gxe_effects.tmp$HI <- gxe_effects.tmp$Estimate+1.96*gxe_effects.tmp$`Std. Error`
gxe_effects.tmp$E <- 'Diet'
gxe_effects.tmp$Val <- c('Low BMI','Med BMI','High BMI')
gxe_effects <- rbind(gxe_effects,gxe_effects.tmp)

f.out <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_results/FTO_marginal_effects.txt'
fwrite(gxe_effects,f.out,quote = F,sep = '\t',na='NA',row.names = F,col.names = T)

