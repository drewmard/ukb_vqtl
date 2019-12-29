# library(data.table)
# library(stringr)
# 
# phenoName <- 'lymphocyte.count.rint.ALL'
# 
# # read in scores for all individuals
# for (i in 1:22) {
#   print(i)
#   # f <- paste0('/home/kulmsc/athena/forAndrew/clumping/ukbb.',i,'.profile')
#   # f <- paste0('/home/kulmsc/athena/forAndrew/noCorrection/ukbb.',i,'.profile')
#   f <- paste0('/home/kulmsc/athena/forAndrew/clumping_impute/ukbb.',i,'.profile')
#   df <- fread(f,data.table = F,stringsAsFactors = F)
#   if (i == 1) {
#     score.save <- df
#   } else {
#     score.save$SCORE <-score.save$SCORE+df$SCORE
#   }
# }
# colnames(score.save)[1] <- 'IID'

library(data.table)
library(stringr)
workdir='/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/imputed/results'

# phenoName <- 'lymphocyte.count.rint.ALL'
# # read in scores for all individuals
# for (i in 1:22) {
#   print(i)
#   f <- paste0(workdir,'/ukbb.impute.',phenoName,'.muGWAS.chr',i,'.p_0.005.r_0.8.kb_250.clumped.PGS.profile')
#   df <- fread(f,data.table = F,stringsAsFactors = F)
#   if (i == 1) {
#     score.save <- df
#   } else {
#     score.save$SCORESUM <-score.save$SCORESUM+df$SCORESUM
#   }
# }

phenoName <- 'lymphocyte.count.rint.ALL'
# read in scores for all individuals
for (i in 1:22) {
  print(i)
  f <- paste0(workdir,'/ukbb.impute.',phenoName,'.muGWAS.chr',i,'.p_0.005.r_0.8.kb_250.clumped.PGS.profile')
  df <- fread(f,data.table = F,stringsAsFactors = F)
  df$SCORESUM <- -1*df$SCORESUM
  if (i == 1) {
    colnames(df)[6] <- paste0('SCORESUM.',i)
    score.save <- df[,c('IID',paste0('SCORESUM.',i))]
  } else {
    colnames(df)[6] <- paste0('SCORESUM.',i)
    score.save <- merge(score.save,df[,c('IID',paste0('SCORESUM.',i))],by='IID')
  }
}

score.save$Additive <- apply(score.save[,grepl('SCORESUM',colnames(score.save))],1,sum)

# read in full data results for 20% UKB testing set
# s <- '80'
# f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/results/full_data_gxe.',s,'.txt')
# df <- fread(f,data.table = F,stringsAsFactors = F)
s <- '20'
f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/results/full_data_gxe.',s,'.txt')
df2 <- fread(f,data.table = F,stringsAsFactors = F)
# df <- as.data.frame(rbind(df,df2))
x <- strsplit(phenoName,'\\.')[[1]]; phenoName2 <- paste0(c(x[-length(x)],'na'),collapse='.')

# read in vQTL results from 80% UKB discovery set
# s <- '80'
f.res <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG/ukbb.lymphocyte.count.rint.ALL.ALL.sub.GxG.epi.qt')
results <- fread(f.res,data.table = F,stringsAsFactors = F)
res.sub <- subset(results,P < 0.05/nrow(results))
res.sub <- res.sub[order(res.sub$P),]

# identify set of genetic interactions for using an "interaction" score.
mode <- 'one_interaction_only'
mode <- 'correlation_purge'; calculate_correlation_matrix <- FALSE
mode <- 'correlation_purge2';

if (mode=='correlation_purge2') {
  
  x <- matrix(NA,nrow(df),nrow(res.sub))
  colnames(x) <- paste0('Int',1:nrow(res.sub))
  for (j in 1:nrow(res.sub)) {
    SNP1 <- res.sub$SNP1[j]
    SNP2 <- res.sub$SNP2[j]
    SNP1 <- str_replace(SNP2,'-','.')
    SNP2 <- str_replace(SNP2,'-','.')
    x[,j] <- df[,SNP1]*df[,SNP2]
  }
  
  too_correlated <- c()
  keepers <- c()
  for (i in 1:nrow(res.sub)) {
    print(i)
    if (i == 1) {
      res.sub.save <- res.sub[i,]
      cor.mat <- cor(as.numeric(x[,i]),as.matrix(x),use='p')
      ind <- as.numeric(which(cor.mat[1,] > 0.8))
      keepers <- c(keepers,i)
      too_correlated <- c(too_correlated,ind)
    } else if (!(i %in% too_correlated)) {
      res.sub.save <- rbind(res.sub.save,res.sub[i,])
      cor.mat <- cor(as.numeric(x[,i]),as.matrix(x),use='p')
      ind <- as.numeric(which(cor.mat[1,] > 0.8))
      keepers <- c(keepers,i)
      too_correlated <- c(too_correlated,ind)
    }
  }
}
if (mode=='correlation_purge') {
  
  if (calculate_correlation_matrix) {
    x <- matrix(NA,nrow(df),nrow(res.sub))
    colnames(x) <- paste0('Int',1:nrow(res.sub))
    for (j in 1:nrow(res.sub)) {
      SNP1 <- res.sub$SNP1[j]
      SNP2 <- res.sub$SNP2[j]
      SNP1 <- str_replace(SNP2,'-','.')
      SNP2 <- str_replace(SNP2,'-','.')
      x[,j] <- df[,SNP1]*df[,SNP2]
    }
    cor.mat <- cor(x,use='p')
  }
  
  too_correlated <- c()
  keepers <- c()
  for (i in 1:nrow(res.sub)) {
    SNP1 <- res.sub$SNP1[i]
    SNP2 <- res.sub$SNP2[i]
    if (i == 1) {
      res.sub.save <- res.sub[i,]
      ind <- as.numeric(which(cor.mat[i,] > 0.8))
      keepers <- c(keepers,i)
      too_correlated <- c(too_correlated,ind)
    } else if (!(i %in% too_correlated)) {
      res.sub.save <- rbind(res.sub.save,res.sub[i,])
      ind <- as.numeric(which(cor.mat[i,] > 0.8))
      keepers <- c(keepers,i)
      too_correlated <- c(too_correlated,ind)
    }
  }
}
if (mode=='one_interaction_only') {
  SNP_list <- c()
  for (i in 1:nrow(res.sub)) {
    SNP1 <- res.sub$SNP1[i]
    SNP2 <- res.sub$SNP2[i]
    if (i == 1) {
      res.sub.save <- res.sub[i,]
      SNP_list <- c(SNP_list,SNP1,SNP2)
    } else if (!(SNP1 %in% SNP_list) & !(SNP2 %in% SNP_list)) {
      res.sub.save <- rbind(res.sub.save,res.sub[i,])
      SNP_list <- c(SNP_list,SNP1,SNP2)
    }
  }
}
# res.sub.save is the interaction set

fwrite(res.sub.save,'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG/ukbb.lymphocyte.count.rint.ALL.ALL.sub.GxG.epi.correlation_purge2.txt')

# calculate interaction score:
score.interaction <- 0
cor.vec <- c()
# res.sub.save <- res.sub
for (i in 1:nrow(res.sub.save)) {
  # for (i in 1:5) {
  BETA.INT <- res.sub.save$BETA_INT[i]
  SNP1 <- str_replace_all(res.sub.save$SNP1[i],'-','.')
  SNP2 <- str_replace_all(res.sub.save$SNP2[i],'-','.')
  score.tmp <- BETA.INT*df[,SNP1]*df[,SNP2]
  score.tmp[is.na(score.tmp)] <- mean(score.tmp,na.rm=T)
  score.interaction <- score.interaction + score.tmp
  cor.vec <- c(cor.vec,cor(score.interaction,df[,phenoName2],use='p'))
}

dataf <- data.frame(IID=df$IID,Interaction=score.interaction,Phenotype=df[,phenoName2])

dataf <- merge(dataf,score.save,by='IID') # andrew's sep chr scores
# dataf <- merge(dataf,score.save[,c('IID','SCORE')],by='IID') # scott's scores
# dataf <- merge(dataf,score.save[,c('IID','SCORESUM')],by='IID') # andrew's scores

# colnames(dataf)[ncol(dataf)] <- 'Additive'
dataf$Joint <- dataf$Interaction+dataf$Additive

for (i in 1:22) {
  dataf[,paste0('SCORESUM.',i,'.scaled')] <- scale(dataf[,paste0('SCORESUM.',i)])
}
fwrite(dataf,'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/downstream/ukbb.lymphocyte.count.rint.ALL.results.PGS.data.txt',quote = F,sep = '\t',na = 'NA',row.names = F,col.names = T)

#############
# if reading in data
f <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/downstream/ukbb.lymphocyte.count.rint.ALL.results.PGS.data.txt'
dataf<- fread(f,data.table = F,stringsAsFactors = F)
read2=TRUE
if (read2) {
  pheno <- fread('/home/kulmsc/athena/ukbiobank/phenotypes/ukb26867.csv.gz',data.table=F,stringsAsFactors = F)
}
pheno2 <- pheno[,c('eid','21000-0.0','30120-0.0','21003-0.0','22001-0.0',
                   '22009-0.1','22009-0.2','22009-0.3','22009-0.4','21001-0.0')]
colnames(pheno2) <- c('IID','Ethnicity','Lymphocyte.Count','Age','Sex',
                      'PC1','PC2','PC3','PC4','BMI')
###################

res=cor.test(dataf$SCORESUM.3,dataf$Phenotype,use='p'); res[c('estimate','p.value')]


i <- sample(1:nrow(dataf),floor(nrow(dataf)/2),replace = F) # to identify train/test

# 22 diff chromosome scores
mod.formula <- paste0('Phenotype','~','SCORESUM.1+SCORESUM.2+SCORESUM.3+SCORESUM.4+
                      SCORESUM.5+SCORESUM.6+SCORESUM.7+SCORESUM.8+
                      SCORESUM.9+SCORESUM.10+SCORESUM.11+SCORESUM.12+
                      SCORESUM.13+SCORESUM.14+SCORESUM.15+SCORESUM.16+
                      SCORESUM.17+SCORESUM.18+SCORESUM.19+SCORESUM.20+
                      SCORESUM.21+SCORESUM.22')
mod.formula <- formula(mod.formula)
mod.test <- lm(mod.formula,data=dataf[i,])
summary(mod.test)
pred <- predict(mod.test,dataf[-i,])

# just one score over all chromosomes
mod.formula <- paste0('Phenotype','~','Additive.scaled')
mod.formula <- formula(mod.formula)
mod.test <- lm(mod.formula,data=dataf[i,])
summary(mod.test)
pred2 <- predict(mod.test,dataf[-i,])

cor(pred,dataf[-i,'Phenotype'],use='p')
cor(pred2,dataf[-i,'Phenotype'],use='p')

# mod.test <- lm(mod.formula,data=dataf)
# summary(mod.test)



dataf$Interaction.scaled <- scale(dataf$Interaction)
dataf$Additive.scaled <- scale(dataf$Additive)
dataf$Joint.scaled <- scale(dataf$Joint)
res=cor.test(dataf$Additive.scaled,dataf$Interaction.scaled); res[c('estimate','p.value')]
res=cor.test(dataf$Additive.scaled,dataf$Phenotype,use='p'); res[c('estimate','p.value')]
res=cor.test(dataf$Interaction.scaled,dataf$Phenotype,use='p'); res[c('estimate','p.value')]
res=cor.test(dataf$Additive.scaled,dataf$Joint.scaled,use='p'); res[c('estimate','p.value')]

mod.formula <- paste0('Phenotype','~','Additive.scaled+Interaction.scaled')
mod.formula <- formula(mod.formula)
mod.int_add <- lm(mod.formula,data=dataf)
# summary(mod.epi)
mod.formula <- paste0('Phenotype','~','Additive.scaled')
mod.formula <- formula(mod.formula)
mod.add <- lm(mod.formula,data=dataf)
mod.formula <- paste0('Phenotype','~','Interaction.scaled')
mod.formula <- formula(mod.formula)
mod.int <- lm(mod.formula,data=dataf)
mod.formula <- paste0('Phenotype','~','Joint.scaled')
mod.formula <- formula(mod.formula)
mod.joint <- lm(mod.formula,data=dataf)
summary(mod.int_add)$adj.r.squared/summary(mod.add)$adj.r.squared
summary(mod.int_add)
summary(mod.int)
summary(mod.add)
summary(mod.joint)

afss=anova(mod.int_add)$"Sum Sq"; x=afss/sum(afss)*100; x=x[1]+x[2];x
afss=anova(mod.add)$"Sum Sq"; y=afss/sum(afss)*100; y=y[1];y
afss=anova(mod.int)$"Sum Sq"; z=afss/sum(afss)*100; z=z[1];z
afss=anova(mod.joint)$"Sum Sq"; w=afss/sum(afss)*100; w=w[1];w

x/y
w/y

# mod.formula <- paste0(phenoName2,'~','PGS.epi+SCORE+age+as.factor(menopause2) ')
# mod.formula <- formula(mod.formula)
# mod <- lm(mod.formula,data=df.save)
# summary(mod)




