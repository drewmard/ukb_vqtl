library(data.table)
library(stringr)

phenoName <- 'lymphocyte.count.rint.ALL'

# read in scores for all individuals
for (i in 1:22) {
  print(i)
  # f <- paste0('/home/kulmsc/athena/forAndrew/clumping/ukbb.',i,'.profile')
  # f <- paste0('/home/kulmsc/athena/forAndrew/noCorrection/ukbb.',i,'.profile')
  f <- paste0('/home/kulmsc/athena/forAndrew/clumping_impute/ukbb.',i,'.profile')
  df <- fread(f,data.table = F,stringsAsFactors = F)
  if (i == 1) {
    score.save <- df
  } else {
    score.save$SCORE <-score.save$SCORE+df$SCORE
  }
}

score.save <- fread('/athena/elementolab/scratch/kulmsc/forAndrew/clumping_impute/ukbb.all.profile',data.table = F,stringsAsFactors = F)
colnames(score.save)[1] <- 'IID'

# read in full data results for 20% UKB testing set
s <- '20'
f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/results/full_data_gxe.',s,'.txt')
df <- fread(f,data.table = F,stringsAsFactors = F)
x <- strsplit(phenoName,'\\.')[[1]]; phenoName2 <- paste0(c(x[-length(x)],'na'),collapse='.')

# read in vQTL results from 80% UKB discovery set
s <- '80'
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
dataf <- merge(dataf,score.save[,c('IID','SCORE')],by='IID')
colnames(dataf)[ncol(dataf)] <- 'Additive'
dataf$Joint <- dataf$Interaction+dataf$Additive

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




