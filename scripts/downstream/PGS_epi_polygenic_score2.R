library(data.table)
library(stringr)
workdir='/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/imputed/results'

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
s <- '20'
f <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/results/full_data_gxe.',s,'.txt')
df <- fread(f,data.table = F,stringsAsFactors = F)
x <- strsplit(phenoName,'\\.')[[1]]; phenoName2 <- paste0(c(x[-length(x)],'na'),collapse='.')

# read in vQTL results from 80% UKB discovery set
# s <- '80'
f.res <- paste0('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG/ukbb.lymphocyte.count.rint.ALL.ALL.sub.GxG.epi.qt')
results <- fread(f.res,data.table = F,stringsAsFactors = F)
res.sub <- subset(results,P < 0.05/nrow(results))
res.sub <- res.sub[order(res.sub$P),]
res.sub <- subset(res.sub,SNP1!='rs146125856' & SNP2!='rs146125856')
# identify set of genetic interactions for using an "interaction" score.
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
# res.sub.save is the interaction set

f <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG/ukbb.lymphocyte.count.rint.ALL.ALL.sub.GxG.epi.correlation_purge2.txt'
fwrite(res.sub.save,f,col.names = T,row.names = F,quote = F,na='NA',sep = '\t')


# library(data.table)
# library(stringr)
# f <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG/ukbb.lymphocyte.count.rint.ALL.ALL.sub.GxG.epi.correlation_purge2.txt'
res.sub.save <- fread(f,data.table = F,stringsAsFactors = F)
# calculate interaction score:
score.interaction <- 0
cor.vec <- c()
for (i in 1:nrow(res.sub.save)) {
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
dataf <- merge(dataf,df[,c('IID','age','sex','PC1','PC2','PC3','PC4','genotyping.array','lymphocyte.count.na')],by='IID')
dataf$Joint <- dataf$Interaction+dataf$Additive
fwrite(dataf,'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/downstream/ukbb.lymphocyte.count.rint.ALL.results.PGS.data.txt',quote = F,sep = '\t',na = 'NA',row.names = F,col.names = T)

library(data.table)
f <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/downstream/ukbb.lymphocyte.count.rint.ALL.results.PGS.data.txt'
dataf<- fread(f,data.table = F,stringsAsFactors = F)
# read2=TRUE
# if (read2) {
#   pheno <- fread('/home/kulmsc/athena/ukbiobank/phenotypes/ukb26867.csv.gz',data.table=F,stringsAsFactors = F)
# }
# pheno2 <- pheno[,c('eid','21000-0.0','30120-0.0','21003-0.0','22001-0.0',
#                    '22009-0.1','22009-0.2','22009-0.3','22009-0.4','21001-0.0')]
# colnames(pheno2) <- c('IID','Ethnicity','Lymphocyte.Count','Age','Sex',
#                       'PC1','PC2','PC3','PC4','BMI')
# dataf2 <- merge(dataf,pheno2,by='IID')


# if reading in data
# f <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/downstream/ukbb.lymphocyte.count.rint.ALL.results.PGS.data.txt'
# dataf<- fread(f,data.table = F,stringsAsFactors = F)

# dataf <- subset(dataf,IID %in% indiv)


res=cor.test(dataf$Additive,dataf$Phenotype,use='p'); res[c('estimate','p.value')]
res
res=cor.test(dataf$Interaction,dataf$Phenotype,use='p'); res[c('estimate','p.value')]
res

cor(dataf$Additive,dataf$Interaction)

res=cor.test(dataf$Additive,dataf$lymphocyte.count.na,use='p'); res[c('estimate','p.value')]
res
mod.formula <- paste0('lymphocyte.count.na','~','Additive+Interaction+age+sex+PC1+PC2+PC3+PC4+genotyping.array')
mod.formula <- formula(mod.formula)
mod.int_add <- lm(mod.formula,data=dataf)
mod.int_add <- lm(mod.formula,data=dataf[i,])
# summary(mod.epi)
mod.formula <- paste0('lymphocyte.count.na','~','Additive','+age+sex+PC1+PC2+PC3+PC4+genotyping.array')
mod.formula <- formula(mod.formula)
mod.add <- lm(mod.formula,data=dataf)
mod.formula <- paste0('lymphocyte.count.na','~','age+sex+PC1+PC2+PC3+PC4+genotyping.array')
mod.formula <- formula(mod.formula)
mod.no_genet <- lm(mod.formula,data=dataf)
summary(mod.int_add)$adj.r.squared/summary(mod.add)$adj.r.squared
summary(mod.int_add)$adj.r.squared/summary(mod.no_genet)$adj.r.squared
summary(mod.add)$adj.r.squared/summary(mod.no_genet)$adj.r.squared
summary(mod.int_add)$adj.r.squared; summary(mod.add)$adj.r.squared;summary(mod.no_genet)$adj.r.squared

i <- sample(1:nrow(dataf),floor(nrow(dataf)/2),replace = F) # to identify train/test
mod.formula <- paste0('lymphocyte.count.na','~','Additive+age+sex+PC1+PC2+PC3+PC4+genotyping.array')
mod.formula <- formula(mod.formula)
mod.test <- lm(mod.formula,data=dataf[i,])
summary(mod.test)
afss=anova(mod.test)$"Sum Sq"; x=afss/sum(afss)*100; x=x[1];x
pred <- predict(mod.test,dataf[-i,])
cor.add <- cor(pred,dataf[-i,'Phenotype'],use='p')
mod.formula <- paste0('lymphocyte.count.na','~','Additive+Interaction+age+sex+PC1+PC2+PC3+PC4+genotyping.array')
mod.formula <- formula(mod.formula)
mod.test <- lm(mod.formula,data=dataf[i,])
summary(mod.test)
afss=anova(mod.test)$"Sum Sq"; x=afss/sum(afss)*100; x=x[1]+x[2];x
pred <- predict(mod.test,dataf[-i,])
cor.add_int <- cor(pred,dataf[-i,'Phenotype'],use='p')
mod.formula <- paste0('lymphocyte.count.na','~','age+sex+PC1+PC2+PC3+PC4+genotyping.array')
mod.formula <- formula(mod.formula)
mod.test <- lm(mod.formula,data=dataf[i,])
summary(mod.test)
pred <- predict(mod.test,dataf[-i,])
cor.no_genet <- cor(pred,dataf[-i,'Phenotype'],use='p')
mod.formula <- paste0('lymphocyte.count.na','~','Interaction+age+sex+PC1+PC2+PC3+PC4+genotyping.array')
mod.formula <- formula(mod.formula)
mod.test <- lm(mod.formula,data=dataf[i,])
summary(mod.test)
pred <- predict(mod.test,dataf[-i,])
cor.int <- cor(pred,dataf[-i,'Phenotype'],use='p')
cor.int
cor.no_genet
cor.add
cor.add_int






mod.formula <- paste0('lymphocyte.count.na','~','Additive+Interaction')
mod.formula <- formula(mod.formula)
mod.int_add <- lm(mod.formula,data=dataf)
# summary(mod.epi)
mod.formula <- paste0('lymphocyte.count.na','~','Additive')
mod.formula <- formula(mod.formula)
mod.add <- lm(mod.formula,data=dataf)
summary(mod.int_add)$adj.r.squared/summary(mod.add)$adj.r.squared


# res=cor.test(dataf2$Additive,dataf2$Lymphocyte.Count,use='p'); res[c('estimate','p.value')]
# res
# mod.formula <- paste0('Lymphocyte.Count','~','Additive+Interaction+Age+Sex+PC1+PC2+PC3+PC4')
# mod.formula <- formula(mod.formula)
# mod.int_add <- lm(mod.formula,data=dataf2)
# # summary(mod.epi)
# mod.formula <- paste0('Lymphocyte.Count','~','Additive','+Age+Sex+PC1+PC2+PC3+PC4')
# mod.formula <- formula(mod.formula)
# mod.add <- lm(mod.formula,data=dataf2)
# summary(mod.int_add)$adj.r.squared/summary(mod.add)$adj.r.squared

mod.formula <- paste0('Phenotype','~','Additive+Interaction')
mod.formula <- formula(mod.formula)
mod.int_add <- lm(mod.formula,data=dataf)
# summary(mod.epi)
mod.formula <- paste0('Phenotype','~','Additive')
mod.formula <- formula(mod.formula)
mod.add <- lm(mod.formula,data=dataf)
mod.formula <- paste0('Phenotype','~','Interaction')
mod.formula <- formula(mod.formula)
mod.int <- lm(mod.formula,data=dataf)
mod.formula <- paste0('Phenotype','~','Joint')
mod.formula <- formula(mod.formula)
mod.joint <- lm(mod.formula,data=dataf)
summary(mod.int_add)$adj.r.squared/summary(mod.add)$adj.r.squared
summary(mod.int_add)
summary(mod.int)
summary(mod.add)
summary(mod.joint)

# score prediction correlation
f <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/downstream/ukbb.lymphocyte.count.rint.ALL.results.PGS.data.txt'
dataf<- fread(f,data.table = F,stringsAsFactors = F)
i <- sample(1:nrow(dataf),floor(nrow(dataf)/2),replace = F) # to identify train/test
i.test <- (1:nrow(dataf))[-i]

mod.formula <- paste0('Phenotype','~','Additive')
mod.formula <- formula(mod.formula)
mod.add <- lm(mod.formula,data=dataf[i,])

mod.formula <- paste0('Phenotype','~','Additive+Interaction')
mod.formula <- formula(mod.formula)
mod.add_int <- lm(mod.formula,data=dataf[i,])

nsim <- 100
diff <- rep(NA,nsim)
r.add.vec <- rep(NA,nsim)
r.int.vec <- rep(NA,nsim)
for (j in 1:nsim) {
  i.test.bs <- sample(i.test,length(i.test),replace = T)
  pred <- predict(mod.add,dataf[i.test.bs,])
  r.add <- cor(pred,dataf[i.test.bs,'Phenotype'],use='p'); r.add
  
  pred <- predict(mod.add_int,dataf[i.test.bs,])
  r.int <- cor(pred,dataf[i.test.bs,'Phenotype'],use='p'); r.int
  
  r.add.vec[j] <- r.add.vec
  r.int.vec[j] <- r.int.vec
  diff[j] <- r.int - r.add
  
}





##########

set.seed(03191995)
nsim <- 1000
diff <- rep(NA,nsim)
r.add.vec <- rep(NA,nsim)
r.int.vec <- rep(NA,nsim)
r.add_int.vec <- rep(NA,nsim)
r.baseline.vec <- rep(NA,nsim)
for (j in 1:nsim) {
  if (j %% 10 == 0) {print(j)}
  i <- sample(1:nrow(dataf),floor(nrow(dataf)/2),replace = F) # to identify train/test
  i.test.bs <- (1:nrow(dataf))[-i]
  
  mod.formula <- paste0('Phenotype','~','Additive+age+sex+PC1+PC2+PC3+PC4+genotyping.array')
  mod.formula <- formula(mod.formula)
  mod.add <- lm(mod.formula,data=dataf[i,])
  
  mod.formula <- paste0('Phenotype','~','Additive+Interaction+age+sex+PC1+PC2+PC3+PC4+genotyping.array')
  mod.formula <- formula(mod.formula)
  mod.add_int <- lm(mod.formula,data=dataf[i,])

  mod.formula <- paste0('Phenotype','~','Interaction+age+sex+PC1+PC2+PC3+PC4+genotyping.array')
  mod.formula <- formula(mod.formula)
  mod.int <- lm(mod.formula,data=dataf[i,])
  
  mod.formula <- paste0('Phenotype','~','age+sex+PC1+PC2+PC3+PC4+genotyping.array')
  mod.formula <- formula(mod.formula)
  mod.base <- lm(mod.formula,data=dataf[i,])
  
  pred <- predict(mod.add,dataf[i.test.bs,])
  r.add <- cor(pred,dataf[i.test.bs,'Phenotype'],use='p'); r.add
  
  pred <- predict(mod.add_int,dataf[i.test.bs,])
  r.add_int <- cor(pred,dataf[i.test.bs,'Phenotype'],use='p'); r.int
  
  pred <- predict(mod.int,dataf[i.test.bs,])
  r.int <- cor(pred,dataf[i.test.bs,'Phenotype'],use='p'); r.int
  
  pred <- predict(mod.base,dataf[i.test.bs,])
  r.baseline <- cor(pred,dataf[i.test.bs,'Phenotype'],use='p'); r.int
  
  r.add_int.vec[j] <- r.add_int
  r.int.vec[j] <- r.int
  r.add.vec[j] <- r.add
  r.baseline.vec[j] <- r.baseline
  
  diff[j] <- r.int - r.add
  
}
quantile(r.add_int.vec,c(0.025,0.975))
quantile(r.add.vec,c(0.025,0.975))
quantile(r.int.vec,c(0.025,0.975))
quantile(r.baseline.vec,c(0.025,0.975))
mean(r.add_int.vec)
mean(r.add.vec)
mean(r.int.vec)
mean(r.baseline.vec)

mean(r.int.vec)/mean(r.add.vec)


table(dataf$genotyping.array)


i <- which(dataf$genotyping.array=='UKBB')
mod.formula <- paste0('Phenotype','~','Additive+age+sex+PC1+PC2+PC3+PC4')
mod.formula <- formula(mod.formula)
mod.add <- lm(mod.formula,data=dataf[i,])

mod.formula <- paste0('Phenotype','~','Additive+Interaction+age+sex+PC1+PC2+PC3+PC4')
mod.formula <- formula(mod.formula)
mod.add_int <- lm(mod.formula,data=dataf[i,])

mod.formula <- paste0('Phenotype','~','Interaction+age+sex+PC1+PC2+PC3+PC4')
mod.formula <- formula(mod.formula)
mod.int <- lm(mod.formula,data=dataf[i,])

mod.formula <- paste0('Phenotype','~','age+sex+PC1+PC2+PC3+PC4')
mod.formula <- formula(mod.formula)
mod.base <- lm(mod.formula,data=dataf[i,])

pred <- predict(mod.add,dataf[i.test.bs,])
r.add <- cor(pred,dataf[i.test.bs,'Phenotype'],use='p'); r.add

pred <- predict(mod.add_int,dataf[i.test.bs,])
r.add_int <- cor(pred,dataf[i.test.bs,'Phenotype'],use='p'); r.add_int

pred <- predict(mod.int,dataf[i.test.bs,])
r.int <- cor(pred,dataf[i.test.bs,'Phenotype'],use='p'); r.int

pred <- predict(mod.base,dataf[i.test.bs,])
r.baseline <- cor(pred,dataf[i.test.bs,'Phenotype'],use='p'); r.baseline
