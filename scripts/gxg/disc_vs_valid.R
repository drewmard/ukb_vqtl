library(data.table)
df1 <- fread('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.bmi.merged_subset2.GxG.FULL.txt',data.table = F,stringsAsFactors = F)
df2 <- fread('/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/sig_results/bmi.sig.txt',data.table = F,stringsAsFactors = F)
df2.sub <- df2[,c('rs','P.MEAN','P.VAR.RAW','P.VAR.RINT','dispersion_pval',
       'BETA.MEAN','BETA.VAR.RAW','BETA.VAR.RINT','dispersion',
       'Mean.QTL','Raw.vQTL','Rint.vQTL','dQTL')]
df2.sub.tmp <- df2.sub; colnames(df2.sub.tmp) <- paste0(colnames(df2.sub.tmp),'.1')
df.mg <- merge(df1,df2.sub.tmp,by.x='SNP1',by.y='rs.1')
df2.sub.tmp <- df2.sub; colnames(df2.sub.tmp) <- paste0(colnames(df2.sub.tmp),'.2')
df.mg <- merge(df.mg,df2.sub.tmp,by.x='SNP2',by.y='rs.2')

# df.mg$Mean.QTL <- df.mg$Mean.QTL.1 + df.mg$Mean.QTL.2
# df.mg$Raw.vQTL <- df.mg$Raw.vQTL.1 + df.mg$Raw.vQTL.2
# df.mg$Rint.vQTL <- df.mg$Rint.vQTL.1 + df.mg$Rint.vQTL.2
# df.mg$dQTL <- df.mg$dQTL.1 + df.mg$dQTL.2
# summary(lm(BETA_INT.80~Mean.QTL*Raw.vQTL*Rint.vQTL*dQTL,data=df.mg))

i <- which(df.mg$P.80 < 0.001)
i <- which(df.mg$SNP1=='rs56094641' | df.mg$SNP2=='rs56094641')

df <- list()
j=1

i <- which(df.mg$P.MEAN.1 < 5e-8 | df.mg$P.MEAN.2 < 5e-8)
res <- cor.test(df.mg$BETA_INT.20[i],df.mg$BETA_INT.80[i])
df[[j]] <- data.frame(SUB='1+ muQTL',cor=as.numeric(res$estimate),LOW=res$conf.int[1],HI=res$conf.int[2],p=res$p.value)
j=j+1

i <- which(df.mg$P.MEAN.1 < 5e-8 & df.mg$P.MEAN.1 < 5e-8)
res <- cor.test(df.mg$BETA_INT.20[i],df.mg$BETA_INT.80[i])
df[[j]] <- data.frame(SUB='2 muQTL',cor=as.numeric(res$estimate),LOW=res$conf.int[1],HI=res$conf.int[2],p=res$p.value)
j=j+1

i <- which(df.mg$P.VAR.RAW.1 < 5e-8 | df.mg$P.VAR.RAW.2 < 5e-8)
res <- cor.test(df.mg$BETA_INT.20[i],df.mg$BETA_INT.80[i])
df[[j]] <- data.frame(SUB='1+ raw vQTL',cor=as.numeric(res$estimate),LOW=res$conf.int[1],HI=res$conf.int[2],p=res$p.value)
j=j+1

i <- which(df.mg$P.VAR.RAW.1 < 5e-8 & df.mg$P.VAR.RAW.2 < 5e-8)
res <- cor.test(df.mg$BETA_INT.20[i],df.mg$BETA_INT.80[i])
df[[j]] <- data.frame(SUB='2 raw vQTL',cor=as.numeric(res$estimate),LOW=res$conf.int[1],HI=res$conf.int[2],p=res$p.value)
j=j+1

i <- which(df.mg$P.VAR.RINT.1 < 1e-5 | df.mg$P.VAR.RINT.2 < 1e-5)
res <- cor.test(df.mg$BETA_INT.20[i],df.mg$BETA_INT.80[i])
df[[j]] <- data.frame(SUB='1+ RINT vQTL',cor=as.numeric(res$estimate),LOW=res$conf.int[1],HI=res$conf.int[2],p=res$p.value)
j=j+1

i <- which(df.mg$P.VAR.RINT.1 < 1e-5 & df.mg$P.VAR.RINT.2 < 1e-5)
res <- cor.test(df.mg$BETA_INT.20[i],df.mg$BETA_INT.80[i])
df[[j]] <- data.frame(SUB='2 RINT vQTL',cor=as.numeric(res$estimate),LOW=res$conf.int[1],HI=res$conf.int[2],p=res$p.value)
j=j+1

i <- which(df.mg$P.VAR.RINT.1 < 5e-8 | df.mg$P.VAR.RINT.2 < 5e-8)
res <- cor.test(df.mg$BETA_INT.20[i],df.mg$BETA_INT.80[i])
df[[j]] <- data.frame(SUB='1+ RINT vQTL (P<5E-8)',cor=as.numeric(res$estimate),LOW=res$conf.int[1],HI=res$conf.int[2],p=res$p.value)
j=j+1

i <- which(df.mg$dispersion_pval.1 < 1e-5 | df.mg$dispersion_pval.2 < 1e-5)
res <- cor.test(df.mg$BETA_INT.20[i],df.mg$BETA_INT.80[i])
df[[j]] <- data.frame(SUB='1+ dQTL',cor=as.numeric(res$estimate),LOW=res$conf.int[1],HI=res$conf.int[2],p=res$p.value)
j=j+1

i <- which(df.mg$dispersion_pval.1 < 1e-5 & df.mg$dispersion_pval.2 < 1e-5)
res <- cor.test(df.mg$BETA_INT.20[i],df.mg$BETA_INT.80[i])
df[[j]] <- data.frame(SUB='2 dQTL',cor=as.numeric(res$estimate),LOW=res$conf.int[1],HI=res$conf.int[2],p=res$p.value)
j=j+1

i <- which(df.mg$dispersion_pval.1 < 5e-8 | df.mg$dispersion_pval.2 < 5e-8)
res <- cor.test(df.mg$BETA_INT.20[i],df.mg$BETA_INT.80[i])
df[[j]] <- data.frame(SUB='1+ dQTL (P<5E-8)',cor=as.numeric(res$estimate),LOW=res$conf.int[1],HI=res$conf.int[2],p=res$p.value)
j=j+1

i <- which((df.mg$P.MEAN.1 < 5e-8 & df.mg$P.VAR.RAW.1 < 5e-8) | (df.mg$P.MEAN.2 < 5e-8 & df.mg$P.VAR.RAW.2 < 5e-8))
res <- cor.test(df.mg$BETA_INT.20[i],df.mg$BETA_INT.80[i])
df[[j]] <- data.frame(SUB='1+ muQTL AND raw vQTL',cor=as.numeric(res$estimate),LOW=res$conf.int[1],HI=res$conf.int[2],p=res$p.value)
j=j+1

i <- which((df.mg$P.MEAN.1 < 5e-8 & df.mg$P.VAR.RAW.1 < 5e-8) & (df.mg$P.MEAN.2 < 5e-8 & df.mg$P.VAR.RAW.2 < 5e-8))
res <- cor.test(df.mg$BETA_INT.20[i],df.mg$BETA_INT.80[i])
df[[j]] <- data.frame(SUB='2 muQTL AND raw vQTL',cor=as.numeric(res$estimate),LOW=res$conf.int[1],HI=res$conf.int[2],p=res$p.value)
j=j+1

i <- which((df.mg$P.VAR.RINT.1 < 1e-5 & df.mg$P.VAR.RAW.1 < 5e-8) | (df.mg$P.VAR.RINT.2 < 1e-5 & df.mg$P.VAR.RAW.2 < 5e-8))
res <- cor.test(df.mg$BETA_INT.20[i],df.mg$BETA_INT.80[i])
df[[j]] <- data.frame(SUB='1+ muQTL AND RINT vQTL',cor=as.numeric(res$estimate),LOW=res$conf.int[1],HI=res$conf.int[2],p=res$p.value)
j=j+1

i <- which(df.mg$SNP1=='rs56094641' | df.mg$SNP2=='rs56094641')
res <- cor.test(df.mg$BETA_INT.20[i],df.mg$BETA_INT.80[i])
df[[j]] <- data.frame(SUB='FTO region (rs56094641)',cor=as.numeric(res$estimate),LOW=res$conf.int[1],HI=res$conf.int[2],p=res$p.value)
j=j+1

i <- which(df.mg$P.80 < 0.001)
res <- cor.test(df.mg$BETA_INT.20[i],df.mg$BETA_INT.80[i])
df[[j]] <- data.frame(SUB='P < 0.001',cor=as.numeric(res$estimate),LOW=res$conf.int[1],HI=res$conf.int[2],p=res$p.value)
j=j+1

df.save <- do.call(rbind,df)
f <- '/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/GxG_correlation_between_disc_and_valid.txt'
fwrite(df.save,f,quote = F,sep = '\t',na = 'NA',row.names = F,col.names = T)


# 
# cor.test(apply(df.mg[i,c('BETA.VAR.RAW.1','BETA.VAR.RAW.2')],1,mean),df.mg$BETA_INT.80[i])
# cor.test(apply(df.mg[i,c('BETA.VAR.RINT.1','BETA.VAR.RINT.2')],1,mean),df.mg$BETA_INT.80[i])
# cor.test(apply(df.mg[i,c('dispersion.1','dispersion.2')],1,mean),df.mg$BETA_INT.80[i])
# 
# cor.test(apply(df.mg[i,c('BETA.MEAN.1','BETA.MEAN.2')],1,mean),df.mg$BETA_INT.80[i])
# cor.test(apply(df.mg[i,c('BETA.VAR.RAW.1','BETA.VAR.RAW.2')],1,mean),df.mg$BETA_INT.80[i])
# cor.test(apply(df.mg[i,c('BETA.VAR.RINT.1','BETA.VAR.RINT.2')],1,mean),df.mg$BETA_INT.80[i])
# cor.test(apply(df.mg[i,c('dispersion.1','dispersion.2')],1,mean),df.mg$BETA_INT.80[i])
# 
# cor.test(-log10(df.mg$P.VAR.RAW.1[i]),df.mg$BETA_INT.80[i])
# cor.test(-log10(df.mg$P.VAR.RAW.2[i]),df.mg$BETA_INT.80[i])
# 
# 
# cor.test(df.mg$BETA.VAR.RAW.1,df.mg$BETA_INT.80)
