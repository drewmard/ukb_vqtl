ind <- which.min(gwas$P.MEAN.i)
gwas[ind,]
i.to_subset <- gwas[ind,'i']

df.mg[i.to_subset,]
df.mg[order(df.mg$true)[1:20],]

gwas.sub <- subset(gwas,(i==i.to_subset | j==i.to_subset) & (causal.i==1 & causal.j==1))
i.to_subset=333
# i.to_subset=48
gwas.sub <- subset(gwas,(i==i.to_subset | j==i.to_subset))
min(gwas.sub$P.GxG)
table(sign(gwas.sub$BETA.GxG))
# table(sign(gwas.sub$true_beta))
sum(gwas.sub$true_beta)
sum(gwas.sub$BETA.GxG)

ggplot(gwas.sub,aes(x=BETA.GxG,y=true_beta)) + geom_point() + geom_smooth(method='lm',se=F,col='red') + theme_bw()

subset(df.mg,i==i.to_subset)
