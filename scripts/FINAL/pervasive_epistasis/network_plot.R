library(parallel)
library(ggplot2)
library(cowplot)
library(stringr)
library(data.table)
set.seed(03191995)

# Initialize
Nindiv=10000
MAF=0.25
Mtotsnp=100
Mcausalsnp=Mtotsnp
Ninteraction <- Mtotsnp * (Mtotsnp-1) / 2

print('Generating genotypes...')
G <- matrix(rbinom(Nindiv*Mtotsnp,2,MAF),nrow = Nindiv,ncol = Mtotsnp)

print('Generating interactions file...')
tmp <- model.matrix( ~.^2, data=as.data.frame(G[,1:Mcausalsnp]))[,-1]
print('Generating phenotypes...')
beta <- rnorm(ncol(tmp))
pheno.interaction <- as.numeric(tmp[,(Mcausalsnp+1):ncol(tmp)] %*% beta[(Mcausalsnp+1):ncol(tmp)])

s <- str_replace_all(colnames(tmp[,(Mcausalsnp+1):ncol(tmp)]),'V','')
SNP1=unlist(lapply(strsplit(s,':'),function(x)x[1]))
SNP2=unlist(lapply(strsplit(s,':'),function(x)x[2]))
BETA <- beta[(Mcausalsnp+1):ncol(tmp)]
df <- data.frame(SNP1,SNP2,BETA)
g <- graph_from_data_frame(df, directed=F)

df.copied <- df
val <- df.copied$SNP1; df.copied$SNP1 <- df.copied$SNP2; df.copied$SNP2 <- val
df.mg <- rbind(df,df.copied)
df.mg.aggre <- aggregate(df.mg$BETA,list(df.mg$SNP1),mean)

E(g)$color <- sapply(sign(df$BETA),function(x) ifelse(x > 0,'red','blue'))
E(g)$weight <- abs(df$BETA^4)
V(g)$color <- sapply(sign(df.mg.aggre$x),function(x) ifelse(x > 0,'red','blue')) 
V(g)$weight <- df.mg.aggre$x
plot(g,edge.width=abs(df$BETA),vertex.label=NA)
