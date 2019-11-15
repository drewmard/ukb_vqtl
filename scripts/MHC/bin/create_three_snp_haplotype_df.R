create_three_snp_haplotype_df <- function(df.mg,SNP1,SNP2,SNP3) {

# SNP1, SNP2, SNP3 are in chromosomal order
  
# SNP1='rs1265158'
# SNP2='rs3131643'
# SNP3='rs2516491'
  
i1 <- grep(SNP1,colnames(df.mg))[1]
i2 <- grep(SNP2,colnames(df.mg))[1]
i3 <- grep(SNP3,colnames(df.mg))[1]

df.mg[,i1] <- as.character(df.mg[,i1])
x <- strsplit(unlist(lapply(strsplit(df.mg[,i1],':'),function(x) x[1])),'|')
df.mg[,paste0(SNP1,'_1')] <- unlist(lapply(x,function(x) x[1]))
df.mg[,paste0(SNP1,'_2')] <- unlist(lapply(x,function(x) x[3]))

df.mg[,i2] <- as.character(df.mg[,i2])
x <- strsplit(unlist(lapply(strsplit(df.mg[,i2],':'),function(x) x[1])),'|')
df.mg[,paste0(SNP2,'_1')] <- unlist(lapply(x,function(x) x[1]))
df.mg[,paste0(SNP2,'_2')] <- unlist(lapply(x,function(x) x[3]))

df.mg[,i3] <- as.character(df.mg[,i3])
x <- strsplit(unlist(lapply(strsplit(df.mg[,i3],':'),function(x) x[1])),'|')
df.mg[,paste0(SNP3,'_1')] <- unlist(lapply(x,function(x) x[1]))
df.mg[,paste0(SNP3,'_2')] <- unlist(lapply(x,function(x) x[3]))

df.mg[,'Haplotype1'] <- paste0(df.mg[,paste0(SNP1,'_1')],df.mg[,paste0(SNP2,'_1')],df.mg[,paste0(SNP3,'_1')])
df.mg[,'Haplotype2'] <- paste0(df.mg[,paste0(SNP1,'_2')],df.mg[,paste0(SNP2,'_2')],df.mg[,paste0(SNP3,'_2')])
df.mg[,'Haplotype'] <- paste0(df.mg[,'Haplotype1'],'/',df.mg[,'Haplotype2'])

df.mg[,'000'] <- apply(df.mg[,c('Haplotype1','Haplotype2')],1,function(x) {sum(x=='000')})
df.mg[,'100'] <- apply(df.mg[,c('Haplotype1','Haplotype2')],1,function(x) {sum(x=='100')})
df.mg[,'110'] <- apply(df.mg[,c('Haplotype1','Haplotype2')],1,function(x) {sum(x=='110')})
df.mg[,'111'] <- apply(df.mg[,c('Haplotype1','Haplotype2')],1,function(x) {sum(x=='111')})
df.mg[,'101'] <- apply(df.mg[,c('Haplotype1','Haplotype2')],1,function(x) {sum(x=='101')})
df.mg[,'011'] <- apply(df.mg[,c('Haplotype1','Haplotype2')],1,function(x) {sum(x=='011')})
df.mg[,'001'] <- apply(df.mg[,c('Haplotype1','Haplotype2')],1,function(x) {sum(x=='001')})
df.mg[,'010'] <- apply(df.mg[,c('Haplotype1','Haplotype2')],1,function(x) {sum(x=='010')})

create_haplotype <- function(x) {
  sorting.index <- x[3] > x[4]
  if (sorting.index) {
    res <- paste0(x[1],'/',x[2])
  } else {
    res <- paste0(x[2],'/',x[1])
  }
}

df.mg$Haplotype1.ind <- as.numeric(as.factor(df.mg$Haplotype1))
df.mg$Haplotype2.ind <- as.numeric(as.factor(df.mg$Haplotype2))
df.mg$HaplotypeCode <- apply(df.mg[,c('Haplotype1','Haplotype2','Haplotype1.ind','Haplotype2.ind')],1,create_haplotype)

df.mg <- subset(df.mg,(!is.na(lymphocyte.count.rint.ALL)))

return(df.mg)

}