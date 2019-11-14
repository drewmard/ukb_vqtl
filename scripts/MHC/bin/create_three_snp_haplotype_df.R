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
# df.mg[,'HaplotypeCode'] <- NA
# df.mg[,'HaplotypeCode'][df.mg[,'Haplotype']%in%c('00/00')] <- '00/00'
# df.mg[,'HaplotypeCode'][df.mg[,'Haplotype']%in%c('10/10')] <- '10/10'
# df.mg[,'HaplotypeCode'][df.mg[,'Haplotype']%in%c('01/01')] <- '01/01'
# df.mg[,'HaplotypeCode'][df.mg[,'Haplotype']%in%c('11/11')] <- '11/11'
# df.mg[,'HaplotypeCode'][df.mg[,'Haplotype']%in%c('10/00','00/10')] <- '10/00'
# df.mg[,'HaplotypeCode'][df.mg[,'Haplotype']%in%c('01/00','00/01')] <- '01/00'
# df.mg[,'HaplotypeCode'][df.mg[,'Haplotype']%in%c('11/00','00/11')] <- '11/00'
# df.mg[,'HaplotypeCode'][df.mg[,'Haplotype']%in%c('10/01','01/10')] <- '10/01'
# df.mg[,'HaplotypeCode'][df.mg[,'Haplotype']%in%c('11/01','01/11')] <- '11/01'
# df.mg[,'HaplotypeCode'][df.mg[,'Haplotype']%in%c('11/10','10/11')] <- '11/10'

df.mg <- subset(df.mg,(!is.na(lymphocyte.count.rint.ALL)))

return(df.mg)

}