library(data.table)
f <- '~/Downloads/2020-04-22_PublicID_iTwinID_to_Andrew.txt'
df <- fread(f,data.table=F,stringsAsFactors = F)
f <- '~/Documents/Research/vQTL/ukb_vqtl/scripts/UKT/2020-04-02_16S_annotation.txt'
annot <- fread(f,data.table = F,stringsAsFactors = F)


sum(duplicated(annot$i.TwinID))
sum(duplicated(annot$i.IndividualID))
sum(annot$i.TwinID!=annot$i.IndividualID)

tail(annot[,c('i.TwinID','i.IndividualID')])
