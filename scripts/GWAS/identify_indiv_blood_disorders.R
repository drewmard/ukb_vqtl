library(data.table)
pheno <- fread('/home/kulmsc/athena/ukbiobank/phenotypes/ukb26867.csv.gz',data.table=F,stringsAsFactors = F)
diseases <- pheno[,c('eid','20001-0.0','20002-0.0','40001-0.0','40013-0.0','40011-0.0')]
colnames(diseases)[2:ncol(diseases)] <- c('cancer.code','noncancer.illness.code','ICD10','ICD9','tumor.histology')

i1 <- diseases$noncancer.illness.code %in% c(1658, #myelofibrosis
                                        1495,1448, #lymphoedema
                                        # 1065:1094, #cardiovascular
                                       # 1220:1223, #diabetes
                                        1327:1340,
                                        1546, #essential thrombocytosis
                                       1438, # polycythaemia vera	
                                       1445:1451 # haematology
                                        )
#1461-1463 for IBD
i2 <- diseases$cancer.code %in% c(
  1047:1058
)

i3 <- diseases$tumor.histology %in% c(
  9590:9999
)

i4 <- diseases$ICD10 %in% c(
  paste0('B',20:24),paste0('B',200:238), #HIV
  paste0('C',81:96),paste0('C',810:969),
  paste0('D',70:73),paste0('D',720:729),
  paste0('D',75:77),paste0('D',750:763),
  paste0('D',80:84),paste0('D',800:849), #immune
  'P611','D45','D751','P583',
  'D46','D471',
  paste0('D',45:47),paste0('D',460:479),
  paste0('D',57:64,570:649)
)


i.all <- which(i1 | i2 | i3 | i4)

fwrite(data.frame(row=i.all,eid=diseases$eid[i.all]),'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/blood_disease.indiv_id.txt',col.names = T,row.names = F,quote=F,sep = '\t',na='NA')


