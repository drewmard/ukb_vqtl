library(data.table)

# read in data
# pheno1 <- fread('/home/kulmsc/athena/ukbiobank/phenotypes/ukb26867.csv.gz',data.table=F,stringsAsFactors = F)
pheno2 <- fread('/home/kulmsc/athena/ukbiobank/setup_morePhenos/ukb33822.csv.gz',data.table=F,stringsAsFactors = F)

# extract specific data
x <- c('22001-0.0','21022-0.0',
         '864-0.0','874-0.0','884-0.0','894-0.0','904-0.0','914-0.0',
       '1090-0.0','1080-0.0','1070-0.0',
       '1239-0.0','1249-0.0')
i <- which(x %in% colnames(pheno2))
COL_NAMES <- c('sex','age',
  'DayW','DurW','DayM','DurM','DayV','DurV',
  'TimeD','TimeC','TimeTV',
  'CurS','PastS')
cov <- pheno2[,c('eid',
                x[i]
)]
colnames(cov)[2:ncol(cov)] <- COL_NAMES[i]

envir_data <- data.frame(eid=cov$eid)

# physical activity
METW=3.3*cov$DayW*cov$DurW
METM=4.0*cov$DayM*cov$DurM
METV=8.0*cov$DayV*cov$DurV
METT=METW+METM+METV
envir_data$PA <- rep(1,nrow(cov))
envir_data$PA[which((cov$DayV >= 3 & cov$DurV >= 20) | 
           (cov$DayM >= 5 & cov$DurM >= 30) |
           (cov$DayW >= 5 & cov$DurW >= 30) |
           (cov$DayM + cov$DayM + cov$DayV >= 5 & METT >= 600)
)] <- 2
envir_data$PA[which((cov$DayV >= 3 & METT >= 1500) | 
        (cov$DayW + cov$DayM + cov$DayV >= 7 & METT >= 3000)
)] <- 3
table(envir_data$PA)

# sedentary behavior:
TimeD <- cov$TimeD
TimeC <- cov$TimeC
TimeTV <- cov$TimeTV

TimeD[cov$TimeD==-10] <- 0
TimeD[cov$TimeD==-1] <- NA
TimeD[cov$TimeD==-3] <- NA
TimeD[which(is.na(cov$TimeD))] <- median(TimeD,na.rm = T)

TimeC[cov$TimeC==-10] <- 0
TimeC[cov$TimeC==-1] <- NA
TimeC[cov$TimeC==-3] <- NA
TimeC[which(is.na(cov$TimeC))] <- median(TimeC,na.rm = T)

TimeTV[cov$TimeTV==-10] <- 0
TimeTV[cov$TimeTV==-1] <- NA
TimeTV[cov$TimeTV==-3] <- NA
TimeTV[which(is.na(cov$TimeTV))] <- median(TimeTV,na.rm = T)

envir_data$SB <- TimeD + TimeC + TimeTV
envir_data$SB[which(abs(scale(envir_data$SB)) > 5)] <- NA
table(envir_data$SB)

# smoking:
envir_data$Smoking.E <- rep(NA,nrow(cov))
envir_data$Smoking.E[which( cov$CurS==0 & (cov$PastS==3 | cov$PastS==4) )] <- 0
envir_data$Smoking.E[which( (cov$CurS==1 | cov$CurS==2) | (cov$PastS==1 | cov$PastS==2) )] <- 1


########################################################################

x <- c(
       '1558-0.0'
)
COL_NAMES <- c(
               'alcohol.freq'
)

cov <- pheno2[,c('eid',
                    x
)]
colnames(cov)[2:ncol(cov)] <- COL_NAMES

envir_data$Alcohol_intake_frequency <- cov$alcohol.freq;
envir_data$Alcohol_intake_frequency[which(envir_data$Alcohol_intake_frequency %in% c(-3))] <- NA

x <- c('1289-0.0','1299-0.0','1309-0.0','1319-0.0','1329-0.0','1339-0.0',
       '1349-0.0','1359-0.0','1369-0.0','1379-0.0','1389-0.0',
       '1408-0.0','1418-0.0','1428-0.0','1438-0.0','1448-0.0','1458-0.0',
       '1478-0.0','1488-0.0','1498-0.0'
)

COL_NAMES <- c('cooked.veggie','raw.veggie','fresh.fruit','dried.fruit','oily.fish','non.oily.fish',
               'processed.meat','poultry','beef','lamb','pork',
               'cheese','milk','spread','bread.intake','bread.type','cereal.intake',
               'salt.added.to.food','tea','coffee'
)
i <- which(x %in% colnames(pheno2))
cov.v3 <- pheno2[,c('eid',
                    x[i]
)]

colnames(cov.v3)[2:ncol(cov.v3)] <- COL_NAMES[i]


##################

envir_data <- merge(envir_data,cov.v3,by='eid')
fwrite(envir_data,'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/envir_data.txt',quote=F,col.names = T,row.names = F,sep = '\t',na='NA')


