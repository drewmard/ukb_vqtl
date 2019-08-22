library(data.table)

# read in data
pheno1 <- fread('/home/kulmsc/athena/ukbiobank/phenotypes/ukb26867.csv.gz',data.table=F,stringsAsFactors = F)
pheno2 <- fread('/home/kulmsc/athena/ukbiobank/setup_morePhenos/ukb33822.csv.gz',data.table=F,stringsAsFactors = F)

# extract specific data
x <- c('22001-0.0','21022-0.0',
         '864-0.0','874-0.0','884-0.0','894-0.0','904-0.0',
         '914-0.0','1090-0.0','1080-0.0','1070-0.0','1239-0.0',
         '1249-0.0')
i <- which(x %in% colnames(pheno2))
COL_NAMES <- c('sex','age',
  'DayW','DurW','DayM','DurM','DayV',
  'DurV','TimeD','TimeC','TimeTV','CurS','PastS')
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

# smoking:
envir_data$Smoking.E <- rep(NA,nrow(cov))
envir_data$Smoking.E[which( cov$CurS==0 & (cov$PastS==3 | cov$PastS==4) )] <- 0
envir_data$Smoking.E[which( (cov$CurS==1 | cov$CurS==2) | (cov$PastS==1 | cov$PastS==2) )] <- 1

# need to impute medians



########################################################################

x <- c('1160-0.0',
       '1170-0.0',
       '1180-0.0',
       '1190-0.0',
       '1050-0.0','1060-0.0',
       '1269-0.0','1279-0.0',
       '1558-0.0',
       '1737-0.0',
       '3436-0.0','2867-0.0',
       '2814-0.0',
       '6144-0.0'
)
COL_NAMES <- c('sleep.duration',
               'getting.up.morning',
               'morning.evening.person',
               'nap.during.day',
               'time.spent.outdoors.summer','time.spent.outdoors.winter',
               'tobacco.smoke.home','tobacco.smoke.nothome',
               'alcohol.freq',
               'childhood.sunburn.occasions',
               'age.started.smoking.current','age.started.smoking.former',
               'hormone.replacement.therapy',
               'never.eat'
)

cov <- pheno2[,c('eid',
                    x
)]
colnames(cov)[2:ncol(cov)] <- COL_NAMES

envir_data$sleep.duration <- cov$sleep.duration
envir_data$sleep.duration[envir_data$sleep.duration %in% c(-1,-3)] <- NA
envir_data$sleep.duration[which(abs(scale(envir_data$sleep.duration)) > 5)] <- NA

# envir_data$sleep.duration[which(is.na(envir_data$sleep.duration))] <- median(envir_data$sleep.duration,na.rm = T)

envir_data$getting.up.morning <- cov$getting.up.morning
envir_data$getting.up.morning[envir_data$getting.up.morning %in% c(-1,-3)] <- NA
# envir_data$getting.up.morning[which(is.na(envir_data$getting.up.morning))] <- median(envir_data$getting.up.morning,na.rm = T)

envir_data$nap.during.day <- cov$nap.during.day
envir_data$nap.during.day[envir_data$nap.during.day %in% c(-3)] <- NA

envir_data$time.spent.outdoors.summer <- cov$time.spent.outdoors.summer
envir_data$time.spent.outdoors.summer[envir_data$time.spent.outdoors.summer == -10] <- 0
envir_data$time.spent.outdoors.summer[envir_data$time.spent.outdoors.summer %in% c(-1,-3)] <- NA
# envir_data$time.spent.outdoors.summer[which(abs(scale(envir_data$time.spent.outdoors.summer)) > 5)] <- NA

envir_data$time.spent.outdoors.winter <- cov$time.spent.outdoors.winter
envir_data$time.spent.outdoors.winter[envir_data$time.spent.outdoors.winter == -10] <- 0
envir_data$time.spent.outdoors.winter[envir_data$time.spent.outdoors.winter %in% c(-1,-3)] <- NA
# envir_data$time.spent.outdoors.winter[which(abs(scale(envir_data$time.spent.outdoors.winter)) > 5)] <- NA

envir_data$time.spent.outdoors <- apply(envir_data[,c('time.spent.outdoors.summer','time.spent.outdoors.summer')],1,mean,na.rm=T)
envir_data$time.spent.outdoors[which(is.na(cov$time.spent.outdoors.summer) & is.na(cov$time.spent.outdoors.winter))] <- NA
# envir_data$time.spent.outdoors[which(abs(scale(envir_data$time.spent.outdoors)) > 5)] <- NA

cov[,'tobacco.smoke.home'][which(cov[,'tobacco.smoke.home'] %in% c(-1,-3))] <- NA
cov[,'tobacco.smoke.nothome'][which(cov[,'tobacco.smoke.nothome'] %in% c(-1,-3))] <- NA
envir_data$tobacco.smoke.exposure <- apply(cov[,c('tobacco.smoke.home','tobacco.smoke.nothome')],1,sum,na.rm=T)
envir_data$tobacco.smoke.exposure <- as.numeric(envir_data$tobacco.smoke.exposure > 0)

envir_data$alcohol.freq.E <- cov$alcohol.freq;
envir_data$alcohol.freq.E[which(envir_data$alcohol.freq.E %in% c(-3))] <- NA
# envir_data$alcohol.freq.dummy <- as.numeric(is.na(envir_data$alcohol.freq))
# envir_data$alcohol.freq[which(is.na(envir_data$alcohol.freq))] <- median(envir_data$alcohol.freq,na.rm=T)

envir_data$childhood.sunburn.occasions <- cov$childhood.sunburn.occasions
envir_data$childhood.sunburn.occasions[which(envir_data$childhood.sunburn.occasions %in% c(-1,-3))] <- NA
envir_data$childhood.sunburn.occasions[which(abs(scale(envir_data$childhood.sunburn.occasions)) > 5)] <- NA

cov[,'age.started.smoking.current'][which(cov[,'age.started.smoking.current'] %in% c(-1,-3))] <- NA
cov[,'age.started.smoking.former'][which(cov[,'age.started.smoking.former'] %in% c(-1,-3))] <- NA
envir_data$age.started.smoking <- apply(cov[,c('age.started.smoking.current','age.started.smoking.former')],1,mean,na.rm=T)

envir_data$hormone.replacement.therapy <- cov$hormone.replacement.therapy
envir_data$hormone.replacement.therapy[which(envir_data$hormone.replacement.therapy %in% c(-1,-3))] <- NA



# envir_data$never.eat <- cov$envir_data
# (envir_data$never.eat==)


############################################################################################

# x <- c('1289-0.0','1299-0.0','1309-0.0','1319-0.0','1329-0.0','1339-0.0',
#        '1349-0.0','1359-0.0','1369-0.0','1379-0.0','1389-0.0',
#        '1408-0.0','1418-0.0','1428-0.0','1438-0.0','1448-0.0','1458-0.0',
#        '1478-0.0','1488-0.0','1498-0.0'
# )
# 
# COL_NAMES <- c('cooked.veggie','raw.veggie','fresh.fruit','dried.fruit','oily.fish','non.oily.fish',
#                'processed.meat','poultry','beef','lamb','pork',
#                'cheese','milk','spread','bread.intake','bread.type','cereal.intake',
#                'salt.added.to.food','tea','coffee'
# )
# 
# # i <- which(x %in% colnames(pheno1))
# i <- which(x %in% colnames(pheno2))
# 
# cov.v3 <- pheno2[,c('eid',
#                     x[i]
# )]
# 
# colnames(cov.v3)[2:ncol(cov.v3)] <- COL_NAMES[i]




########################################################################

# Pull 2

# x <- c('20458-0.0',
#               '20161-0.0',
#               '20404-0.0','20406-0.0','20415-0.0',
#               '22614-0.0',
#               '22615-0.0',
#               '41202-0.0',
#               '135-0.0',
#               '134-0.0',
#               '1558-0.0',
#               '2453-0.0', '2443-0.0',
#               '2306-0.0',
#               '2463-0.0',
#               '2178-0.0',
#        '2473-0.0',
#        '23104-0.0'
#               )
# COL_NAMES <- c('happiness',
#                'smoking.years',
#                'alcohol.phys.dependent','alcohol.addict','alcohol.curr.addict',
#                'worked.with.pesticides',
#                'diesel.exhaust',
#                'diseases',
#                'number.noncancer.illnesses',
#                'number.cancers',
#                'alcohol.intake.freq',
#                'have.had.cancer','diabetes',
#                'weight.change',
#                'broken.bones.last.5.yr',
#                'overall.health',
#                'serious.medical.condition.diagnosed',
#                'bmi'
#                )
# i <- which(x %in% colnames(pheno1))
# # i <- which(x %in% colnames(pheno2))
# 
# cov.v2 <- pheno1[,c('eid',
#                  x[i]
# )]
# 
# colnames(cov.v2)[2:ncol(cov.v2)] <- COL_NAMES[i]
# 
# cov.v2.2 <- data.frame(eid=cov.v2$eid)
# #
# cov.v2.2$happiness <- NA;  cov.v2.2$happiness[cov.v2$happiness >= 1] <- 3; cov.v2.2$happiness[cov.v2$happiness == 3] <- 2; cov.v2.2$happiness[cov.v2$happiness>=4] <- 1;
# #
# cov.v2.2$smoking.years <- cov.v2$smoking.years; cov.v2.2$smoking.years[which(abs(scale(cov.v2.2$smoking.years)) > 5)] <- NA
# #
# cov.v2.2$alcohol.addiction <- 0; 
# cov.v2.2$alcohol.addiction[which(cov.v2$alcohol.phys.dependent==1 | cov.v2$alcohol.addict==1 | cov.v2$alcohol.curr.addict==1)] <- 1
# cov.v2.2$alcohol.addiction[which(is.na(cov.v2$alcohol.phys.dependent) & is.na(cov.v2$alcohol.addict) & is.na(cov.v2$alcohol.curr.addict))] <- NA
# #
# cov.v2.2$worked.with.pesticides <- cov.v2$worked.with.pesticides
# cov.v2.2$worked.with.pesticides[cov.v2$worked.with.pesticides %in% c(-141,-131)] <- 1
# cov.v2.2$worked.with.pesticides[cov.v2$worked.with.pesticides == -121] <- NA
# #
# cov.v2.2$diesel.exhaust <- cov.v2$diesel.exhaust
# cov.v2.2$diesel.exhaust[cov.v2$diesel.exhaust %in% c(-141,-131)] <- 1
# cov.v2.2$diesel.exhaust[cov.v2$diesel.exhaust == -121] <- NA
# #
# head(cov.v2.2$diseases)
# 
# summary(cov.v2$smoking.years)
# 
# #
# 

##################

# envir_data <- merge(cov2,envir_data,by='eid')
fwrite(envir_data,'/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/envir_data.txt',quote=F,col.names = T,row.names = F,sep = '\t',na='NA')


