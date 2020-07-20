library(data.table)
f <- '/Users/andrewmarderstein/Documents/Research/vQTL/bmi.80.txt'
df <- fread(f,data.table = F,stringsAsFactors = F,header = T)
# df <- subset(df,Alcohol_intake_frequency %in% 1:6)
# table(df$Alcohol_intake_frequency)

val=1
mod <- list()
for (val in 1:4) {
  mod[[val]]<-lm(bmi.na~SNP,data = subset(df,Alcohol_intake_frequency==val))
}
mod[[5]]<-lm(bmi.na~SNP,data = subset(df,Alcohol_intake_frequency%in%c(5,6)))
gxe_effects <- as.data.frame(do.call(rbind,lapply(1:5,function(i){summary(mod[[i]])$coef[2,c(1,2)]})))
gxe_effects$E <- 'Alc'
gxe_effects$LOW <- gxe_effects$Estimate-1.96*gxe_effects$`Std. Error`
gxe_effects$HI <- gxe_effects$Estimate+1.96*gxe_effects$`Std. Error`
x <- c('Daily','3-4x/wk','1-2x/wk','1-3x/mo','Rarely/never')
gxe_effects$Val <- factor(x,levels = x)

# table(floor(df$age/10)*10)
mod <- list()
mod[[1]] <- lm(bmi.na~SNP,data = subset(df,age<=49))
mod[[2]] <- lm(bmi.na~SNP,data = subset(df,age>=50 & age<=59))
mod[[3]] <- lm(bmi.na~SNP,data = subset(df,age>=60))
gxe_effects.tmp <- as.data.frame(do.call(rbind,lapply(1:3,function(i){summary(mod[[i]])$coef[2,c(1,2)]})))
gxe_effects.tmp$E <- 'Age'
gxe_effects.tmp$LOW <- gxe_effects.tmp$Estimate-1.96*gxe_effects.tmp$`Std. Error`
gxe_effects.tmp$HI <- gxe_effects.tmp$Estimate+1.96*gxe_effects.tmp$`Std. Error`
gxe_effects.tmp$Val <- c('37-49','50-59','60-72')
gxe_effects <- rbind(gxe_effects,gxe_effects.tmp)

mod <- list()
mod[[1]] <- lm(bmi.na~SNP,data = subset(df,PA==1))
mod[[2]] <- lm(bmi.na~SNP,data = subset(df,PA==2))
mod[[3]] <- lm(bmi.na~SNP,data = subset(df,PA==3))
gxe_effects.tmp <- as.data.frame(do.call(rbind,lapply(1:3,function(i){summary(mod[[i]])$coef[2,c(1,2)]})))
gxe_effects.tmp$LOW <- gxe_effects.tmp$Estimate-1.96*gxe_effects.tmp$`Std. Error`
gxe_effects.tmp$HI <- gxe_effects.tmp$Estimate+1.96*gxe_effects.tmp$`Std. Error`
gxe_effects.tmp$E <- 'PA'
gxe_effects.tmp$Val <- c('Low','Mod','High')
gxe_effects <- rbind(gxe_effects,gxe_effects.tmp)

mod <- list()
mod[[1]] <- lm(bmi.na~SNP,data = subset(df,Smoking.E==0))
mod[[2]] <- lm(bmi.na~SNP,data = subset(df,Smoking.E==1))
gxe_effects.tmp <- as.data.frame(do.call(rbind,lapply(1:2,function(i){summary(mod[[i]])$coef[2,c(1,2)]})))
gxe_effects.tmp$LOW <- gxe_effects.tmp$Estimate-1.96*gxe_effects.tmp$`Std. Error`
gxe_effects.tmp$HI <- gxe_effects.tmp$Estimate+1.96*gxe_effects.tmp$`Std. Error`
gxe_effects.tmp$E <- 'Smok'
gxe_effects.tmp$Val <- c('No','Yes')
gxe_effects <- rbind(gxe_effects,gxe_effects.tmp)

mod <- list()
mod[[1]] <- lm(bmi.na~SNP,data = subset(df,sex==0))
mod[[2]] <- lm(bmi.na~SNP,data = subset(df,sex==1))
gxe_effects.tmp <- as.data.frame(do.call(rbind,lapply(1:2,function(i){summary(mod[[i]])$coef[2,c(1,2)]})))
gxe_effects.tmp$LOW <- gxe_effects.tmp$Estimate-1.96*gxe_effects.tmp$`Std. Error`
gxe_effects.tmp$HI <- gxe_effects.tmp$Estimate+1.96*gxe_effects.tmp$`Std. Error`
gxe_effects.tmp$E <- 'Sex'
gxe_effects.tmp$Val <- c('Female','Male')
gxe_effects <- rbind(gxe_effects,gxe_effects.tmp)

mod <- list()
mod[[1]] <- lm(bmi.na~SNP,data = subset(df,SB>=0 & SB<=2))
mod[[2]] <- lm(bmi.na~SNP,data = subset(df,SB>=3 & SB<=4))
mod[[3]] <- lm(bmi.na~SNP,data = subset(df,SB>=5 & SB<=6))
mod[[4]] <- lm(bmi.na~SNP,data = subset(df,SB>=7))
gxe_effects.tmp <- as.data.frame(do.call(rbind,lapply(1:4,function(i){summary(mod[[i]])$coef[2,c(1,2)]})))
gxe_effects.tmp$LOW <- gxe_effects.tmp$Estimate-1.96*gxe_effects.tmp$`Std. Error`
gxe_effects.tmp$HI <- gxe_effects.tmp$Estimate+1.96*gxe_effects.tmp$`Std. Error`
gxe_effects.tmp$E <- 'SB'
gxe_effects.tmp$Val <- c('0-2','3-4','5-6','7+')
gxe_effects <- rbind(gxe_effects,gxe_effects.tmp)


# ggplot(data=gxe_effects, aes(x=E, y=Estimate, ymin=LOW, ymax=HI,fill=as.factor(Val),col=as.factor(Val))) +

  #
# gxe_effects$Val <- as.factor(gxe_effects$Val)
gxe_effects <- gxe_effects[order(gxe_effects$E),]
ggCOL <- c(
  'lightsalmon','lightsalmon1','lightsalmon2','lightsalmon3','lightsalmon4',
  'lightskyblue1','lightskyblue3','lightskyblue4',
  'tomato1','tomato2','tomato4',
  'midnightblue','mediumvioletred',
  'cadetblue1','cadetblue4',
  'bisque1','bisque2','bisque3','bisque4',
  'steelblue1','steelblue2','steelblue3'
)
g <- ggplot(data=gxe_effects, aes(x=as.factor(Val), y=Estimate, ymin=LOW, ymax=HI,fill=as.factor(Val),col=as.factor(Val))) +
  geom_linerange(size=5,position=position_dodge(width = 0.5)) +
  geom_point(size=3, shape=21, colour="white",fill='black', stroke = 0.5,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, lty=2,col='red') +
  coord_flip() +
  labs(x="Environmental Factor",y=expression('Avg BMI increase per rs56094641 G allele ('~kg/m^2~')'),col='QTL Type') +
  theme_bw()  +
  theme(axis.text = element_text(family = 'Helvetica'),legend.position = 'none') +
  facet_grid(E~.,scales='free',space='free') +
  scale_colour_manual(values=ggCOL)

x <- 1
png('~/Documents/Research/vQTL/ukb_vqtl/output/GxE/GxE_results/fto_gxe.png',width=5000,height=2200,res=500)
print(g)
dev.off()

# 
# tmp <- aggregate(gxe_effects$Val,by=list(gxe_effects$E),function(x) {length(unique(x))})


