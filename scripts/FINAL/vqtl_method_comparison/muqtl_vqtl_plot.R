library(ggplot2)
library(data.table)
k=1
results.list <- list()
for (k in 1:2) {
  set.seed(03191995)
  phenotype_noise <- phenotype_noise.vec[k]
  simulation_results <- lapply(1:length(genetic_variance_explained.vec),runSimulation,type=simulation_type)
  results <- do.call(rbind,lapply(simulation_results,function(x){x[[1]]}))
  results <- as.data.frame(results,stringsAsFactors = FALSE)
  
  #1
  results$LR_p.reject <- (results$LR_p<0.05)
  results$DRM_p.reject <- (results$DRM_p<0.05)
  results.list[[k]] <- results
}
results.save <- do.call(rbind,results.list)

tmp <- list(); i = 0
for (pheno_noise in c('NORMAL','CHISQ4')) {
  for (val in c(0.01,0.02,0.03)) {
    i=i+1
    results.sub <- subset(results.save,h==val & noise==pheno_noise)
    tab <- table(results.sub[,c('LR_p.reject','DRM_p.reject')])
    P=(fisher.test(tab)$p.value)
    POWER1=tab[1,2]/(tab[1,1]+tab[1,2])
    POWER2=tab[2,2]/(tab[2,1]+tab[2,2])
    tmp[[i]] <- data.frame(noise=pheno_noise,h=val,P,POWER1,POWER2)
  }
}
tmp.full <- do.call(rbind,tmp)
library(reshape2)
tmp.save <- melt(tmp.full[,c(1,2,4,5)],id.vars=c('noise','h'))
tmp.save$Name <- paste(tmp.save$noise,tmp.save$h,sep=': h = ')
tmp.save$Name <- factor(tmp.save$Name,tmp.save$Name[1:6])

g <- ggplot(tmp.save,aes(x=Name,y=value,fill=variable)) + 
  geom_bar(stat='identity',position = 'dodge',col='black') +
  theme_bw() + 
  labs(x='Simulation',y='vQTL power',fill='muQTL?') +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(hjust=1,angle=30)) +
  scale_fill_manual(values=c('black','steelblue2'),labels=c('No','Yes')) +
  scale_x_discrete(labels=c(expression('Normal, '*V[G]*'=0.01'),
                            expression('Normal, '*V[G]*'=0.02'),
                            expression('Normal, '*V[G]*'=0.03'),
                            expression('Non-normal,'*V[G]*'=0.01'),
                            expression('Non-normal,'*V[G]*'=0.02'),
                            expression('Non-normal,'*V[G]*'=0.03')
  )) +
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),limits=c(0,1.1))
  # lims(y=c(0,1.1));
png('~/Documents/Research/vQTL/ukb_vqtl/output/simulation/barplot_muqtl_vqtl.png',width=4150,height=2660,res=600)
print(g)
dev.off()

# res <- aggregate(results[,c('LR_p.reject','DRM_p.reject')],by = list(h=results$h),mean)
# colnames(res)[2:3] <- c('LR','DRM')
# png('~/Documents/Research/vQTL/ukb_vqtl/output/simulation/muqtl_vs_vqtl.png',width = 2700,height = 2700,res=650)
# ggplot(melt(res,id.vars='h'),aes(x=h,y=value,col=variable)) + geom_line() + geom_point() + theme_bw() + theme(panel.grid=element_blank()) + labs(col='Method',x='Variance explained by GxG',y='Power')
# dev.off()
#2
val <- 0.01
results.sub <- subset(results,h==val)
tab <- table(results.sub[,c('LR_p.reject','DRM_p.reject')])
fisher.test(tab)$p.value

#3
val <- 0.02
results.sub <- subset(results,h==val)
tab <- table(results.sub[,c('LR_p.reject','DRM_p.reject')])
fisher.test(tab)$p.value
tab[1,2]/(tab[1,1]+tab[1,2])
tab[2,2]/(tab[2,1]+tab[2,2])

#4
val <- 0.03
results.sub <- subset(results,h==val)
tab <- table(results.sub[,c('LR_p.reject','DRM_p.reject')])
fisher.test(tab)$p.value
tab[1,2]/(tab[1,1]+tab[1,2])
tab[2,2]/(tab[2,1]+tab[2,2])

