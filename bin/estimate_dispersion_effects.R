require(MASS)
r_var_mean=rlm(var~0+add,data=results)

mean_noise=mean(results$add_se^2,na.rm=T)
noise_adjustment=1+mean_noise/(var(results$add,na.rm=T)-mean_noise)

r_av=r_var_mean$coefficients[1]*noise_adjustment

results$dispersion=results$var-r_av*results$add
results$dispersion_se=sqrt(results$var_se^2+(r_av^2)*results$add_se^2)
results$dispersion_t=results$dispersion/results$dispersion_se
results$dispersion_pval=-log10(pchisq(results$dispersion_t^2,1,lower.tail=F))
