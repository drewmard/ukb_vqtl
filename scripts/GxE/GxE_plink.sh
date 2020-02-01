source activate vQTL

phenoName=neutrophil.count.rint.ALL

for covariate in PA SB Smoking.E sleep.duration time.spent.outdoors tobacco.smoke.exposure alcohol.freq.E age sex;
# for covariate in age;
do

if [ "$covariate" == "Smoking.E" ]; then
	num=2
elif [ "$covariate" == "alcohol.freq.E" ]; then
	num=3
elif [ "$covariate" == "age" ]; then
	num=4
	echo $covariate 
	echo "woo!"
else
	num=1
fi

# num=1
# covariate=sex
suffix=resid.$num
phenoAnalysisName=$phenoName.$suffix

pheno=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/EnvData/pheno.80.txt
indir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_GxG
prefix=ukbb.${phenoName}.ALL.sub

mkdir -p /athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_res
outFile=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/GxE_res/$prefix.$covariate.GxE.80

covarData=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxE/EnvData/$covariate.80.txt
plink --bfile $indir/$prefix --pheno $pheno --pheno-name $phenoAnalysisName --linear --covar $covarData --interaction --allow-no-sex --out $outFile

done
