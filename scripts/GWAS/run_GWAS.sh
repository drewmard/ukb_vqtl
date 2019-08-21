echo "Activating environment..."
source activate vQTL

echo "Starting..."
spack load -r r@3.5.0
Rscript /home/anm2868/scripts/UKB/061719_Rserve.R

pheno=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/phenotypes_processed.80.txt
dir=/home/kulmsc/athena/ukbiobank/calls
outdir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/results

for CHR in {1..22};
do

echo "Chromosome $CHR analysis..."
prefix=ukbb.$CHR
# phenoName=$2

# 1 : freq
echo "MAF measurements..."
outdir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/MAF
mkdir -p $outdir
plink --bfile $dir/$prefix --freq --out $outdir/$prefix

# GWAS: ##########
outdir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/results

# 2 : mean
echo "muQTL testing..."
mkdir -p $outdir
# plink --bfile $dir/$prefix --pheno $pheno --pheno-name $phenoName --assoc --out $outdir/$prefix.$phenoName.muGWAS
plink --bfile $dir/$prefix --pheno $pheno --all-pheno --assoc --out $outdir/$prefix.muGWAS

done


echo "Complete."

