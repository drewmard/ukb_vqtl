#!/bin/bash -l
#SBATCH -J vQTL
#SBATCH -n 8
#SBATCH --mem=64G

phenoName=$1
CHR=$2

echo "sbatch run_vGWAS_4.sh $1 $2"

echo "Activating environment..."
source activate vQTL

echo "Starting..."
spack load -r r@3.5.0
Rscript /home/anm2868/scripts/UKB/061719_Rserve.R

pheno=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/phenotypes_processed.80.txt
dir=/home/kulmsc/athena/ukbiobank/calls
outdir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/results
workdir=/athena/elementolab/scratch/anm2868/vQTL

# for CHR in {1..22};
# do

echo "Chromosome $CHR analysis..."
prefix=ukbb.$CHR

# 3 : variance
echo "vQTL testing..."
mkdir -p $outdir
myscript=${workdir}/ukb_vqtl/scripts/GWAS/vGWAS.R
plink --bfile $dir/$prefix --pheno $pheno --pheno-name $phenoName --R $myscript --threads 1 --out $outdir/$prefix.$phenoName.vGWAS
# plink --bfile $dir/$prefix --pheno $pheno --all-pheno --R $myscript --threads 1 --out $outdir/$prefix.vGWAS

# done

echo "Complete."
