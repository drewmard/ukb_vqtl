#!/bin/bash -l
#SBATCH -J vQTL
#SBATCH -n 16
#SBATCH --mem=258G

CHR=$1
echo "Chromosome $CHR analysis..."

echo "Activating environment..."
source activate HLMM

echo "Starting..."
spack load -r r@3.5.0
Rscript /home/anm2868/scripts/UKB/061719_Rserve.R

# input
workdir=/athena/elementolab/scratch/anm2868/vQTL

pheno=${workdir}/ukb_vqtl/output/GWAS/$phenoFile # /athena/elementolab/scratch/anm2868/vQTL/UKB/phenotypes_blood.2.txt
dir=/home/kulmsc/athena/ukbiobank/calls # genotype file location
outdir=${workdir}/ukb_vqtl/output/GWAS/results # /athena/elementolab/scratch/anm2868/vQTL/UKB/results
prefix=ukbb.$CHR
phenoName=$2 # bmi
myscript=${workdir}/ukb_vqtl/scripts/GWAS/vGWAS.R

mkdir -p $outdir

# 1 : freq
echo "MAF measurements..."
plink --bfile $dir/$prefix --freq --out $outdir/$prefix

# 2 : mean
echo "muQTL testing..."
plink --bfile $dir/$prefix --pheno $pheno --pheno-name $phenoName --assoc --out $outdir/$prefix.$phenoName.muGWAS

# 3 : variance
echo "vQTL testing..."
myscript=/athena/elementolab/scratch/anm2868/vQTL/UKB/Rscript/vGWAS.R
# plink --bfile $dir/$prefix --pheno $pheno --pheno-name $phenoName --R $myscript --threads 1 --out $outdir/$prefix.$phenoName.vGWAS
plink --bfile $dir/$prefix --pheno $pheno --pheno-name $phenoName --R $myscript --threads 1 --out $outdir/$prefix.$phenoName.vGWAS

echo "Complete."
