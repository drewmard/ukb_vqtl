#!/bin/bash -l
#SBATCH -J vQTL
#SBATCH --mem=128G
#SBATCH --array=6-6:1


# #SBATCH --array=2,3,4,5,6,7,8,11,14,15,17,18,22
# #SBATCH --array=1-22:1

phenoName=$1
CHR=$SLURM_ARRAY_TASK_ID

echo "sbatch run_vGWAS_5.sh $phenoName $CHR"

echo "Activating environment..."
source activate vQTL

echo "Starting..."
spack load -r r@3.5.0
Rscript /home/anm2868/scripts/UKB/061719_Rserve.R

pheno=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/phenotypes_processed.80.txt
dir=/home/kulmsc/athena/ukbiobank/calls
outdir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/results
workdir=/athena/elementolab/scratch/anm2868/vQTL

echo "Chromosome $CHR analysis..."
prefix=ukbb.$CHR

# 3 : variance
echo "vQTL testing..."
# mkdir -p $outdir
myscript=${workdir}/ukb_vqtl/scripts/GWAS/vGWAS.R
plink --bfile $dir/$prefix --pheno $pheno --pheno-name $phenoName --R $myscript --threads 1 --memory 96000 --out $outdir/$prefix.$phenoName.vGWAS

echo "Complete."
