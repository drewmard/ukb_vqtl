#!/bin/bash -l
#SBATCH -J vQTL
#SBATCH --mem=128G
#SBATCH --array=1-22:1

echo "Activating environment..."
source activate vQTL

echo "Starting..."
spack load -r r@3.5.0
Rscript /home/anm2868/scripts/UKB/061719_Rserve.R
pheno=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/phenotypes_processed.80.txt
dir=/athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/andrew_copies

outdir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/imputed
mkdir -p $outdir
mkdir -p $outdir/MAF
mkdir -p $outdir/results

# phenoName=lymphocyte.count.ALL
phenoName=$1
# phenoName=bmi.ALL
CHR=$SLURM_ARRAY_TASK_ID

# for CHR in {1..22};
# do

echo "Chromosome $CHR analysis..."
prefix=ukbb.$CHR.impute

# 1 : freq
# echo "MAF measurements..."
# outdir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/imputed/MAF
# plink --bfile $dir/$prefix --freq --allow-no-sex --memory 120000 --out $outdir/$prefix

# GWAS: ##########

# 2 : mean
echo "muQTL testing..."
outdir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/imputed/results
plink --bfile $dir/$prefix --pheno $pheno --pheno-name $phenoName --assoc --out $outdir/$prefix.$phenoName.muGWAS

# for some reason, isn't running all phenotypes!
# plink --bfile $dir/$prefix --pheno $pheno --all-pheno --assoc --allow-no-sex --memory 120000 --out $outdir/$prefix.muGWAS

# done


echo "Complete."

# done < /athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/phenotype_names.txt


