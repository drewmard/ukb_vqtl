#!/bin/bash -l
#SBATCH -J epi
#SBATCH --mem=64G

source activate vQTL

# 4:
phenoName=lymphocyte.count.rint.ALL
outdir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG
prefix=ukbb.${phenoName}.ALL.sub
merged_outFile3=$outdir/ukbb.${phenoName}.results.var.clumped.cut.2.txt
pheno=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/phenotypes_processed.80.txt
plink --bfile $outdir/$prefix --pheno $pheno --pheno-name $phenoName --epistasis set-by-set --set ${merged_outFile3} --epi1 1 --allow-no-sex --threads 1 --memory 56000 --out $outdir/$prefix.GxG


