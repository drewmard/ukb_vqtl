#!/bin/bash -l
#SBATCH -J clump
#SBATCH --mem=128G
#SBATCH --array=1-22:1

source activate vQTL

CHR=$SLURM_ARRAY_TASK_ID
Pthres=0.005
r2thres=0.8
kbthres=250

genodir=/athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/andrew_copies
geno=$genodir/ukbb.${CHR}.impute

assocdir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/imputed/results
assoc=$assocdir/ukbb.lymphocyte.count.rint.ALL.results.txt
out=$assocdir/ukbb.impute.lymphocyte.count.rint.ALL.muGWAS.chr${CHR}.p_${Pthres}.r_${r2thres}.kb_${kbthres}

plink --bfile $geno --clump $assoc --clump-p1 0.005 --clump-r2 0.8 --clump-kb 250 --out $out