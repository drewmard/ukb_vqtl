#!/bin/bash -l
#SBATCH -J clump
#SBATCH --mem=128G
#SBATCH --array=21-22:1

source activate vQTL

CHR=$SLURM_ARRAY_TASK_ID

genodir=/athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/andrew_copies
geno=$genodir/ukbb.${CHR}.impute

assocdir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/imputed/results
assoc=$assocdir/ukbb.lymphocyte.count.rint.ALL.results.txt
out=$assocdir/ukbb.impute.lymphocyte.count.rint.ALL.muGWAS

Pthres=0.005
r2thres=0.8
kbthres=250

plink --bfile $geno --clump $assoc --clump-p1 0.005 --clump-r2 0.8 --clump-kb 250 --out $out.p_${Pthres}.r_${r2thres}.kb_{kbthres}

