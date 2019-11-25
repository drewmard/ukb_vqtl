#!/bin/bash -l
#SBATCH -J PGS
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
gwasFile=$assocdir/ukbb.lymphocyte.count.rint.ALL.results.txt
clumpFile=$assocdir/ukbb.impute.lymphocyte.count.rint.ALL.muGWAS.chr${CHR}.p_${Pthres}.r_${r2thres}.kb_${kbthres}.clumped

awk 'NR==FNR{A[$3];next} $1 in A' $clumpFile $gwasFile > $clumpFile.txt

out=$clumpFile.PGS

plink --bfile $geno --score $clumpFile.txt 1 11 5 sum --out $out

