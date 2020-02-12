#!/bin/bash -l
#SBATCH -J vQTL
#SBATCH -n 1
#SBATCH --mem=16G
#SBATCH --array=1-3:1
# #SBATCH --array=123,411,413
# #SBATCH --array=1-1857:1

arg1=$SLURM_ARRAY_TASK_ID
# arg1=1
phenoName=$1
# phenoName=lymphocyte.count.rint.ALL
numCores=4

echo "sbatch run_HLMM_subset.sh $phenoName"
echo "$arg1"

echo "Activating environment..."
source activate HLMM

echo "Starting..."
# spack load -r r@3.6.0

x1=$( awk 'NR == $$arg1 {print $1}' /athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/andrew_copies/subset/ID.impute.txt )
x2=$( awk 'NR == $$arg1 {print $2}' /athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/andrew_copies/subset/ID.impute.txt )
prefix=/athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/andrew_copies/subset/ukbb.$x1.$x2.impute
BED=$prefix.bed
FAM=$prefix.fam
BIM=$prefix.bim
START=0 
END=$( wc -l $BIM | cut -f1 -d' ' )
PHENO=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/phenotypes_processed.80.NA.txt
OUT=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/HLMM_results/ukbb.$x1.$x2.impute.HLMM_results.txt

# OUT=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/HLMM_results/ukbb.$x1.$x2.impute.HLMM_results.tmp.txt
# START=4980
# END=5000

# 
i=57 # hard coded for now
#
/home/anm2868/hlmm/bin/hlmm_chr.py $BED $START $END $PHENO $OUT --phen_index $i --min_maf 0 --missing_char NA --max_missing 100
