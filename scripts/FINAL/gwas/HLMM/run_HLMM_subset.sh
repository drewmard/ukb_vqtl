#!/bin/bash -l
#SBATCH -J HLMM
#SBATCH -n 1
#SBATCH --mem=32G
#SBATCH --array=1-1857:1
# #SBATCH --array=220-229:1,88-97:1,542
# #SBATCH --array=123,411,413
# #SBATCH --array=1-1857:1

arg1=$SLURM_ARRAY_TASK_ID
phenoName=$1
numCores=4

echo "sbatch run_HLMM_subset.sh $phenoName"
echo "$arg1"

echo "Activating environment..."
source activate HLMM

echo "Starting..."
# spack load -r r@3.6.0

x1=$( awk -v var="$arg1" 'NR == var {print $1}' /athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/andrew_copies/subset/ID.impute.txt )
x2=$( awk -v var="$arg1" 'NR == var {print $2}' /athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/andrew_copies/subset/ID.impute.txt )
prefix=/athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/andrew_copies/subset/ukbb.$x1.$x2.impute
BED=$prefix.bed
FAM=$prefix.fam
BIM=$prefix.bim
START=0 
END=$( wc -l $BIM | cut -f1 -d' ' )
PHENO=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/phenotypes_processed.80.NA.txt
OUT=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/HLMM_results/ukbb.$x1.$x2.impute.$phenoName.HLMM_results.txt

# 
i=$(head -1 /athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/phenotypes_processed.80.NA.txt | tr -s '\t' '\n' | nl -nln |  grep $phenoName | cut -f1)
i=$(($i-2))
#
echo "python2.7 /home/anm2868/hlmm/bin/hlmm_chr.py $BED $START $END $PHENO $OUT --phen_index $i --min_maf 0 --missing_char NA --max_missing 100"
python2.7 /home/anm2868/hlmm/bin/hlmm_chr.py $BED $START $END $PHENO $OUT --phen_index $i --min_maf 0 --missing_char NA --max_missing 100
