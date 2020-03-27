#!/bin/bash -l
#SBATCH -J vQTL
#SBATCH -n 4
#SBATCH --mem=32G
#SBATCH --array=1-1857:1
# #SBATCH --array=123,411,413
# #SBATCH --array=1-1857:1

arg1=$SLURM_ARRAY_TASK_ID
phenoName=neutrophil.count.rint.ALL
phenoName=$1
numCores=4

echo "sbatch run_vGWAS_subset.sh $arg1 $phenoName $numCores"

echo "Activating environment..."
source activate vQTL

echo "Starting..."
# spack load -r r@3.6.0

# 3 : variance
echo "vQTL testing..."
x1="$(cut -f1 /athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/andrew_copies/subset/ID.impute.txt | head -${arg1} | tail -1)"
x2="$(cut -f2 /athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/andrew_copies/subset/ID.impute.txt | head -${arg1} | tail -1)"
FILE=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_subset/ukbb.$x1.$x2.$phenoName.txt
# if [ ! -f "$FILE" ]; then
	echo "Let it run!"
	Rscript /athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/scripts/GWAS/subset/run_vGWAS.R $arg1 $phenoName $numCores
# fi
echo "Complete."
