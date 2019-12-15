#!/bin/bash -l
#SBATCH -J vQTL
#SBATCH -n 4
#SBATCH --mem=32G
#SBATCH --array=1-5:1

arg1=$SLURM_ARRAY_TASK_ID
phenoName=lymphocyte.count.rint.ALL
numCores=4

echo "sbatch run_vGWAS_subset.sh $arg1 $phenoName $numCores"

echo "Activating environment..."
source activate vQTL

echo "Starting..."
spack load -r r@3.6.0
Rscript /home/anm2868/scripts/UKB/061719_Rserve.R

# 3 : variance
echo "vQTL testing..."
sbatch run_vGWAS_subset.sh $arg1 $phenoName $numCores

echo "Complete."
