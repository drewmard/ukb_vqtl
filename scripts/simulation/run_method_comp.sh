#!/bin/bash -l
#SBATCH -J vQTLcomp
#SBATCH --mem=64G
#SBATCH -N 8
#SBATCH -n 8

source activate vQTL
Rscript /athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/scripts/simulation/vQTL_method_compare_simulation.R
