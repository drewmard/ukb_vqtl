#!/bin/bash -l
#SBATCH -J sim
#SBATCH -n 16
#SBATCH --mem=48G

spack load /kvarx7m # car
spack load /gmwzft7 # quantreg
Rscript /athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/scripts/vqtl_method_comparison/vQTL_method_compare_simulation_trans.R