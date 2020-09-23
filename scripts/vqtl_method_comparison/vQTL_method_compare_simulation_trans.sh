#!/bin/bash -l
#SBATCH -J sim
#SBATCH -n 64
#SBATCH --mem=128G

spack load /kvarx7m # car
spack load /gmwzft7 # quantreg
Rscript /athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/scripts/vqtl_method_comparison/vQTL_method_compare_simulation_trans.R