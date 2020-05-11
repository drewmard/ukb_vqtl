#!/bin/bash -l
#SBATCH -J gxe
#SBATCH -n 32
#SBATCH --mem=64G

Rscript /athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/scripts/FINAL/gxe/GxE_updated.matched.R
