BATCHID=5011143
dir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/scripts/GWAS/subset


grep Killed $dir/slurm-${BATCHID}_*out | cut -f1 -d':' | cut -f3 -d'_' | cut -f1 -d'.' > $dir/Killed.txt

while read -r line; do
# echo $line
sbatch /athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/scripts/FINAL/gwas/HLMM/errors/run_HLMM_subset_specific.sh bmi.rint.ALL $line
done < "$dir/Killed.txt"