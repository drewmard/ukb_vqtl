# run vGWAS/muGWAS on curie.pbtech:

for i in {1..22}
do
PHENO=lymphocyte.count
sbatch /home/anm2868/scripts/UKB/072019_run_GWAS.slurm.sh $i $PHENO
done
