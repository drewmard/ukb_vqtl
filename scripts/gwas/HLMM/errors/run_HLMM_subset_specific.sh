#!/bin/bash -l
#SBATCH -J HLMM
#SBATCH -n 1
#SBATCH --mem=32G

arg1=$2

# arg1=1
phenoName=$1
# phenoName=lymphocyte.count.rint.ALL
numCores=4

echo "sbatch run_HLMM_subset.sh $phenoName"
echo "$arg1"

echo "Activating environment..."
source activate HLMM

echo "Starting..."

# initialize
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
rm ${OUT}.models.gz

# correct phenotype ($phenoName)
i=$(head -1 /athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/phenotypes_processed.80.NA.txt | tr -s '\t' '\n' | nl -nln |  grep $phenoName | cut -f1)
i=$(($i-2))

# Run HLMM
echo "python2.7 /home/anm2868/hlmm/bin/hlmm_chr.py $BED $START $END $PHENO $OUT --phen_index $i --min_maf 0 --missing_char NA --max_missing 100"
python2.7 /home/anm2868/hlmm/bin/hlmm_chr.py $BED $START $END $PHENO $OUT --phen_index $i --min_maf 0 --missing_char NA --max_missing 100
