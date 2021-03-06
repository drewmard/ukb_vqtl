# input
ID=/athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/andrew_copies/subset/ID.impute.txt
# ID=/athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/andrew_copies/subset/ID.impute.txt.tmp # testing!

phenotype=$1
phenotype=bmi.rint.ALL # hard input

#initialize
outFile=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/HLMM_results/ukbb.$phenotype.HLMM.txt
rm $outFile
head -1 /athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/HLMM_results/ukbb.1.1.impute.$phenotype.HLMM_results.txt.models.gz > $outFile

# loop
i=0
while IFS=$'\t' read -r -a myArray
do
 # loop initialize
 x1=${myArray[0]} #chr
 x2=${myArray[1]} #iterator
 i=$(($i+1))

 f=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/HLMM_results/ukbb.$x1.$x2.impute.$phenotype.HLMM_results.txt.models.gz
 
 # print to script
 if [ ! -f $f ]; then
   echo "ERROR: $i ($x1, $x2)"
   sbatch /athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/scripts/FINAL/gwas/HLMM/errors/run_HLMM_subset_specific.sh $phenotype $i
 else
   tail -n +2 $f >> $outFile
 fi
 
done < $ID


wc -l /athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/HLMM_results/ukbb.*.*.impute.$phenotype.HLMM_results.txt.models.gz > /athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/HLMM_results/wordcounts.txt
