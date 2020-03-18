# input
ID=/athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/andrew_copies/subset/ID.impute.txt
# ID=/athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/andrew_copies/subset/ID.impute.txt.tmp # testing!

phenoName=$1
# phenoName=monocyte.count.rint.ALL # hard input

#initialize
outFile=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_subset/ukbb.$phenoName.vGWAS.txt
rm $outFile
head -1 /athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_subset/ukbb.1.1.$phenoName.txt > $outFile

# loop
i=0
while IFS=$'\t' read -r -a myArray
do
 # loop initialize
 x1=${myArray[0]} #chr
 x2=${myArray[1]} #iterator
 i=$(($i+1))


#  x1=11
#  x2=29
#  f=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_subset/ukbb.$x1.$x2.$phenoName.txt
#  # print to script
#  if [ $(wc -l $f | cut -f1 -d' ') -gt 5001 ];
#  then
#    echo "error"
#  fi
#  
 
 if [ ! -f $f ] || [ $(wc -l $f | cut -f1 -d' ') -gt 5001 ]; 
 then
   echo "ERROR: $i ($x1, $x2)"
   sbatch /athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/scripts/GWAS/subset/run_vGWAS_subset_specific.sh $phenoName $i
 else
   tail -n +2 $f >> $outFile
 fi
 
done < $ID

