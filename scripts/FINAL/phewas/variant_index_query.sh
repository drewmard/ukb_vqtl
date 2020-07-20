# for uk biobank snps

SNPLIST_OLD=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/sig_results/bmi.SNP_sig_list.txt
SNPLIST=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/sig_results/bmi.SNP_sig_list2.txt
cat $SNPLIST_OLD | tail -n +2 | head -n -1 > $SNPLIST
output=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/sig_results/bmi.SNP_sig_list.query.txt

rm $output
while read -r SNP; do
if [[ "$SNP" == *":"* ]]; then
SNP=$(echo $SNP | sed 's/:/_/g')
SNP=${SNP}_b37
else
SNP=${SNP},
fi
query=$(grep -m 1 $SNP /athena/elementolab/scratch/anm2868/open_targets/csv_files/variant_index.csv)
echo $query >> $output
done < $SNPLIST

