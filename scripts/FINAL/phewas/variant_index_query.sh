# for uk biobank snps

pheno=bmi
QTLtype=var.raw
QTLtype=mean
# SNPLIST=/athena/elementolab/scratch/anm2868/open_targets/query/$pheno.$QTLtype.txt
# output=/athena/elementolab/scratch/anm2868/open_targets/query/$pheno.$QTLtype.query.txt

SNPLIST_OLD=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/sig_results/bmi.SNP_sig_list.txt
SNPLIST=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/sig_results/bmi.SNP_sig_list2.txt
cat $SNPLIST_OLD | tail -n +2 | head -n -1 > $SNPLIST
output=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/sig_results/bmi.SNP_sig_list.query.txt

# i=0
rm $output
while read -r SNP; do
# SNP=6:163043152_CACTTT_C
# SNP=rs10423928
# i=$($i+1)
if [[ "$SNP" == *":"* ]]; then
SNP=$(echo $SNP | sed 's/:/_/g')
SNP=${SNP}_b37
else
SNP=${SNP},
fi
# echo $SNP
query=$(grep -m 1 $SNP /athena/elementolab/scratch/anm2868/open_targets/csv_files/variant_index.csv)
# query=$(echo $query | cut -f1,2 -d',')
echo $query >> $output
done < $SNPLIST


# grep -m 1 141202978 /athena/elementolab/scratch/anm2868/open_targets/csv_files/variant_index.csv


# grep -m 1 rs543124501 /athena/elementolab/scratch/anm2868/open_targets/csv_files/variant_index.csv