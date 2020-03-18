# GWAS_to_clumped_SNP_list.sh

pheno=$1
suffix1=$2
suffix2=$3
thres=$4

r2thres=0.5

# initialize directories
genodir=/athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/andrew_copies
if [ "$suffix2" == "var" ]; then
  outdir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_clump # /${phenoName}
  results=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_subset/ukbb.${phenoName}.vGWAS.txt
  clump_snp_field=rs
elif [ "$suffix2" == "mean" ]; then
  outdir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/muGWAS_clump
  results=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/imputed/results/ukbb.${phenoName}.results.txt
  clump_snp_field=SNP
fi

echo 'Clumping!'
# perform clumping
for CHR in {1..22}
do
echo ${CHR}
geno=$genodir/ukbb.${CHR}.impute
outFile=ukbb.${phenoName}.${CHR}.${suffix2}.clumped
plink --bfile $geno --maf 0.05 --clump $results --clump-p1 $thres --clump-p2 $thres --clump-r2 $r2thres --clump-field P --clump-snp-field $clump_snp_field --out $outdir/$outFile
done

# merge diff chr clump files
dir=$outdir
merged_outFile=$dir/ukbb.${phenoName}.${suffix2}.clumped.clumped.combined.txt
head -1 $dir/ukbb.${phenoName}.${CHR}.${suffix2}.clumped.clumped > $merged_outFile
for CHR in {1..22}
do
echo $CHR
geno=ukbb.${CHR}
outFile=$dir/ukbb.${phenoName}.${CHR}.${suffix2}.clumped.clumped
tail -n +2 $outFile | head -n -2 >> $merged_outFile
done

# List of SNPs...
dir=$outdir
merged_outFile=$dir/ukbb.${phenoName}.${suffix2}.clumped.clumped.combined.txt 
merged_outFile2=$dir/ukbb.${phenoName}.${suffix2}.clumped.clumped.combined.cut.txt
awk '{print $1,$3}' $merged_outFile | tail -n +2 > $merged_outFile2

# list of snps part 2...
dir=$outdir
outdir=$outdir # /athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_GxG
mkdir -p $outdir
merged_outFile2=$dir/ukbb.${phenoName}.${suffix2}.clumped.clumped.combined.cut.txt
merged_outFile3=$dir/ukbb.${phenoName}.${suffix2}.clumped.clumped.combined.cut.2.txt
echo $phenoName > $merged_outFile3
# cut -f2 -d' ' $dir/$merged_outFile2 >> $outdir/$merged_outFile3
cut -f2 -d' ' $merged_outFile2 >> $merged_outFile3
echo "END" >> $merged_outFile3
