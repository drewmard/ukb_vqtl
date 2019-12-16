#!/bin/bash -l
#SBATCH -J clump
#SBATCH --mem=128G
#SBATCH --tasks-per-node=1
#SBATCH -N 1
#SBATCH -n 1

# initialize
phenoName=lymphocyte.count.rint.ALL
# thres="5e-8"
thres="0.001"
genodir=/athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/andrew_copies

# make important directories
# mkdir -p ${indir}/output/GWAS/subset
# mkdir -p ${indir}/output/GWAS/subset/$phenoName
# mkdir -p $indir/output/GWAS/results2/
# mkdir -p $indir/output/GWAS/results2/$phenoName

# clump
echo 'Clumping!'
source activate vQTL
outdir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_clump # /${phenoName}
results=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_subset/ukbb.${phenoName}.vGWAS.txt
for CHR in {1..22}
do
echo ${CHR}
geno=$genodir/ukbb.${CHR}.impute
outFile=ukbb.${phenoName}.${CHR}.results.var
plink --bfile $geno --maf 0.05 --clump $results --clump-p1 $thres --clump-p2 $thres --clump-r2 0.8 --clump-field P --clump-snp-field rs --out $outdir/$outFile

done


# Merge clump 

echo 'Merging!'
# echo "Doing the means..."
# merged_outFile=ukbb.${phenoName}.results.mean.clumped.txt
# head -1 ukbb.${phenoName}.1.results.mean.clumped > $merged_outFile
# for CHR in {1..22}
# do
# echo $CHR
# geno=ukbb.${CHR}
# outFile=ukbb.${phenoName}.${CHR}.results.mean.clumped
# tail -n +2 $indir/$outFile | head -n -2 >> $indir/$merged_outFile
# done

echo "Doing the variances..."
dir=$indir/output/GWAS/results2/${phenoName}
merged_outFile=$dir/ukbb.${phenoName}.results.var.clumped.txt
head -1 $dir/ukbb.${phenoName}.1.results.var.clumped > $merged_outFile
for CHR in {1..22}
do
echo $CHR
geno=ukbb.${CHR}
outFile=$dir/ukbb.${phenoName}.${CHR}.results.var.clumped
tail -n +2 $outFile | head -n -2 >> $merged_outFile
done

# ok if error: may be no significant hits


# pull out SNPs into list

# echo "List of mean SNPs..."
# line=1
# merged_outFile=ukbb.${phenoName}.results.mean.clumped.txt
# merged_outFile2=ukbb.${phenoName}.results.mean.clumped.cut.txt
# awk '{print $1,$3}' $indir/$merged_outFile | tail -n +2 > $indir/$merged_outFile2

echo "List of var SNPs..."
# line=1
dir=$indir/output/GWAS/results2/${phenoName}
merged_outFile=ukbb.${phenoName}.results.var.clumped.txt
merged_outFile2=ukbb.${phenoName}.results.var.clumped.cut.txt
awk '{print $1,$3}' $dir/$merged_outFile | tail -n +2 > $dir/$merged_outFile2


# create raw genotype files

echo "Extract SNP time."
dir=$indir/output/GWAS/results2/${phenoName}
outdir=$indir/output/GWAS/subset/$phenoName
merged_outFile2=ukbb.${phenoName}.results.var.clumped.cut.txt
mkdir -p $outdir
while IFS=$' ' read -r -a myArray
do
 CHR=${myArray[0]}
 snp=${myArray[1]}
 prefix=ukbb.$CHR
 plink --bfile $genodir/$prefix --snp $snp --recodeA --out $outdir/$prefix.$snp
done < $dir/$merged_outFile2 # if running pipeline

