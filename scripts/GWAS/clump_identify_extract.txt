# clump

indir=/athena/elementolab/scratch/anm2868/vQTL/UKB/results
phenoName=lymphocyte.count
results=$indir/ukbb.${phenoName}.results.txt
genodir=/home/kulmsc/athena/ukbiobank/calls
mkdir -p /athena/elementolab/scratch/anm2868/vQTL/UKB/subset/phenoName
for CHR in {1..22}
do
echo ${CHR}
geno=ukbb.${CHR}
outFile=ukbb.${phenoName}.${CHR}.results.var
plink --bfile $genodir/$geno --maf 0.05 --clump $results --clump-p1 5e-8 --clump-p2 5e-8 --clump-r2 0.01 --clump-kb 5000 --clump-field P.y --out $indir/$outFile
# outFile=ukbb.${phenoName}.${CHR}.results.mean
# plink --bfile $genodir/$geno --maf 0.05 --clump $results --clump-p1 5e-8 --clump-p2 5e-8 --clump-r2 0.01 --clump-kb 5000 --clump-field P.x --out $indir/$outFile
done


# Merge clump 

merged_outFile=ukbb.${phenoName}.results.mean.clumped.txt
head -1 ukbb.${phenoName}.1.results.mean.clumped > $merged_outFile
for CHR in {1..22}
do
geno=ukbb.${CHR}
outFile=ukbb.${phenoName}.${CHR}.results.mean.clumped
tail -n +2 $indir/$outFile | head -n -2 >> $indir/$merged_outFile
done

merged_outFile=ukbb.${phenoName}.results.var.clumped.txt
head -1 ukbb.${phenoName}.1.results.var.clumped > $merged_outFile
for CHR in {1..22}
do
geno=ukbb.${CHR}
outFile=ukbb.${phenoName}.${CHR}.results.var.clumped
tail -n +2 $indir/$outFile | head -n -2 >> $indir/$merged_outFile
done

# ok if error: may be no significant hits


# pull out SNPs into list

line=1
merged_outFile=ukbb.${phenoName}.results.mean.clumped.txt
merged_outFile2=ukbb.${phenoName}.results.mean.clumped.cut.txt
awk '{print $1,$3}' $indir/$merged_outFile | tail -n +2 > $indir/$merged_outFile2

line=1
merged_outFile=ukbb.${phenoName}.results.var.clumped.txt
merged_outFile2=ukbb.${phenoName}.results.var.clumped.cut.txt
awk '{print $1,$3}' $indir/$merged_outFile | tail -n +2 > $indir/$merged_outFile2


# create raw genotype files

dir=/home/kulmsc/athena/ukbiobank/calls
outdir=/athena/elementolab/scratch/anm2868/vQTL/UKB/subset/$phenoName
merged_outFile2=ukbb.${phenoName}.results.mean.clumped.cut.txt
mkdir -p $outdir
while IFS=$' ' read -r -a myArray
do
 CHR=${myArray[0]}
 snp=${myArray[1]}
 prefix=ukbb.$CHR
 plink --bfile $dir/$prefix --snp $snp --recodeA --out $outdir/$prefix.$snp
done < $indir/$merged_outFile2 # if running pipeline
