# phenoName=lymphocyte.count.rint.ALL
# phenoName=monocyte.count.rint.ALL
# phenoName=neutrophil.count.rint.ALL
# phenoName=wbc.leukocyte.count.rint.ALL
# phenoName=rbc.erythrocyte.count.rint.ALL
# phenoName=monocyte.count.ALL
phenoName=lymphocyte.count.ALL
phenoName=eosinophil.count.rint.ALL

# Run vQTL screen
sbatch run_vGWAS_subset.sh $phenoName

##########################################################

# Merge together results
./merge_vGWAS_subset.sh $phenoName
mv /athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_subset/ukbb.$phenoName.vGWAS.txt /athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_subset/ukbb.$phenoName.vGWAS.old.txt
Rscript merge_vGWAS_subset_2.R $phenoName

# Rscript merge_vGWAS_subset.R $phenoName

##########################################################

thres="0.001"
genodir=/athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/andrew_copies

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

echo "Doing the variances..."
dir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_clump
merged_outFile=$dir/ukbb.${phenoName}.results.var.clumped.txt
head -1 $dir/ukbb.${phenoName}.${CHR}.results.var.clumped > $merged_outFile
for CHR in {1..22}
do
echo $CHR
geno=ukbb.${CHR}
outFile=$dir/ukbb.${phenoName}.${CHR}.results.var.clumped
tail -n +2 $outFile | head -n -2 >> $merged_outFile
done

echo "List of var SNPs..."
# line=1
dir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_clump/
merged_outFile=ukbb.${phenoName}.results.var.clumped.txt
merged_outFile2=ukbb.${phenoName}.results.var.clumped.cut.txt
awk '{print $1,$3}' $dir/$merged_outFile | tail -n +2 > $dir/$merged_outFile2


##########################################################

# 1: form list of SNPs
dir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_clump
outdir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_GxG
mkdir -p $outdir
merged_outFile2=ukbb.${phenoName}.results.var.clumped.cut.txt
merged_outFile3=ukbb.${phenoName}.results.var.clumped.cut.2.txt
echo $phenoName > $outdir/$merged_outFile3
cut -f2 -d' ' $dir/$merged_outFile2 >> $outdir/$merged_outFile3
echo "END" >> $outdir/$merged_outFile3

# 2: extract SNPs from files
genodir=/athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/andrew_copies

for CHR in {1..22}
do
echo ${CHR}
geno=$genodir/ukbb.${CHR}.impute
outFile=ukbb.${phenoName}.${CHR}
plink --bfile $geno --extract $outdir/$merged_outFile3 --make-bed --out $outdir/$outFile
done

# 3: merge different CHR subsets into 1 file

genoSub_next=$outdir/ukbb.${phenoName}.2.sub
plink --bfile $outdir/ukbb.${phenoName}.1 --bmerge $outdir/ukbb.${phenoName}.2 --make-bed --out $genoSub_next
if [ ! -f $genoSub_next ]; then
	cp $outdir/ukbb.${phenoName}.1.bed ${genoSub_next}.bed
	cp $outdir/ukbb.${phenoName}.1.fam ${genoSub_next}.fam
	cp $outdir/ukbb.${phenoName}.1.bim ${genoSub_next}.bim
	# cp $genoSub_prev $genoSub_next
fi

for CHR in {3..22}
do
echo ${CHR}
geno=$outdir/ukbb.${phenoName}.${CHR}
genoSub_prev=$outdir/ukbb.${phenoName}.$(($CHR-1)).sub
genoSub_next=$outdir/ukbb.${phenoName}.${CHR}.sub
if [ -f $geno.fam ]; then
plink --bfile $genoSub_prev --bmerge $geno --make-bed --out $genoSub_next
else 
cp $genoSub_prev.bed ${genoSub_next}.bed
cp $genoSub_prev.fam ${genoSub_next}.fam
cp $genoSub_prev.bim ${genoSub_next}.bim
fi
done

mv $genoSub_next.bim $outdir/ukbb.${phenoName}.ALL.sub.bim
mv $genoSub_next.fam $outdir/ukbb.${phenoName}.ALL.sub.fam
mv $genoSub_next.bed $outdir/ukbb.${phenoName}.ALL.sub.bed

for CHR in {2..22}
do
genoSub_next=$outdir/ukbb.${phenoName}.${CHR}.sub
rm ${genoSub_next}*
done

# 4:
outdir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_GxG
prefix=ukbb.${phenoName}.ALL.sub
merged_outFile3=$outdir/ukbb.${phenoName}.results.var.clumped.cut.2.txt
pheno=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/phenotypes_processed.80.txt
plink --bfile $outdir/$prefix --pheno $pheno --pheno-name $phenoName --epistasis set-by-set --set ${merged_outFile3} --epi1 1 --allow-no-sex --out $outdir/$prefix.GxG.80

pheno=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/phenotypes_processed.20.txt
plink --bfile $outdir/$prefix --pheno $pheno --pheno-name $phenoName --epistasis set-by-set --set ${merged_outFile3} --epi1 1 --allow-no-sex --out $outdir/$prefix.GxG.20




