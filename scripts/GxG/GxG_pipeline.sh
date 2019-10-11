
# 1: form list of SNPs
# phenoName=lymphocyte.count.rint.ALL
phenoName=monocyte.count.rint.ALL
dir=$indir/output/GWAS/results2/${phenoName}
outdir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG
merged_outFile2=ukbb.${phenoName}.results.var.clumped.cut.txt
merged_outFile3=ukbb.${phenoName}.results.var.clumped.cut.2.txt
echo $phenoName > $outdir/$merged_outFile3
cut -f2 -d' ' $dir/$merged_outFile2 >> $outdir/$merged_outFile3
echo "END" >> $outdir/$merged_outFile3
# 2: extract SNPs from files

for CHR in {1..22}
do
echo ${CHR}
geno=$genodir/ukbb.${CHR}
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
plink --bfile $genoSub_prev --bmerge $geno --make-bed --out $genoSub_next
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
outdir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG
prefix=ukbb.${phenoName}.ALL.sub
merged_outFile3=$outdir/ukbb.${phenoName}.results.var.clumped.cut.2.txt
pheno=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/phenotypes_processed.80.txt
plink --bfile $outdir/$prefix --pheno $pheno --pheno-name $phenoName --epistasis set-by-set --set ${merged_outFile3} --epi1 1 --allow-no-sex --out $outdir/$prefix.GxG
