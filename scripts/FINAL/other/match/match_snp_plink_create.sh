
# 2: extract SNPs from files
QTLtype=QTL

pheno=bmi
snp_file=/athena/elementolab/scratch/anm2868/open_targets/query/bmi.$QTLtype.matched_snp.txt
outdir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2
genodir=/athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/andrew_copies
for CHR in {1..22}
do
echo ${CHR}
geno=$genodir/ukbb.${CHR}.impute
outFile=ukbb.${pheno}.${CHR}.${QTLtype}_matched_snp
plink --bfile $geno --extract $snp_file --make-bed --out $outdir/$outFile
done


########################

# 3: merge different CHR subsets into 1 file
###################
# merge subsets together

rm $outdir/ukbb.${pheno}.${QTLtype}_matched_snp.merged_subset.{bim,bed,fam}
genoSub_next=$outdir/ukbb.${pheno}.${QTLtype}_matched_snp.merged_subset
plink --bfile $outdir/ukbb.${pheno}.1.${QTLtype}_matched_snp --bmerge $outdir/ukbb.${pheno}.2.${QTLtype}_matched_snp --make-bed --out $genoSub_next
# if [ ! -f $genoSub_next.bed ]; then
# 	cp $outdir/ukbb.${pheno}.1.${QTLtype}_matched_snp.bed ${genoSub_next}.bed
# 	cp $outdir/ukbb.${pheno}.1.${QTLtype}_matched_snp.fam ${genoSub_next}.fam
# 	cp $outdir/ukbb.${pheno}.1.${QTLtype}_matched_snp.bim ${genoSub_next}.bim
# fi

for CHR in {3..22}
do
echo ${CHR}
geno=$outdir/ukbb.${pheno}.${QTLtype}_matched_snp.merged_subset
genoSub=$outdir/ukbb.${pheno}.${CHR}.${QTLtype}_matched_snp
if [ -f $genoSub.fam ]; then
plink --bfile $genoSub --bmerge $geno --make-bed --out $geno
fi
done
