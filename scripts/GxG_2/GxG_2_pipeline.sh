source activate vQTL
pheno=lymphocyte.count
pheno=bmi

#############

phenoName=${pheno}.rint.ALL
SCRIPTDIR=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/scripts/GWAS/subset

#########
# vGWAS
sbatch $SCRIPTDIR/run_vGWAS_subset.sh $phenoName

# Merge together results
# note: if error on the first run, will RE-RUN!!!
./merge_vGWAS_subset.sh $phenoName

# simple rename & re-adjusting datasets for vGWAS
mv /athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_subset/ukbb.$phenoName.vGWAS.txt /athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_subset/ukbb.$phenoName.vGWAS.old.txt
Rscript merge_vGWAS_subset_2.R $phenoName


#################

source activate HLMM

sbatch $SCRIPTDIR/run_HLMM_subset.sh $phenoName

$SCRIPTDIR/merge_HLMM_subset.sh $phenoName # this does merge & re-run for HLMM


#################
source activate vQTL

# a: trim vGWAS on RINT
pheno=$pheno
suffix1="rint.ALL"
suffix2="var"
thres="0.001"
phenoName=${pheno}.${suffix1}

/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/scripts/GxG_2/GWAS_to_clumped_SNP_list.sh $pheno $suffix1 $suffix2 $thres

# b: trim vGWAS on raw
pheno=$pheno
suffix1="ALL"
suffix2="var"
thres="1e-5"
phenoName=${pheno}.${suffix1}

/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/scripts/GxG_2/GWAS_to_clumped_SNP_list.sh $pheno $suffix1 $suffix2 $thres


# c: trim vGWAS on mean-based RINT
### need to adjust for mean-based GWAS files
pheno=$pheno
suffix1="rint.ALL"
suffix2="mean"
thres="5e-8"
phenoName=${pheno}.${suffix1}

/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/scripts/GxG_2/GWAS_to_clumped_SNP_list.sh $pheno $suffix1 $suffix2 $thres

# d: trim vGWAS on mean-based raw
### need to adjust for mean-based GWAS files
pheno=$pheno
suffix1="ALL"
suffix2="mean"
thres="5e-8"
phenoName=${pheno}.${suffix1}

/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/scripts/GxG_2/GWAS_to_clumped_SNP_list.sh $pheno $suffix1 $suffix2 $thres

######################

# merge into SNPs into single file...
suffix1="rint.ALL"
suffix2="var"
phenoName=${pheno}.${suffix1}
if [ "$suffix2" == "var" ]; then
  dir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_clump # /${phenoName}
elif [ "$suffix2" == "mean" ]; then
  dir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/muGWAS_clump
fi
file1=$dir/ukbb.${phenoName}.${suffix2}.clumped.clumped.combined.cut.2.txt

suffix1="ALL"
suffix2="var"
phenoName=${pheno}.${suffix1}
if [ "$suffix2" == "var" ]; then
  dir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_clump # /${phenoName}
elif [ "$suffix2" == "mean" ]; then
  dir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/muGWAS_clump
fi
file2=$dir/ukbb.${phenoName}.${suffix2}.clumped.clumped.combined.cut.2.txt

suffix1="ALL"
suffix2="mean"
phenoName=${pheno}.${suffix1}
if [ "$suffix2" == "var" ]; then
  dir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_clump # /${phenoName}
elif [ "$suffix2" == "mean" ]; then
  dir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/muGWAS_clump
fi
file3=$dir/ukbb.${phenoName}.${suffix2}.clumped.clumped.combined.cut.2.txt

fileMerged=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.${pheno}.SNP_list.txt
echo $pheno > $fileMerged
sed '1d;$d' $file1 >> $fileMerged
sed '1d;$d' $file2 >> $fileMerged
sed '1d;$d' $file3 >> $fileMerged
echo "END" >> $fileMerged

###########################

# 2: extract SNPs from files

snp_file=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.${pheno}.SNP_list.txt
outdir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2
genodir=/athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/andrew_copies
for CHR in {1..22}
do
echo ${CHR}
geno=$genodir/ukbb.${CHR}.impute
outFile=ukbb.${pheno}.${CHR}
plink --bfile $geno --extract $snp_file --make-bed --out $outdir/$outFile
done


########################

# 3: merge different CHR subsets into 1 file

genoSub_next=$outdir/ukbb.${pheno}.2.sub
plink --bfile $outdir/ukbb.${pheno}.1 --bmerge $outdir/ukbb.${pheno}.2 --make-bed --out $genoSub_next
if [ ! -f $genoSub_next ]; then
	cp $outdir/ukbb.${pheno}.1.bed ${genoSub_next}.bed
	cp $outdir/ukbb.${pheno}.1.fam ${genoSub_next}.fam
	cp $outdir/ukbb.${pheno}.1.bim ${genoSub_next}.bim
	# cp $genoSub_prev $genoSub_next
fi

for CHR in {3..22}
do
echo ${CHR}
geno=$outdir/ukbb.${pheno}.${CHR}
genoSub_prev=$outdir/ukbb.${pheno}.$(($CHR-1)).sub
genoSub_next=$outdir/ukbb.${pheno}.${CHR}.sub
if [ -f $geno.fam ]; then
plink --bfile $genoSub_prev --bmerge $geno --make-bed --out $genoSub_next
else 
cp $genoSub_prev.bed ${genoSub_next}.bed
cp $genoSub_prev.fam ${genoSub_next}.fam
cp $genoSub_prev.bim ${genoSub_next}.bim
fi
done

mv $genoSub_next.bim $outdir/ukbb.${pheno}.ALL.sub.bim
mv $genoSub_next.fam $outdir/ukbb.${pheno}.ALL.sub.fam
mv $genoSub_next.bed $outdir/ukbb.${pheno}.ALL.sub.bed

for CHR in {2..22}
do
genoSub_next=$outdir/ukbb.${pheno}.${CHR}.sub
rm ${genoSub_next}*
done

# 4:
# pheno=lymphocyte.count

outdir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2
prefix=ukbb.${pheno}.ALL.sub
merged_outFile3=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.${pheno}.SNP_list.txt
pheno=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/phenotypes_processed.80.txt
phenoName=lymphocyte.count.ALL
plink --bfile $outdir/$prefix --pheno $pheno --pheno-name $phenoName --epistasis set-by-set --set ${merged_outFile3} --epi1 1 --allow-no-sex --out $outdir/$prefix.GxG.80
plink --bfile $outdir/$prefix --r2 dprime yes-really --ld-window-kb 50000 --ld-window-r2 0 --ld-window 1000 --out $outdir/$prefix.LD
wc -l $outdir/$prefix.LD.ld

plink --bfile $outdir/$prefix --ld 6:31242257_AC_A rs2394982
--ld-window-r2 0.1 --r2 dprime


pheno=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/phenotypes_processed.20.txt
plink --bfile $outdir/$prefix --pheno $pheno --pheno-name $phenoName --epistasis set-by-set --set ${merged_outFile3} --epi1 1 --allow-no-sex --out $outdir/$prefix.GxG.20

plink --bfile $outdir/$prefix --snps 19:16533781_AG_A,rs73511187,rs12461448,rs61420814 --recodeA


plink --bfile $outdir/$prefix --freq --out $outdir/$prefix.GxG


# 2: prune
# 250 kb
# r < 0.5



