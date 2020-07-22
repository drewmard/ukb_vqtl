source activate vQTL
pheno=bmi
SCRIPTDIR=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/scripts/FINAL/gwas/run

#############

# vGWAS on raw
phenoName=${pheno}.ALL
sbatch $SCRIPTDIR/run_vGWAS_subset.sh $phenoName

# vGWAS on rint
phenoName=${pheno}.rint.ALL
sbatch $SCRIPTDIR/run_vGWAS_subset.sh $phenoName

# muGWAS on raw
phenoName=${pheno}.ALL
sbatch $SCRIPTDIR/run_GWAS.impute.sh $phenoName

# muGWAS on rint
phenoName=${pheno}.rint.ALL
sbatch $SCRIPTDIR/run_GWAS.impute.sh $phenoName

# Run HLMM/DET
# see ukb_vqtl/scripts/FINAL/gwas/HLMM/README.md

###############

# Merge together results
# note: if error on the first run, will RE-RUN!!!

# vGWAS on raw
phenoName=${pheno}.ALL
$SCRIPTDIR/merge_vGWAS_subset.sh $phenoName
mv /athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_subset/ukbb.$phenoName.vGWAS.txt /athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_subset/ukbb.$phenoName.vGWAS.old.txt # simple rename & re-adjusting datasets for vGWAS
Rscript $SCRIPTDIR/merge_vGWAS_subset_2.R $phenoName

# vGWAS on rint
phenoName=${pheno}.rint.ALL
$SCRIPTDIR/merge_vGWAS_subset.sh $phenoName
mv /athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_subset/ukbb.$phenoName.vGWAS.txt /athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_subset/ukbb.$phenoName.vGWAS.old.txt # simple rename & re-adjusting datasets for vGWAS
Rscript $SCRIPTDIR/merge_vGWAS_subset_2.R $phenoName

# muGWAS on raw
phenoName=${pheno}.ALL
Rscript /athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/scripts/GWAS/mergeResults_impute.R $phenoName

# muGWAS on rint
phenoName=${pheno}.rint.ALL
Rscript /athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/scripts/GWAS/mergeResults_impute.R $phenoName

# HLMM/DET merge: see ukb_vqtl/scripts/FINAL/gwas/HLMM/README.md

#################

source activate vQTL

# a: trim vGWAS on RINT
pheno=$pheno
suffix1="rint.ALL"
suffix2="var"
# thres="0.001" # old
thres="1e-5"
phenoName=${pheno}.${suffix1}

/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/FINAL/other/GWAS_to_clumped_SNP_list.sh $pheno $suffix1 $suffix2 $thres

# b: trim vGWAS on raw
pheno=$pheno
suffix1="ALL"
suffix2="var"
# thres="1e-5" # old
thres="5e-8"
phenoName=${pheno}.${suffix1}

/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/FINAL/other/GWAS_to_clumped_SNP_list.sh $pheno $suffix1 $suffix2 $thres


# c: trim vGWAS on mean-based RINT
### need to adjust for mean-based GWAS files
# pheno=$pheno
# suffix1="rint.ALL"
# suffix2="mean"
# thres="5e-8"
# phenoName=${pheno}.${suffix1}
# 
# /athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/FINAL/other/GWAS_to_clumped_SNP_list.sh $pheno $suffix1 $suffix2 $thres

# d: trim vGWAS on mean-based raw
### need to adjust for mean-based GWAS files
pheno=$pheno
suffix1="ALL"
suffix2="mean"
thres="5e-8"
phenoName=${pheno}.${suffix1}

/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/FINAL/other/GWAS_to_clumped_SNP_list.sh $pheno $suffix1 $suffix2 $thres


# e: trim dGWAS on dispersion-based rint
### need to adjust for mean-based GWAS files
pheno=$pheno
suffix1="rint.ALL"
suffix2="DET"
thres="1e-5"
phenoName=${pheno}.${suffix1}

/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/FINAL/other/GWAS_to_clumped_SNP_list.sh $pheno $suffix1 $suffix2 $thres


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

suffix1="rint.ALL"
suffix2="DET"
phenoName=${pheno}.${suffix1}
dir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/dGWAS_clump
file4=$dir/ukbb.${phenoName}.${suffix2}.clumped.clumped.combined.cut.2.txt

fileMerged=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.${pheno}.SNP_list.txt
echo $pheno > $fileMerged
sed '1d;$d' $file1 >> $fileMerged
sed '1d;$d' $file2 >> $fileMerged
sed '1d;$d' $file3 >> $fileMerged
sed '1d;$d' $file4 >> $fileMerged
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
###################
# merge subsets together

genoSub_next=$outdir/ukbb.${pheno}.merged_subset
plink --bfile $outdir/ukbb.${pheno}.1 --bmerge $outdir/ukbb.${pheno}.2 --make-bed --out $genoSub_next
if [ ! -f $genoSub_next.bed ]; then
	cp $outdir/ukbb.${pheno}.1.bed ${genoSub_next}.bed
	cp $outdir/ukbb.${pheno}.1.fam ${genoSub_next}.fam
	cp $outdir/ukbb.${pheno}.1.bim ${genoSub_next}.bim
fi

for CHR in {3..22}
do
echo ${CHR}
geno=$outdir/ukbb.${pheno}.merged_subset
genoSub=$outdir/ukbb.${pheno}.${CHR}
if [ -f $genoSub.fam ]; then
plink --bfile $genoSub --bmerge $geno --make-bed --out $geno
fi
done


##################################

snp_file=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/sig_results/${pheno}.SNP_sig_list.txt
outdir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2

geno=$outdir/ukbb.${pheno}.merged_subset
outFile=$outdir/ukbb.${pheno}.merged_subset2
plink --bfile $geno --extract $snp_file --make-bed --out $outFile



###################################



# 4:
pheno=bmi
outdir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2
prefix=ukbb.${pheno}.merged_subset # use before LD trimming
prefix=ukbb.${pheno}.merged_subset2
merged_outFile3=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GxG_2/ukbb.${pheno}.SNP_list.txt
phenoFile=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/phenotypes_processed.80.txt
phenoName=$pheno.ALL
plink --bfile $outdir/$prefix --pheno $phenoFile --pheno-name $phenoName --epistasis set-by-set --set ${merged_outFile3} --epi1 1 --allow-no-sex --out $outdir/$prefix.GxG.80
plink --bfile $outdir/$prefix --r2 dprime yes-really --ld-window-kb 50000 --ld-window-r2 0 --ld-window 1000 --out $outdir/$prefix.LD

pheno=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/GWAS/preprocess/phenotypes_processed.20.txt
plink --bfile $outdir/$prefix --pheno $pheno --pheno-name $phenoName --epistasis set-by-set --set ${merged_outFile3} --epi1 1 --allow-no-sex --out $outdir/$prefix.GxG.20

