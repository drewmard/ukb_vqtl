source activate vQTL
pheno=lymphocyte.count

#############
# vGWAS!

phenoName=${pheno}.rint.ALL
# Run vQTL screen
sbatch run_vGWAS_subset.sh $phenoName

# Merge together results
# note: if error on the first run, will RE-RUN!!!
./merge_vGWAS_subset.sh $phenoName

# simple rename & re-adjusting datasets
mv /athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_subset/ukbb.$phenoName.vGWAS.txt /athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_subset/ukbb.$phenoName.vGWAS.old.txt
Rscript merge_vGWAS_subset_2.R $phenoName

#################
source activate vQTL

# a: trim vGWAS on RINT

# input:
pheno=$pheno
suffix1="rint.ALL"
suffix2="var"
thres="0.001"
phenoName=${pheno}.${suffix1}

/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/scripts/GxG_2/GWAS_to_clumped_SNP_list.sh $pheno $suffix1 $suffix2 $thres

# b: trim vGWAS on raw
thres="1e-5"
phenoName=${pheno}.ALL
outdir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_clump # /${phenoName}
results=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_subset/ukbb.${phenoName}.vGWAS.txt

echo 'Clumping!'
for CHR in {1..22}
do
echo ${CHR}
geno=$genodir/ukbb.${CHR}.impute
outFile=ukbb.${phenoName}.${CHR}.results.var
plink --bfile $geno --maf 0.05 --clump $results --clump-p1 $thres --clump-p2 $thres --clump-r2 0.8 --clump-field P --clump-snp-field rs --out $outdir/$outFile
done


# c: trim vGWAS on mean-based RINT
### need to adjust for mean-based GWAS files
thres="5e-8"
phenoName=${pheno}.rint.ALL
outdir=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_clump # /${phenoName}
results=/athena/elementolab/scratch/anm2868/vQTL/ukb_vqtl/output/vGWAS_subset/ukbb.${phenoName}.vGWAS.txt

echo 'Clumping!'
for CHR in {1..22}
do
echo ${CHR}
geno=$genodir/ukbb.${CHR}.impute
outFile=ukbb.${phenoName}.${CHR}.results.var
plink --bfile $geno --maf 0.05 --clump $results --clump-p1 $thres --clump-p2 $thres --clump-r2 0.8 --clump-field P --clump-snp-field rs --out $outdir/$outFile
done


# 2: prune
# 250 kb
# r < 0.5



