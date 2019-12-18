# Trim SNPs #1:
# Exclude MAF < 0.05
# Exclude Imputation INFO score < 0.3
# trim dataset using bgenix
# then convert to plink BED format

# (1)
# Identify SNPs:
for CHR in {1..22}; do
  Rscript ./scripts1/variants1.R $CHR
done

# (2)
# Use bgenix to trim:
for CHR in {1..22}; do
  for j in {1..20}; do
  echo ${CHR}
  inFile=/home/kulmsc/athena/ukbiobank/imputed/ukbb.$CHR.bgen
  rsid=/athena/elementolab/scratch/anm2868/vQTL/UKB/imputed/variants1/ukb.$CHR.imputed.variants1.txt
  outFile=/athena/elementolab/scratch/anm2868/vQTL/UKB/imputed/bgen_variants1/ukbb.$CHR.variants1.bgen
  rsid=/athena/elementolab/scratch/anm2868/vQTL/UKB/imputed/variants1/ukb.$CHR.imputed.variants1.$j.txt
  outFile=/athena/elementolab/scratch/anm2868/vQTL/UKB/imputed/bgen_variants1/ukbb.$CHR.variants1.$j.bgen

  /home/kulmsc/bin/bgenix -g $inFile -incl-rsids $rsid > $outFile
  done
done

# (3)
# Convert to plink BED using plink2 and encode hard-call genotypes (--hard-call 0.1)
for CHR in {1..22}; do
  for j in {1..20}; do

  echo ${CHR}
  inBGEN=/athena/elementolab/scratch/anm2868/vQTL/UKB/imputed/bgen_variants1/ukbb.$CHR.variants1.bgen
  inSAMPLE=/home/kulmsc/athena/ukbiobank/imputed/ukbb.${CHR}.sample
  outFile=/athena/elementolab/scratch/anm2868/vQTL/UKB/imputed/plink_variants1/ukbb.$CHR.variants1
  
  inBGEN=/athena/elementolab/scratch/anm2868/vQTL/UKB/imputed/bgen_variants1/ukbb.$CHR.variants1.$j.bgen
  inSAMPLE=/home/kulmsc/athena/ukbiobank/imputed/ukbb.${CHR}.sample
  outFile=/athena/elementolab/scratch/anm2868/vQTL/UKB/imputed/plink_variants1/ukbb.$CHR.variants1.$j

  plink2 --bgen $inBGEN --sample $inSAMPLE --make-bed --hard-call-threshold 0.1 --out $outFile

  # plink2 --bgen $inBGEN --sample $inSAMPLE --make-bed --hard-call-threshold 0.1 --out $outFile
  done

done

# Trim SNPs #2:
# Exclude HWE test p-values < 10^-16
# Exclude missing genotype rate > 0.05

for CHR in {1..22}
do
  for j in {1..20}; do
  echo ${CHR}
  geno=/athena/elementolab/scratch/anm2868/vQTL/UKB/imputed/plink_variants1/ukbb.$CHR.variants1
  genoNew=/athena/elementolab/scratch/anm2868/vQTL/UKB/imputed/plink_variants2/ukbb.$CHR.variants2
  
  geno=/athena/elementolab/scratch/anm2868/vQTL/UKB/imputed/plink_variants1/ukbb.$CHR.variants1.$j
  genoNew=/athena/elementolab/scratch/anm2868/vQTL/UKB/imputed/plink_variants2/ukbb.$CHR.variants2.$j

  plink --bfile $geno --hwe 1e-16 --geno 0.05 --make-bed --out $genoNew
  done
done


