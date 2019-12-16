# vQTL analysis in UK Biobank

## UKB GWAS of blood cell phenotypes
Directory: ./scripts/GWAS/


#### pre-processing:
split_80_20.R \
identify_indiv_blood_disorders.R \
sample_qc.R \
gen_cov1.R \
gen_cov2.R \
gen_pheno.R \
gen_full_data.R \
preprocess.R

#### Script to perform L-MOD variance testing in Plink
vGWAS.R

#### run GWAS
This is a standard GWAS analysis, where SNPs are tested one-by-one for association with the marginal means of a quantitative phenotype.

run_GWAS.sh (run on cannes; local)

### run vGWAS
In the vGWAS setting, SNPs are tested one-by-one for association with the marginal variance of a quantitative phenotype. In our testing framework, we developed a Deviation Regression Model for variance testing, which is a linear model modification to the  Levene's and Brown-Forsythe's classical tests. The DRM assesses whether additive differences in the variance exist across genotypes. The test is also robust to any genetic effects on the mean of the phenotype.

#### run vGWAS using plink's R plug-in
array job: sbatch run_vGWAS_5.sh <phenotype> (run on cannes locally/curie.pbtech cluster for all chr) \
(sbatch) run_vGWAS_4.sh <phenotype> <CHR> (run on cannes locally/curie.pbtech cluster for spec chr) \
run_vGWAS_3.sh <phenotype> (run on cannes locally in a loop) \

#### run vGWAS using R's BEDMatrix library



#### merge diff assoc files together:
mergeResults.R ########## note: might need to change directories for mean gwas analyses
mergeTransform.R

#### clump to significant hits, identify top snps, extract from genotype file
clump_identify_extract.sh



## UKB GxE testing (post-vQTL screening)
Directory: ./scripts/GxE


#### generate environmental factors
gen_envir_factors.R

#### generate genotype files
gen_geno.R (need to run clump_identify_extract.sh first)

#### Run GxE
GxE.R (might need to clean up?) \
(need to switch var hits file to be automated) \
(switch to calculating residuals, identical to preprocess; maybe just extract preprocess) \
(need to make it look cleaner)


## UKB downstream testing (post-GxE)
Directory: ./scripts/downstream

#### new phenotypes
gen_pheno_other.R

#### do PGS calculations
create_PGS.R # for creating PGS file for scott to run \
PGS.R


