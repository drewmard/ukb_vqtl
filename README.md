# vQTL analysis in UK Biobank

## Linear Modeling of Deviations (L-MOD)
Directory: ./scripts/L-MOD/


#### Script 1
< insert script 1 >

< insert script 1 > \
< insert >



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
run_GWAS.sh (run on cannes; local) \

#### run vGWAS
(sbatch) run_vGWAS_4.sh <phenotype> <CHR> (run on cannes locally/curie.pbtech cluster for spec chr) 
run_vGWAS_3.sh <phenotype> (run on cannes locally in a loop) 

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
PGS.R


