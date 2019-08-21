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
preprocess.R \

#### Script to perform L-MOD variance testing in Plink
vGWAS.R

#### run GWAS & vGWAS
run_GWAS.sh (run on cannes; local) \
run_vGWAS.sh (run on curie.pbtech SLURM cluster) 

#### merge diff assoc files together:
mergeResults.R ########## note: might need to change directories for mean gwas analyses

#### clump to significant hits, identify top snps, extract from genotype file
clump_identify_extract.sh



## UKB GxE testing (post-vQTL screening)
Directory: _______


#### generate environmental factors
< insert script 1 >

#### Run GxE
GxE.R


