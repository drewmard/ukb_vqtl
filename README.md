# vQTL analysis in UK Biobank

## Linear Modeling of Deviations
Directory: ./scripts/L-MOD/


#### Script 1
< insert script 1 >

#### Script 1
< insert script 1 > \
< insert >



## UKB GWAS of blood cell phenotypes
Directory: ./scripts/GWAS/


#### pre-processing:
(1) identify_indiv_blood_disorders.R \
(2) ukb_sample_qc.txt \

#### run GWAS
run_GWAS.sh \
(this script runs run_GWAS.slurm.sh)

#### merge diff assoc files together:
mergeResults.R

#### clump to significant hits, identify top snps, extract from genotype file
clump_identify_extract.sh



## UKB GxE testing (post-vQTL screening)
Directory: _______


#### generate environmental factors
< insert script 1 >

#### Run GxE
GxE.R


