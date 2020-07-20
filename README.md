# Leveraging phenotypic variability to identify genetic interactions in human phenotypes

This repository contains all code to conduct the analyses present in the manuscript, "Leveraging phenotypic variability to identify genetic interactions in human phenotypes". This research describes a statistical framework to find SNPs associated with the means and variances of a quantitative phenotype, and then using these SNPs to search for gene-environment interactions. We applied these methods to study the genetic basis of body mass index levels and diabetes risk.

All scripts are written in the R or Bash programming languages. The analysis was performed on a linux system. Plots were created on a macOS Catalina system.

If you can not find the code you are looking for or have any questions, please contact:

Andrew Marderstein \
anm2868@med.cornell.edu


## Deviation Regression Model discovers vQTLs that are due to GxE interactions

Directory: ./scripts/vqtl_method_comparison.R

#### vQTL_method_compare_simulation.R
**Description:** Main script 1. Simulating population genetic data and comparing the false positive rates and power for different variance tests. Results are saved using various parameter settings.

#### vQTL_vs_muQTL_compare_simulation.R
**Description:** Main script 2. Simulating population genetic data and contrasting power for a muQTL test versus a vQTL test. Results are saved using various parameter settings.

#### Other scripts:
**Description:** (1) Supplementary figure of interaction effect size versus variance effect size. (2) Power heatmap (supp fig). (3) muQTL vs vQTL output & plots. (4) vQTL method comparison figure.

	1. beta_vs_variance_explained_boxplots.R
	2. heatmap_vg_vs_N.R
	3. muqtl_vqtl_plot.R
	5. vQTL_method_compare_simulation_save_draw_plots.R




## UKB GWAS of blood cell phenotypes
Directory: ./scripts/GWAS/



#### pre-processing:
Want to add another phenotype? Need to add the phenotype using gen_pheno, re-run gen_full_data and preprocess w/ updated PHENOTYPE_NAMES variable. good luck!

split_80_20.R \
identify_indiv_blood_disorders.R \
sample_qc.R \
gen_cov1.R \
gen_cov2.R \
gen_pheno.R \
gen_full_data.R \
preprocess.R

#### Script to perform DRM variance testing in Plink
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

#### run vGWAS -> GxG pipeline using R's BEDMatrix library
Use ./scripts/GWAS/subset/vQTLs_to_GxG_pipeline.sh


#### merge diff assoc files together:
mergeResults.R ########## note: might need to change directories for mean gwas analyses
mergeTransform.R

#### clump to significant hits, identify top snps, extract from genotype file
clump_identify_extract.sh



## GxG testing
Directory: ./scripts/GxG

#### pipeline
see GxG_pipeline.sh

#### analyze results
GxG_results_analyze.R




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


