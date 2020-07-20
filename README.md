# Leveraging phenotypic variability to identify genetic interactions in human phenotypes

This repository contains all code to conduct the analyses present in the manuscript, "Leveraging phenotypic variability to identify genetic interactions in human phenotypes". This research describes a statistical framework to find SNPs associated with the means and variances of a quantitative phenotype, and then using these SNPs to discover gene-environment interactions. We applied these methods to study the genetic basis of body mass index levels and diabetes risk.

All scripts are written in the R or Bash programming languages. The analysis was performed on a linux system. Plots were created on a macOS Catalina system.

If you can not find the code you are looking for or have any questions, please contact:

Andrew Marderstein \
anm2868@med.cornell.edu


## Deviation Regression Model discovers vQTLs that are due to GxE interactions

Directory: ./scripts/vqtl_method_comparison

#### (A) vQTL_method_compare_simulation.R
**Description:** Main script 1. Simulating population genetic data and comparing the false positive rates and power for different variance tests. Results are saved using various parameter settings.

#### (B) vQTL_vs_muQTL_compare_simulation.R
**Description:** Main script 2. Simulating population genetic data and contrasting power for a muQTL test versus a vQTL test. Results are saved using various parameter settings.

#### (C) Other scripts:
**Description:** (1) Supplementary figure of interaction effect size versus variance effect size. (2) Power heatmap (supp fig). (3) muQTL vs vQTL output & plots. (4) vQTL method comparison figure.

	1. beta_vs_variance_explained_boxplots.R
	2. heatmap_vg_vs_N.R
	3. muqtl_vqtl_plot.R
	4. vQTL_method_compare_simulation_save_draw_plots.R



## Genome-wide association studies in UK Biobank

Directory: ./scripts/gwas

### Directory: preprocess

#### (A) Various scripts
**Description:** Various scripts to partition UKB into discovery and validation cohorts (1), pull study covariates (2, 3, 4), and generate the full dataset for GWAS (5).

	1. split_80_20.R
	2. sample_qc.R
	3. gen_cov1.R
	4. gen_cov2.R
	5. gen_full_data.R
	
#### (B) preprocess.R
**Description:** Adjust phenotype for covariates prior to running a GWAS.


### Directory: run

#### (A) PIPELINE.sh
**Description:** This is the master pipeline script used to perform a GWAS for muQTLs, raw vQTLs, and RINT vQTLs in UKB. It was manually ran piece-by-piece. 

	1. run_GWAS.impute.sh
	2. run_vGWAS_subset.sh
	3. mergeResults_impute.R
	4. merge_vGWAS_subset.sh
	5. merge_vGWAS_subset_2.R
	6. ./errors/run_vGWAS_subset_specific.sh
	
Briefly, it runs a GWAS for muQTLs using (1) and for vQTLs using (2). Next, it merges the different files together, using (3) for muQTLs and (4, 5) for vQTLs. If an error occurs for vQTL script executions, then it will re-run using (6). Finally, the results are clumped and an all-by-all SNP-by-SNP epistasis analysis was performed from the clumped loci.


#### (B) run_vGWAS.R
**Description:** This is the script to perform the Deviation Regression Model on PLINK BED formatted genotype files. It is used within (A2).


### Directory: HLMM

**Description:** See *README.md* within the HLMM directory.

### output

#### (A) generate_sig_results.R
**Description:** This script merges the muQTL, raw vQTL, RINT vQTL, and dQTL results together, and extracts out the significant SNPs.

#### (B) mean_vs_var_plots.R
**Description:** Creates figures that visualize the results from the different GWAS analyses. For example, displaying muQTL effects versus raw vQTL effects.

#### (C) assign_genes.R
**Description:** Pipeline to map SNPs to the most likely gene.



## Discovery and validation of gene-environment interactions 

#### generate environmental factors
gen_envir_factors.R


## Large gene-environment interactions influence BMI


## GxE interactions have pleiotropic effects over BMI and diabetes risk



## Evidence for weak epistatic interactions associated with BMI
Directory: ./scripts/GxG

#### (A) GxG_results_analyze.R
**Description:** Analyzing the GxG results.

#### (B) disc_vs_valid.R
**Description:** Subsetting the tested GxG interactions and comparing the correlation of observed effects between discovery and validation sets.

#### (C) gxg_corr_plot_disc_vs_valid.R
**Description:** Plot of the results from (B).

#### (D) gxg_panel_plot.R
**Description:** Figure panel for the GxG results.

 
## vQTLs are linked to environmentally-influenced pathways and phenotypes
Directory: ./scripts/phewas

#### (A) variant_index_query.sh
**Description:** Maps list of SNPs from UK Biobank's SNP IDs to Open Target's variant IDs.

#### (B) pheWAS.R
**Description:** The pheWAS() function takes a variant ID and finds all associated phenotypes in Open Targets (P < 0.05).

#### (C) pheWAS_QTL_run.R
**Description:** Runs the pheWAS() function for a series of input SNPs and saves the results.

#### (D) pheWAS_enrich.R
**Description:** This script runs the PheWAS enrichment analysis, where the proportion of raw vQTLs associated with some phenotype are compared to the proportion of muQTLs associated with the same phenotype.

#### (E) phewas_plot.R
**Description:** Scatterplot figure displaying the PheWAS results.


## Polygenic heritability analysis implicates stomach cell types in regulating BMI variance

#### Various scripts
**Description:** Three scripts used to generate summary statistics for LDSC input (1), to run LDSC and receive results (2), and to visualize the output (3).

	1. ldsc_SS_generate.R
	2. ldsc_run.txt
	3. ldsc_plot.R

## Other scripts

#### (A) FPR_inflation_vqtl.R
**Description:** Generates a SNP associated with the variance of a phenotype with no true phenotypic effect, and assesses the false positive rate and visualizes the results.

#### (B) population_statistics_QTL_compare.R
**Description:** Comparison of several population attributes between muQTLs and raw vQTLs.

#### (C) variance_of_statistics.R
**Description:** Analysis of the standard error of a mean estimate and a variance estimate.

#### (D) ./match directory
**Description:** Match random genome-wide SNPs to the observed QTLs using population attributes, and create a genotype file in PLINK BED format.



