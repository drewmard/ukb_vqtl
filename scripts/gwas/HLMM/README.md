### 1: Run HLMM in parallel, in 5000-SNP sets

run_HLMM_subset.sh

### 2: Merge HLMM results. Will indicate if any HLMM runs failed and re-run for those that failed.

merge_HLMM_subset.sh

# You can also use the errors directory to check killed jobs: (if having cluster problems)

HLMM_Killed.txt

### 3: Perform DET and generate combined data frame with the vQTL and muQTL results.

spack load -r r@3.5.0
merge_HLMM_subset_2.R

### 4:

plots.R