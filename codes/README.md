# Codes

Codes for replicating all the results presented in the manuscript

## Real Data analysis

We share all the codes used for generating the results of our real data applications. This pipeline can be applied to any GWAS summary data in principle. 

* step0 pre-process.R and step1 prepare-variant.R: prepare the GWAS summary statistics
* step2_0_genetic_cor:  using [LD score regression](https://github.com/bulik/ldsc) to calculate each element of the trait correlation matrix
* step3_prepare_Z_score:  divide GWAS Z score to several pieces to facilitate parallel computing
* step4a_MAT_construct_trait_cor.R: construct the trait correlation used for aMAT
* step4b_MAT.R: conduct multi-trait association tests (including aMAT and other popular ones)
* step5_combine_results.R: combine the results
* multitrait_test_support.R: supporting functions for our analyses. You can use the cleaned version in aMAT software. 

## Simulations

Because we did many simulations with many different settings (all these are just slightly different in terms of settings), we only provide some sample codes. 

* null_hypothesis.R: evaluate empirical type 1 error rates for different methods
* power_comparison.R: compare power for different methods
* time_comparison.R: compare running time for different methods
* submit_power_jobs.R: a wrapping code for submitting jobs more easily
