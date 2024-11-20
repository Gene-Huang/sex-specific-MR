This repository contains R code and summary statistics for implementing the simulations and real data analysis from our work: “A semi-empirical Bayes approach for calibrating weak instrumental bias in sex-specific Mendelian randomization studies”.

The repository is organized into the following folders:

1. "main_function" folder:
This folder contains the primary functions used in our analysis:
(i) SEB_function.R: Implements the proposed semi-empirical Bayes method for computing shrinkage estimation;
(ii) Freq_AW_function.R: Implements adaptive weight (AW) method for computing shrinkage estimation;
(iii) TwoSampleMR_methods_function.R: Computes causal effect of an exposure on an outcome using two-sample Mendelian Randomization (MR) methods.

2. "simulation" folder:
This folder includes the codes for conducting simulation studies to evaluate the proposed methods.

3. "data_analysis" folder:
This folder contains codes for estimating the sex-specific causal effects of sleep-related phenotypes on cardiovascular disease (CVD)-related outcomes.

4. "harmonized_exposure_outcome_sumstat" folder: 
This folder provides the harmonized summary statistics from exposure and outcome genome-wide association studies (GWAS) used in our real data analysis.
