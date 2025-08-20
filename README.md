# ukb_alcohol_gwas
Author: Jeanne Savage (j.e.savage@vu.nl)

This repository includes scripts for various alcohol-related genetic analyses I have carried out on the UK Biobank dataset under application #16406. GWAS summary statistics can be downloaded from https://cncr.nl/research/summary_statistics/.


## Genetic heterogeneity across dimensions of alcohol use behaviors
Results from [Savage et al., 2024, American Journal of Psychiatry](https://www.psychiatryonline.org/doi/10.1176/appi.ajp.20231055). A genomic structural equation model of 18 alcohol-related phenotypes revealed 4 underlying genetic factors which were carried forward to GWAS.

- ['SavageAJP2024_DerivePhenotypes.R'](SavageAJP2024_DerivePhenotypes.R) - Script for deriving phenotypes from UKB survey data and medical records
- ['SavageAJP2024_clinical_alcohol_maps.xlsx'](SavageAJP2024_clinical_alcohol_maps.xlsx) - Medical record phenotype definitions
- ['SavageAJP2024_fieldcodes.txt'](SavageAJP2024_fieldcodes.txt) - Survey phenotype definitions
- ['SavageAJP2024_gsem_script_alc.txt'](SavageAJP2024_gsem_script_alc.txt) - Script for carrying out gSEM factor analysis and GWAS 



## Refining the scope of genetic influences on alcohol misuse through environmental stratification and gene-environment interaction
Results from [Savage et al., 2024, Alcohol: Clinical and Experimental Research](https://onlinelibrary.wiley.com/doi/10.1111/acer.15425). A genome-wide-by-environment-interaction study of multiple alcohol use behaviors.

- ['SavageACER2024_stratified_gwas_trauma_ses.R'](SavageACER2024_stratified_gwas_trauma_ses.R) - Script for cleaning/recoding phenotypes and carrying out stratified GWAS and GWEIS.
- ['SavageACER2024_strat_codes.txt](SavageACER2024_strat_codes.txt) - Phenotype field codes
- ['SavageACER2024_anova2df_gxe_test.R'](SavageACER2024_anova2df_gxe_test.R) - Script to calculate a 2df ANOVA test of the main+interaction effect coefficients for comparison across different types of models.
- ['SavageACER2024_run_simGxE.sh'](SavageACER2024_run_simGxE.sh) - Script to simulate GxE effects and how they are captured under different types of models. 



## Investigating genetically stratified subgroups to better understand the origins of alcohol misuse
Results from [Thijssen et al., 2023, Molecular Psychiatry](https://www.nature.com/articles/s41380-023-02174-0). A latent class analysis and GWAS comparison of 4 subgroups of individuals based on patterns of alcohol use/problems and internalizing and externalizing psychopathology.

- ['ThijssenMP2023_DerivePhenotypes1.R'](ThijssenMP2023_DerivePhenotypes1.R) - Script to derive phenotypes for use in latent class analysis.
- ['ThijssenMP2023_fieldcodes.txt'](ThijssenMP2023_fieldcodes.txt) - Phenotype field codes
- ['ThijssenMP2023_DerivePhenotypes2.inp'](ThijssenMP2023_DerivePhenotypes2.inp) - Mplus script to carry out mixture modelling for deriving the latent classes.

