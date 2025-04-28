<img src="AMRrules_logo.png" width="200" align="left">

# AMRrulesR: Define and Test AMRrules for Interpreting AMR Genotypes


This package is in early development. Its purpose will be to define and test AMRrules, which specify how the presence of an antimicrobial resistance (AMR) marker in the genome of a given bacterial species should be interpreted in terms of the expected phenotype.

An overview of the AMRrules concept, with example data structures and code, is available in the [AMRrules](https://github.com/interpretAMR/AMRrules) repository.

Development of the AMRrules scheme, and curation of organism-specific rule sets, is being undertaken by a working group of [ESGEM, the ESCMID Study Group on Epidemiological Markers](https://www.escmid.org/esgem/), known as the ESGEM-AMR Working Group. Details of this group and their activities can be found [here](https://github.com/interpretAMR/AMRrulescuration), including draft rules in the AMRrules format ([specification v0.5](https://docs.google.com/spreadsheets/d/1F-J-_8Kyo3W0Oh6eDYyd0N8ahqVwiddM2112-Fg1gKc/edit?usp=sharing)).

To support the work of defining AMRrules, we have developed the [AMRgen](https://github.com/interpretAMR/AMRgen) R package to facilitate comparison of AMR genotype and phenotype data and the defining and checking of rules based on quantitative analysis. 

__The `AMRrulesR` package draws on the `AMR` and `AMRgen` package functions to:__
* summarise the available genotype-phenotype data for a given drug, and the available breakpoints with which to interpret phenotypes, using the `summarise_data` function
* automate runnning a series of genotype-phenotype analyses (including calcualting solo marker positive predictive value, generating upset plots, calculating distributions for MIC and disk diffusion assay values for individual markers and combinations, logistic regression of markers vs phenotypes etc), using the `amrrules_analysis` function
* use these quantitative results to draft AMRrules to interpret individual markers or combinations in terms of expected phenotypes (wildtype/nonwildtype) and clinical categories (S/I/R), according to the [AMRrules specification (v0.5)](https://docs.google.com/spreadsheets/d/1F-J-_8Kyo3W0Oh6eDYyd0N8ahqVwiddM2112-Fg1gKc/edit?usp=sharing), using the `makerules` function
* automate generating and saving AMRrules, and saving the quantitative results to files (TSV tables, PDF figures), using the `amrrules_save` function
* perform basic testing/validation of the auto-defined rules by applying them to interpret AMRfinderplus results into wildtype/nonwildtype S/I/R, using the `test_rules_amrfp` function (Note that this is separate from the main [AMRrules interpretation engine](https://github.com/interpretAMR/AMRrules), which is being developed in Python.)


## Getting Started

To install and explore the AMRrulesR package, follow the instructions below:

### Installation
Note that this package requires the latest version of the `AMR` and `AMRgen` packages.

Install the latest version of these with:
```r
install.packages("remotes") # if you haven't already
remotes::install_github("msberends/AMR")
remotes::install_github("interpretAMR/AMRgen")
```

Then install this package
```r
# Install from GitHub
remotes::install_github("interpretAMR/AMRrulesR")
```

It is best to restart R before running the installation. If you didn't do this and/or you encounter issues with the examples below after install, it may help to also restart after the install and start fresh with the examples below.


### Examples
```r

library(AMRrules)
library(tidyverse)

# example data from AMRgen package: E. coli MIC data from NCBI, matching AMRfinderplus data
ecoli_ast
ecoli_geno <- import_amrfp(ecoli_geno_raw, "Name")

# run quantitative analyses
cip_analysis <- amrrules_analysis(ecoli_geno, ecoli_ast, antibiotic="Cipro", drug_class_list=c("Quinolones"), species="E. coli")

cip_analysis$ppv_plot
cip_analysis$upset_mic_plot
cip_analysis$upset_disk_plot

# save tables and plots and generate rules
cip_rules <- amrrules_save(cip_analysis, bp_site="Non-meningitis", dir_path="amrrules", file_prefix="Cipro")

# test rules
calls <- test_rules_amrfp(geno_table = ecoli_geno, rules = cip_rules$rules, species = "Escherichia coli")

# compare these calls to the AST data phenotypes in a separate dataframe, `pheno_table` with SIR phenotypes in `pheno`
calls_vs_pheno <- calls %>% left_join(ecoli_ast, join_by(Name==id))
calls_vs_pheno %>% group_by(pheno, category) %>% count() %>% filter(pheno %in% c("S", "I", "R")) 
```
