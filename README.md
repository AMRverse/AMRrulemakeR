<img src="AMRrules_logo.png" width="200" align="left">

# AMRrulesR: Apply or Define AMRrules for Interpreting AMR Genotypes


This package is in early development. Its purpose will be to apply or define AMRrules, which specify how the presence of an antimicrobial resistance (AMR) marker in the genome of a given bacterial species should be interpreted in terms of the expected phenotype.

An overview of the AMRrules concept, with example data structures and code, is available in the [AMRrules](https://github.com/interpretAMR/AMRrules) repository.

Development of the AMRrules scheme, and curation of organism-specific rule sets, is being undertaken by a working group of [ESGEM, the ESCMID Study Group on Epidemiological Markers](https://www.escmid.org/esgem/), known as the ESGEM-AMR Working Group. Details of this group and their activities can be found [here](https://github.com/interpretAMR/AMRrulescuration), including draft rules in the AMRrules format ([specification v0.5](https://docs.google.com/spreadsheets/d/1F-J-_8Kyo3W0Oh6eDYyd0N8ahqVwiddM2112-Fg1gKc/edit?usp=sharing)).

To support the work of defining AMRrules, we have developed the [AMRgen](https://github.com/interpretAMR/AMRgen) R package to facilitate comparison of AMR genotype and phenotype data and the defining of rules based on quantitative analysis. 

__The `AMRrulesR` package draws on the `AMR` and `AMRgen` package functions to:__
* summarise the available genotype-phenotype data for a given drug, and the available breakpoints with which to interpret phenotypes, using the `summarise_data` function
* automate runnning a series of genotype-phenotype analyses (including calcualting solo marker positive predictive value, generating upset plots, calculating distributions for MIC and disk diffusion assay values for individual markers and combinations, logistic regression of markers vs phenotypes etc), using the `amrrules_analysis` function
* use these quantitative results to draft AMRrules to interpret individual markers or combinations in terms of expected phenotypes (wildtype/nonwildtype) and clinical categories (S/I/R), according to the [AMRrules specification (v0.5)](https://docs.google.com/spreadsheets/d/1F-J-_8Kyo3W0Oh6eDYyd0N8ahqVwiddM2112-Fg1gKc/edit?usp=sharing), using the `makerules` function
* the function `amrrules_save` automates saving the quantitative results to files (TSV tables, PDF figures), and generating and saving AMRrules

During 2025, as the rules are deveoped, functions will be added to interpret AMRfinderplus genotyping results into inferred phenotypes, using the AMRrules defined by the ESGEM-AMR Working Group.


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
remotes::install_github("katholt/AMRrulesR")
```

It is best to restart R before running the installation. If you didn't do this and/or you encounter issues with the examples below after install, it may help to also restart after the install and start fresh with the examples below.


### Examples
TBD
