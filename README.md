<img src="AMRrules_logo.png" width="200" align="left">

# AMRrulemakeR: Define and Test AMRrules for Interpreting AMR Genotypes


This package is in early development. Its purpose will be to define and test AMRrules, which specify how the presence of an antimicrobial resistance (AMR) marker in the genome of a given bacterial species should be interpreted in terms of the expected phenotype.

An overview of the AMRrules concept, with example data structures and code, is available in the [AMRrules](https://github.com/AMRverse/AMRrules) repository.

Development of the AMRrules scheme, and curation of organism-specific rule sets, is being undertaken by a working group of [ESGEM, the ESCMID Study Group on Epidemiological Markers](https://www.escmid.org/esgem/), known as the ESGEM-AMR Working Group. Details of this group and their activities can be found [here](https://github.com/AMRverse/AMRrulescuration), including draft rules in the AMRrules format ([specification v0.6](https://docs.google.com/spreadsheets/d/1t6Lr_p-WAOY0yAXWKzoKk4yb56D2JdSqwImg4RZBvFA/edit?usp=sharing)).

To support the work of defining AMRrules, we have developed the [AMRgen](https://github.com/AMRverse/AMRgen) R package to facilitate comparison of AMR genotype and phenotype data and the defining and checking of rules based on quantitative analysis. 

__The `AMRrulemakeR` package draws on the `AMR` and `AMRgen` package functions to:__
* summarise the available genotype-phenotype data for a given drug, and the available breakpoints with which to interpret phenotypes, using the `summarise_data` function
* automate runnning a series of genotype-phenotype analyses (including calcualting solo marker positive predictive value, generating upset plots, calculating distributions for MIC and disk diffusion assay values for individual markers and combinations, logistic regression of markers vs phenotypes etc), using the `amrrules_analysis` function
* use these quantitative results to draft AMRrules to interpret individual markers or combinations in terms of expected phenotypes (wildtype/nonwildtype) and clinical categories (S/I/R), according to the [AMRrules specification (v0.6)](https://docs.google.com/spreadsheets/d/1t6Lr_p-WAOY0yAXWKzoKk4yb56D2JdSqwImg4RZBvFA/edit?usp=sharing), using the `makerules` function
* automate generating and saving AMRrules, and saving the quantitative results to files (TSV tables, PDF figures), using the `amrrules_save` function
* perform basic testing/validation of the auto-defined rules by applying them to interpret AMRfinderplus results into wildtype/nonwildtype and S/I/R, using the `test_rules_amrfp` function (Note that this is separate from the main [AMRrules interpretation engine](https://github.com/AMRverse/AMRrules), which is being developed in Python.)


## Getting Started

To install and explore the AMRrulemakeR package, follow the instructions below:

### Installation
Note that this package requires the latest version of the `AMR` and `AMRgen` packages.

Install the latest version of these with:
```r
install.packages("remotes") # if you haven't already
remotes::install_github("msberends/AMR")
remotes::install_github("AMRverse/AMRgen")
```

Then install this package
```r
# Install from GitHub
remotes::install_github("AMRverse/AMRrulemakeR")
```

It is best to restart R before running the installation. If you didn't do this and/or you encounter issues with the examples below after install, it may help to also restart after the install and start fresh with the examples below.


### Examples
```r
library(AMRgen)
library(AMRrulemakeR)
library(tidyverse)

# example data included in AMRrulemakeR package
ecoli_ast_ebi
ecoli_afp_atb

# run quantitative analyses for ciprofloxacin phenotypes vs quinolone marker genotypes, using the S/I/R calls made against EUCAST breakpoints
cip_analysis <- amrrules_analysis(geno_table=ecoli_afp_atb, pheno_table=ecoli_ast_ebi, 
                                    antibiotic="Ciprofloxacin", drug_class_list=c("Quinolones"),
                                    sir_col="pheno_eucast", ecoff_col="ecoff",
                                    species="Escherichia coli",
                                    info=ecoli_ast_ebi %>% select(id, source, method))

# check key output plots
cip_analysis$ppv_plot
cip_analysis$ppv_plot_all # this includes data with S/I/R interpretations from EBI but no raw assay values (treated as the 'extended' dataset in the analysis)
cip_analysis$logistic_plot # note this is only used if the marker is not found solo, to support a call of WT S based on lack of association with resistance in the regression
cip_analysis$upset_mic_plot

# save analysis tables and plots, and generate rules using the Non-meningitis breakpoint, save output to 'amrrules/' with filenames starting with 'Ciprofloxacin'
# then use the rules to predicted phenotypes from genotypes and compare to the observed phenotypes (to help try to spot issues with input data and proposed rules)
cip_rules <- amrrules_save(cip_analysis, bp_site="Non-meningitis",
                           dir_path="amrrules", file_prefix="Ciprofloxacin")

# alternatively, call makerules directly on the analysis object without saving outputs or running predictions
cip_rules <- makerules(cip_analysis, bp_site="Non-meningitis")

# view the proposed rules, in AMRrules specification format, with quantitative fields added
view(cip_rules$rules)

# manually apply rules to interpret quinolone marker genotypes
cip_test <- test_rules_amrfp(ecoli_afp_atb %>% filter(drug_class %in% c("Quinolones")),
                                 rules=cip_rules$rules, species="Escherichia coli")

# compare these to the input phenotypes
cip_test %>% left_join(ecoli_ast_ebi, join_by("Name"=="id")) %>% count(category,pheno_eucast)

# make some plots to view predictions vs raw assay values from different methods, and
# explore positive predictive value of the overall ruleset including stratified by method 
compare_pred <- compare_interpretations(pred=cip_test, obs=ecoli_ast_ebi,
                                        antibiotic="Ciprofloxacin",
                                        sir_col="pheno_eucast", ecoff_col="ecoff",
                                        var="method")

# check positive predictive value of rules
compare_pred$pred_ppv$plot_sir

# check positive predictive value of rules, stratified by assay method
compare_pred$pred_ppv_bymethod$plot_sir

# review MIC distribution, stratified by assay method and coloured by prediction
compare_pred$dist_mic_bypred_bymethod$pred

```

# Work in progress - suggested protocol for developing AMRrules using this package

## Collate and format phenotype data
For use with the AMRgen & AMRrulemakeR packages, you need to get the antimicrobial susceptibility AST (phenotype) data into long format (one row per sample/drug result), with some key fields.
Data in NCBI or EBI antibiogram format can be automatically imported to the right dataframe using the `import_ast()` in the AMRgen package.
Or you can format your data manually. 

Example data frame containg data from EBI on E. coli with AST results for five drugs, imported using import_ast():

```
> ecoli_ast_ebi %>% filter(!is.na(mic))
# A tibble: 43,288 × 47
   id        drug_agent   mic  disk pheno_eucast pheno_clsi ecoff guideline method source pheno_provided spp_pheno   
   <chr>     <ab>       <mic> <dsk> <sir>        <sir>      <sir> <chr>     <chr>  <chr>  <sir>          <mo>        
 1 SAMN1302… AMP           >8    NA   R            I          R   EUCAST    BD Ph… 32205…   R            B_ESCHR_COLI
 2 SAMN1302… AMP           >8    NA   R            I          R   EUCAST    BD Ph… 32205…   R            B_ESCHR_COLI
 3 SAMN1302… AMP           >8    NA   R            I          R   EUCAST    BD Ph… 32205…   R            B_ESCHR_COLI
 4 SAMN1302… AMP          <=2    NA   S            S          S   EUCAST    BD Ph… 32205…   S            B_ESCHR_COLI
 5 SAMN1302… AMP           >8    NA   R            I          R   EUCAST    BD Ph… 32205…   R            B_ESCHR_COLI
 6 SAMN1302… AMP           >8    NA   R            I          R   EUCAST    BD Ph… 32205…   R            B_ESCHR_COLI
 7 SAMN1302… AMP           >8    NA   R            I          R   EUCAST    BD Ph… 32205…   R            B_ESCHR_COLI
 8 SAMN1302… AMP          <=2    NA   S            S          S   EUCAST    BD Ph… 32205…   S            B_ESCHR_COLI
 9 SAMN1302… AMP          <=2    NA   S            S          S   EUCAST    BD Ph… 32205…   S            B_ESCHR_COLI
10 SAMN1302… AMP           >8    NA   R            I          R   EUCAST    BD Ph… 32205…   R            B_ESCHR_COLI
```

The key fields needed are:
* `id` - Sample name, must match that in the corresponding genotype file.
    * If using `import_ast()`, this is taken from the biosample field.
* `drug_agent` - Name of the antibiotic, formatted as class 'ab' (using `as.ab()`)
* `mic` - MIC assay measurement (where available), formated as class 'mic' (using `as.mic()`). Note that if your input file has MIC value in one column, and sign (e.g. >, <, <=, etc) in another column (e.g. the NCBI antibiogram format), you will need to paste those two columns together first before applying `as.mic()`.
* `disk` - Disk assay measurement (where available), formated as class 'disk' (using `as.disk()`).
* `pheno_eucast`, `pheno_clsi` (can be any name) - A column of class 'sir', indicating the interpretation of the data in `mic` and `disk` as S/I/R, using clinical breakpoints.
    * If using `import_ast()`, you can set `interpret_eucast=TRUE` and/or `interpret_clsi=TRUE` to generate these fields automatically by interpreting input MIC and disk assay measures against the EUCAST and CLSI breakpoints via the `as.sir()` function in the `AMR` package.
    * This field will be used as the primary S/I/R data to calculate positive predictive value and defining rules, so should only include calls based on assay measurements such that you can ensure they have been called consistently.
    * Consider whether disk content has changed since the assays were done, before re-calling disk data with latest breakpoints.
* `pheno_provided` (can be any name) - A column of class 'sir', indicating S/I/R calls for **all** samples, even those for which the raw assay measures are not available in `mic` or `disk`.
    * This data will be used as a secondary, 'extended' data source when there is insufficient primary data (with raw assay measures) to define a rule. Such data should be included only when it comes from a reliable source, where you know that the assay and interpretations are trustworthy, and the breakpoints for the relevant bug/drug combinations have not changed since the calls were made. For disk data, it is also important to consider whether disk content has changed since the assays were done.
* `ecoff` (can be any name) - A column of class 'sir', indicating the interpretation against ECOFF as S/R
    * If using `import_ast()`, you can set `interpret_ecoff=TRUE` to generate this field automatically by interpreting input MIC and disk assay measures against the EUCAST ECOFF via the `as.sir()` function in the `AMR` package.
* A column indicating the source of each dataset (this can be any kind of indicator you like, such as bioproject accession, PubMed ID, project name) (e.g. `source`)
    * If using `import_ast()`, this is taken from the bioproject field (importing from NCBI format) or pubmed field (importing from EBI format)
* A column indicating the MIC assay method used to generate the measurement (e.g. microbroth dilution, Vitek, Sensititre, BD Phoenix, etc)
    * If using `import_ast()`, this is taken from the 'Laboratory typing platform' field (importing from NCBI format) or 'phenotype-platform' field (importing from EBI format)

### Examples

```
# download NCBI data from https://www.ncbi.nlm.nih.gov/pathogens/ast#taxgroup_name:%22E.coli%20and%20Shigella%22
# import data from file
ncbi_ast <-import_ncbi_ast("NCBI_AST_EcoliShigella.tsv.gz")

# download EBI data from https://www.ebi.ac.uk/amr/data/?view=experiments
# import data from file
ebi_ast <- import_ebi_ast("EBI_CABBAGE_EcoliShigella.csv.gz")

# import study data and format it
study1_ast <- read_tsv("AST.txt") %>%
   rename(id="Strain name") %>%
   mutate(mic = paste0(`Measurement sign`, `MIC (mg/L)`)) %>% # combine measurement sign with MIC value into single column
   mutate(mic = gsub("==", "", mic)) %>%
   mutate(mic = as.mic(mic)) %>% # format MIC values to class 'mic'
   mutate(drug_agent = as.ab(Antibiotic)) %>% # format drug name to class 'ab'
   mutate(spp_pheno = as.mo(`Species name`)) %>% # format organism name to class 'mo' (optional unless you have data from multiple organisms in the same dataframe)
   mutate(method="Sensititre") %>%
   mutate(source="Study1")

# interpret MIC data using EUCAST breakpoints and ECOFFs
study1_ast <- interpret_ast(study1_ast, interpret_ecoff = TRUE, interpret_eucast = TRUE, interpret_clsi = FALSE, species = species, ab = ab)
```

## Collate and format AMRfinderplus genotype data
For use with the AMRrulemakeR package, you need to process the AMRfinderplus data using `import_afp()`, which generates consistent marker labels that will be the units of analysis, and which can ultimately be represented in the AMRrules variant specification format. It also parses refgene Class/Subclass into drug classes recognised by the AMR R packagea and CARD ARO. 

Example data frame included in the AMRgen package:

```
> ecoli_afp_atb
# A tibble: 73,221 × 36
   Name         gene  mutation node  `variation type` marker marker.label drug_agent drug_class status
   <chr>        <chr> <chr>    <chr> <chr>            <chr>  <chr>        <ab>       <chr>      <chr> 
 1 SAMN03177641 blaC… NA       blaC… Gene presence d… blaCM… blaCMY-2     NA         Cephalosp… PASS  
 2 SAMN03177641 sul2  NA       sul2  Gene presence d… sul2   sul2         SSS        Sulfonami… PASS  
 3 SAMN03177641 blaT… NA       blaT… Gene presence d… blaTE… blaTEM-1     NA         Beta-lact… PASS  
 4 SAMN13024031 dfrA5 NA       dfrA5 Gene presence d… dfrA5  dfrA5        NA         Trimethop… PASS  
 5 SAMN13024031 sul2  NA       sul2  Gene presence d… sul2   sul2         SSS        Sulfonami… PASS  
 6 SAMN13024031 blaT… NA       blaT… Gene presence d… blaTE… blaTEM-1     NA         Beta-lact… PASS  
 7 SAMN26308528 sul2  NA       sul2  Gene presence d… sul2   sul2         SSS        Sulfonami… PASS  
 8 SAMD00499563 gyrA  Asp87Asn gyrA  Protein variant… gyrA_… gyrA:Asp87A… NA         Quinolones PASS  
 9 SAMD00499563 gyrA  Ser83Leu gyrA  Protein variant… gyrA_… gyrA:Ser83L… NA         Quinolones PASS  
10 SAMD00499563 parC  Glu84Val parC  Protein variant… parC_… parC:Glu84V… NA         Quinolones PASS  
```

The key fields needed are:
* `Name` - Sample name, must match that in the corresponding phenotype file.
    * If using Allthebacteria genotype files this will be the biosample, facilitating matches to public phenotype data.
* `gene`, `mutation`, `variation type` - Variant specification fields in AMRrules format, inferred from the `Gene symbol` and `Method` fields
* `node` - NCBI reference gene hierarchy node ID corresponding to this hit, taken from `Hierarchy node` if available (otherwise copied from gene)
* `marker` - Original content of `Gene symbol` field, ie matching the entry in NCBI refgene
* `drug_agent` - Specific drug to which this marker applies (if an individual drug is named in the `Subclass` fields)
* `drug_class` - Drug class to which this marker applies, mapped to classes in AMR package (taken from `Class` and `Subclass` fields)

Notes: 
1. Currently the AMR package lacks some key classes of relevance to us, e.g. Fosfomycin. These will be added in future.
2. Currently the AMRrulemakeR package functions extract markers by drug_class, but we should consider updating to exclude for markers that are associated with a different drug of the same class. Or perhaps test with these included vs excluded, as the specificity of drug assignments in NCBI may be imperfect.

### Examples

```
# process AMRfp genotype data from Allthebacteria (this input file is not distributed with the AMRrulemakeR package)

# import AMRfp results file for E. coli
afp <-read_tsv("ATB_AFP_Ecoli_AMR.tsv.gz")

# read species calls so we can ensure only those confirmed as the target species will be included
species_calls <- read_tsv("ATB_species_calls.tsv.gz")
ecoli <- species_calls %>% filter(Species=="Escherichia coli") %>% pull(Sample) # list of confirmed E. coli

# read AMRfp run status for E. coli genomes
afp_status <- read_tsv("genotypes/ATB_AMRFP_status.tsv.gz") %>% filter(sample %in% ecoli)

# AMRfp results including null row for each genome that ran but returned no hits
afp_all <- afp_status %>% rename(Name=sample) %>% left_join(afp) %>% filter(status=="PASS")

# include only those genomes that we have ast data for in our `ast` object
afp_matching <- afp_all %>% filter(Name %in% ast$id)

# import to AMRgen format
afp_matching <- import_afp(afp_matching)
```

## Review the available data and breakpoints

It's a good idea to review the sources of data available, and the breakpoints available to interpret them.

### Examples of tricky breakpoints

There are multiple MIC breakpoints for ciprofloxacin, Meningitis and Non-meningitis:
```
> getBreakpoints(species="Escherichia coli", guide="EUCAST 2025", antibiotic="Ciprofloxacin")
# A tibble: 3 × 14
  guideline   type  host  method site    mo               rank_index ab   ref_tbl  disk_dose breakpoint_S breakpoint_R
  <chr>       <chr> <chr> <chr>  <chr>   <mo>                  <dbl> <ab> <chr>    <chr>            <dbl>        <dbl>
1 EUCAST 2025 human human DISK   Non-me… B_[ORD]_ENTRBCTR          5 CIP  Enterob… 5 mcg           25           22    
2 EUCAST 2025 human human MIC    Non-me… B_[ORD]_ENTRBCTR          5 CIP  Enterob… NA               0.25         0.5  
3 EUCAST 2025 human human MIC    Mening… B_[ORD]_ENTRBCTR          5 CIP  Enterob… NA               0.125        0.125
```

EUCAST guidance notes the Meningitis breakpoint is set lower to try to identify strains with any resistance mechanism, so this is redundant for genotype interpretation.
When we run rules, we can set the breapoint side with `bp_site="Non-meningitis` so we retrieve the right breakpoint. Or we can set the breakpoints manually as `mic_S=0.25` and `mic_R=0.5`.

There are no breakpoints for azithromycin, but there is an MIC ECOFF we can use instead to define rules, 
by setting `mic_S=16`, `mic_R=16` and `use_disk=F` as we have no way to intpret the disk values.
```
> getBreakpoints(species="Escherichia coli", guide="EUCAST 2025", antibiotic="Azithromycin")
# A tibble: 0 × 14
# ℹ 14 variables: guideline <chr>, type <chr>, host <chr>, method <chr>, site <chr>, mo <mo>, rank_index <dbl>, ab <ab>, ref_tbl <chr>, disk_dose <chr>,
#   breakpoint_S <dbl>, breakpoint_R <dbl>, uti <lgl>, is_SDD <lgl>

> getBreakpoints(species="Escherichia coli", guide="CLSI", antibiotic="Azithromycin")
# A tibble: 0 × 14
# ℹ 14 variables: guideline <chr>, type <chr>, host <chr>, method <chr>, site <chr>, mo <mo>, rank_index <dbl>,
#   ab <ab>, ref_tbl <chr>, disk_dose <chr>, breakpoint_S <dbl>, breakpoint_R <dbl>, uti <lgl>, is_SDD <lgl>

> getBreakpoints(species="Escherichia coli", guide="EUCAST 2025", antibiotic="Azithromycin", type_filter = "ECOFF")
# A tibble: 1 × 14
  guideline   type  host  method site  mo           rank_index ab   ref_tbl disk_dose breakpoint_S breakpoint_R uti   is_SDD
  <chr>       <chr> <chr> <chr>  <chr> <mo>              <dbl> <ab> <chr>   <chr>            <dbl>        <dbl> <lgl> <lgl> 
1 EUCAST 2025 ECOFF ECOFF MIC    NA    B_ESCHR_COLI          2 AZM  ECOFF   NA                  16           16 FALSE FALSE
```

### Exploring data from different platforms
We can use the `assay_by_var` function to plot the distribution of MIC or disk measures vs another variable, like assay platform or data source.

For example we can check ciprofloxacin MIC distributions by assay platform, highlighting values that are expressed as ranges.
```
cip_mic_bymethod <- assay_by_var(ecoli_ast_ebi, antibiotic="Ciprofloxacin", measure="mic", var="method",
                           species="Escherichia coli", bp_site="Non-meningitis")

cip_mic_bymethod$plot

# summarise BD Phoenix data
ecoli_ast_ebi %>% filter(drug_agent==as.ab("Ciprofloxacin") & method=="BD Phoenix") %>% count(mic,pheno_eucast)
# A tibble: 6 × 3
    mic pheno_eucast     n
  <mic> <sir>        <int>
1 <=0.5   NI          1588
2   1.0   R            125
3   2.0   R              9
4  >2.0   R            341
5  16.0   R              1
6 >16.0   R              1
```

This highlights that a lot of the automated platforms report MIC data as capped ranges (values highlighted in red). Many of these are still interpretable against the breakpoints, but most of the BD Phoenix data is recorded as '<=0.5', which can't be interpreted against the breakpoints S <=0.25, R >0.5, ECOFF 0.064. These values won't be used in the AMRrules analysis as they are not interpretable (NI).

## Run analyses need to define rules
The function `amrrules_analysis()` takes our phenotype table and extracts the data for a specified drug; then takes our genotype table and extracts the data for the relevant markers; and compares these geno/pheno data using several different `AMRgen` functions.

The function `amrrules_save()` takes these analysis results, and saves key tables and figures. It then uses these results to define rules in AMRrules format via the `makerules()` function, and writes these to an output TSV file. Finally, the rules are applied to the input genotypes to predict S/I/R and wildtype/nonwildtype for each sample using the `test_rules_amrfp()` function, and these results are compared to the input phenotypes and summarised in terms of positive predictive value (including stratified by assay method) to help assess the validity of the rules.

Example command
```
# example data included in AMRrulemakeR package: EBI AST data for 19,797 E. coli tested against at least one of 5 drugs (ampicillin, ceftriaxone, ciprofloxacin, azithromycin, trimethoprim-sulfamethozole), with matching AMRfinderplus data from the Allthebacteria project:
ecoli_ast_ebi
ecoli_afp_atb

# extract the information fields we want to consider in the analyses (source, method)
info_obj <- ecoli_ast_ebi %>% select(id, source, method)

# run the required analyses to compare phenotypes for a specific drug (e.g. ciprofloxacin, excluding BD Phoenix measures)
# with genetic markers associated with the corresponding drug class (e.g. quinolones)
cip_analysis <- amrrules_analysis(geno_table=ecoli_afp_atb,
                                    pheno_table=ecoli_ast_ebi, 
                                    antibiotic="Ciprofloxacin",
                                    drug_class_list=c("Quinolones"),
                                    sir_col="pheno_eucast", ecoff_col="ecoff",
                                    species="Escherichia coli",
                                    minPPV=1, mafLogReg=5, mafUpset=1,
                                    info=info_obj)

# check key output plots
cip_analysis$ppv_plot
cip_analysis$ppv_plot_all # this includes data with S/I/R interpretations from EBI but no raw assay values (treated as the 'extended' dataset in the analysis)
cip_analysis$logistic_plot # note this is only used if the marker is not found solo, to support a call of WT S based on lack of association with resistance in the regression
cip_analysis$upset_mic_plot

# note this has very litte information content as there is very limited public disk data for ciprofloxacin (n=240)
cip_analysis$upset_disk_plot 
ecoli_ast_ebi %>% filter(drug_agent==as.ab("Ciprofloxacin")) %>% filter(!is.na(disk))

# use the results of these analyses to define rules, then apply the rules back to the data to predict phenotypes
# write out files and figures for the analysis, assay distributions, and predicted vs observed phenotypes
cip_rules <- amrrules_save(cip_analysis,
                           dir="amrrules",
                           bp_site="Non-meningitis",
                           ruleID_start=1,
                           use_mic=TRUE, use_disk=TRUE,
                           file_prefix="Ciprofloxacin")

# view the proposed rules, in AMRrules specification format, with quantitative fields added
view(cip_rules$rules)

# check the positive predictive value of the rules for interpreting the input genotypes
cip_rules$predict_vs_obs_stats$plot_sir             # alldata
cip_rules$predict_vs_obs_stats_byMethod$plot_sir    # stratified by assay method
```

### Examples with markers for multiple drug classes

**Combination antibiotics:** Trimethoprim-sulfamethoazole is a commonly used combination drug. To define rules for its resistance, we need to consider markers associated with resistance to both classes included in the combination, via `drug_class_list=c("Trimethoprims", "Sulfonamides")`:
```
trimsulfa_analysis <- amrrules_analysis(geno_table=ecoli_afp_atb,
                                    pheno_table=ecoli_ast_ebi, 
                                    antibiotic="Trimethoprim-Sulfamethoxazole",
                                    drug_class_list=c("Trimethoprims", "Sulfonamides"),
                                    sir_col="pheno_eucast", ecoff_col="ecoff",
                                    species="Escherichia coli",
                                    minPPV=1, mafLogReg=5, mafUpset=1,
                                    info=info_obj)

# solo PPV shows individual markers don't confer clinical resistance
trimsulfa_analysis$ppv_plot

# upset plots show some marker combinations lead to resistance
trimsulfa_analysis$upset_mic_plot
trimsulfa_analysis$upset_disk_plot

# define rules and assess them
trimsulfa_rules <- amrrules_save(trimsulfa_analysis,
                           dir="amrrules",
                           file_prefix="TrimSulfa")

view(trimsulfa_rules$rules)
trimsulfa_rules$predict_vs_obs_stats$plot_sir
trimsulfa_rules$predict_vs_obs_stats_byMethod$plot_sir
```

**Beta-lactam antibiotics:** Beta-lactamases classified as cephalosporinases (class 'Cephalosporins' in NCBI refgene) have activity against narrow-spectrum beta-lactam drugs (such as ampicillin) as well as cephalosporins. Beta-lactamases classified as carbapenemases (class 'Carbapenems' in NCBI refgene) also have activity against narrow-spectrum beta-lactam drugs and cephalosporins as well as carbapenems. So to analyse geno-pheno relationships for beta-lactam drugs like ampicillin, we need to consider all three classes of markers in AMRfinderplus output that have activity against beta-lactams: `drug_class_list=c("Beta-lactams/penicillins", "Cephalosporins", "Carbapenems")`. 
(Note the ecoli_afp_atb example genotypes have already had blaEC calls removed, as this a core gene that doesn't contribute to resistance.)
```
amp_analysis <- amrrules_analysis(geno_table=ecoli_afp_atb,
                                    pheno_table=ecoli_ast_ebi, 
                                    antibiotic="Ampicillin",
                                    drug_class_list=c("Beta-lactams/penicillins", "Cephalosporins", "Carbapenems"),
                                    sir_col="pheno_eucast", ecoff_col="ecoff",
                                    species="Escherichia coli",
                                    minPPV=1, mafLogReg=5, mafUpset=1,
                                    info=info_obj)

# the solo PPV plot confirms that, as expected, ESBL genes (like blaCTX-M-15) or carbapenemases (blaKPC-3) found solo without other narrow-spectrum beta-lactamases are associated with resistance to ampicillin, and therefore need rules specified for ampicillin and other penicillins
amp_analysis$ppv_plot

# define rules and assess them
amp_rules <- amrrules_save(amp_analysis,
                           dir="amrrules",
                           file_prefix="Ampicillin")

view(amp_rules$rules)
amp_rules$predict_vs_obs_stats$plot_sir
amp_rules$predict_vs_obs_stats_byMethod$plot_sir

```

For cephalosporin drugs we also need to consider all three classes, since many carbapenemases are expected to act on cephalosporins, and in principle any beta-lactamases can have some activity against cephalosporins that in combination with other markers may be clinically relevant. It is also important to define rules for beta-lactamases that do NOT confer resistance when found solo, so we can explicitly encode 'WT S' rules for these. The same logic applies for defining rules for carbapenems.

```
cef_analysis <- amrrules_analysis(geno_table=ecoli_afp_atb,
                                    pheno_table=ecoli_ast_ebi, 
                                    antibiotic="Ceftriaxone",
                                    drug_class_list=c("Beta-lactams/penicillins", "Cephalosporins", "Carbapenems"),
                                    sir_col="pheno_eucast", ecoff_col="ecoff",
                                    species="Escherichia coli",
                                    minPPV=1, mafLogReg=5, mafUpset=1,
                                    info=info_obj)
# define rules and assess them
cef_rules <- amrrules_save(cef_analysis,
                           dir="amrrules",
                           file_prefix="Ceftriaxone")

view(cef_rules$rules)
cef_rules$predict_vs_obs_stats$plot_sir
cef_rules$predict_vs_obs_stats_byMethod$plot_sir
```

### Supplying manual breakpoints
If there are no clinical breakpoints, it may make sense to use the ECOFF as the breakpoint to define S/I/R. 

For example, as noted above there are no breakpoints for azithromycin, but there is an MIC ECOFF (16 mg/L) that we can use instead to define rules.

We could either do this by manually defining a phenotype column in our `pheno_table` input dataframe before running the analysis, or we can just supply the parameters `mic_S=16` and `mic_R=16` to the amrrules functions and this will be done automatically for us without modifying the input table.

**Note** if you are supplying your own breakpoints you will need to modify the rules file to record the source of your breakpoints in the `breakpoint standard` field, as the default is to assume the breakpoints you've used are standard ones.
```
> getBreakpoints(species="Escherichia coli", guide="EUCAST 2025", antibiotic="Azithromycin")
# A tibble: 0 × 14
# ℹ 14 variables: guideline <chr>, type <chr>, host <chr>, method <chr>, site <chr>, mo <mo>, rank_index <dbl>, ab <ab>, ref_tbl <chr>, disk_dose <chr>,
#   breakpoint_S <dbl>, breakpoint_R <dbl>, uti <lgl>, is_SDD <lgl>

> getBreakpoints(species="Escherichia coli", guide="EUCAST 2025", antibiotic="Azithromycin", type_filter = "ECOFF")
# A tibble: 1 × 14
  guideline   type  host  method site  mo           rank_index ab   ref_tbl disk_dose breakpoint_S breakpoint_R uti   is_SDD
  <chr>       <chr> <chr> <chr>  <chr> <mo>              <dbl> <ab> <chr>   <chr>            <dbl>        <dbl> <lgl> <lgl> 
1 EUCAST 2025 ECOFF ECOFF MIC    NA    B_ESCHR_COLI          2 AZM  ECOFF   NA                  16           16 FALSE FALSE

azi_analysis <- amrrules_analysis(geno_table=ecoli_afp_atb,
                                    pheno_table=ecoli_ast_ebi, 
                                    antibiotic="Azithromycin",
                                    drug_class_list=c("Macrolides", "Macrolides/lincosamides"),
                                    sir_col="pheno_eucast", ecoff_col="ecoff",
                                    species="Escherichia coli",
                                    minPPV=1, mafLogReg=5, mafUpset=1,
                                    info=info_obj,
                                    mic_S=16, mic_R=16) # manually specify breakpoints to define S/I/R

# note we need to set use_disk=FALSE otherwise the function will stop and tell us it can't find disk breakpoints
# and we need to set expected_R=FALSE, otherwise we assume it is expected R based on the EUCAST general guideline that Enterobacterales are expected resistant to macrolides
azi_rules <- amrrules_save(azi_analysis,
                           dir="amrrules",
                           file_prefix="Azithromycin",
                           mic_S=16, mic_R=16,
                           use_disk=FALSE,
                           expected_R=FALSE)

# update the `breakpoint standard` field to record where the S/R breakpoints came from
azi_rules$rules <- azi_rules$rules %>% mutate(`breakpoint standard`=if_else(`clinical category`!="-", "ECOFF 2025", "-"))

view(azi_rules$rules)

# compare predictions vs observed ECOFF calls
azi_rules$predict_vs_obs_stats$plot_ecoff
azi_rules$predict_vs_obs_stats_byMethod$plot_ecoff
azi_rules$predict_vs_mic_dist_byMethod$pred_ecoff
```

## Running a full analysis for a given species

To draft a full set of rules for a given species, it is suggested to summarise all the available ast + geno data, and use this to create a table of drugs to analyse, with all the necessary run parameters. Then use this parater table to make rules sequentially for each drug in the list, incrementing the ruleIDs for each new analysis so you end up with a set of rules files that have non-overlapping ruleIDs and can ultimately be concatenated to make a single set of proposed.

Suggested code to help do this:

### Summarise available data per drug
```
# data per drug
data_per_drug <- ast %>% 
  mutate(method = if_else(!is.na(disk), "Disk", method)) %>%
  group_by(drug_agent, method) %>% count() %>% 
  pivot_wider(id_cols=drug_agent, names_from=method, values_from = n, values_fill=0) %>%
  mutate(n=sum(c_across(where(is.numeric)), na.rm = TRUE)) %>%
  mutate(n_mic=sum(c_across(c(`BD Phoenix`:`microbroth dilution`, Microscan)))) %>%
  mutate(n_disk=Disk) %>%
  mutate(n_measure=sum(c_across(n_mic:n_disk))) %>%
  mutate(n_SIR=`NA`) %>%
  mutate()

# breakpoints per drug
breakpoints <- tibble()
for (drug in data_per_drug$drug_agent) {
  breakpoints <- bind_rows(breakpoints, getBreakpoints(species="Escherichia coli", antibiotic=drug, type_filter="human", guide="EUCAST 2025"))
  breakpoints <- bind_rows(breakpoints, getBreakpoints(species="Escherichia coli", antibiotic=drug, type_filter="human", guide="CLSI 2025"))
  breakpoints <- bind_rows(breakpoints, getBreakpoints(species="Escherichia coli", antibiotic=drug, type_filter="ECOFF"))
}

# drugs with different S or R breakpoints (i.e. for different sites)
multi_S <- breakpoints %>% group_by(guideline, ab, method) %>% summarise(s=length(unique(breakpoint_S))) %>% filter(s>1)
multi_R <- breakpoints %>% group_by(guideline, ab, method) %>% summarise(r=length(unique(breakpoint_R))) %>% filter(r>1)
  
# full set of breakpoints by site, method and guide
breakpoints_wide <- breakpoints %>% 
  mutate(breakpoint_description=paste(guideline, method, site)) %>% select(ab, breakpoint_description, breakpoint_S, breakpoint_R) %>%
  pivot_longer(cols=breakpoint_S:breakpoint_R) %>%
  mutate(breakpoint_description = paste(breakpoint_description, name)) %>% select(-name) %>%
  pivot_wider(names_from=breakpoint_description, values_from=value)

# simple breakpoints, first per guideline and method

breakpoints_simple_S_MIC <- breakpoints %>% filter(method=="MIC") %>%
    mutate(description=if_else(type=="ECOFF", "ECOFF MIC", paste(guideline, method, "S"))) %>%
    select(ab, description, breakpoint_S) %>% 
    pivot_wider(names_from=description, values_from=breakpoint_S, values_fn=~max(.x, na.rm=T))
  
breakpoints_simple_S_disk <- breakpoints %>% filter(method=="DISK") %>%
    mutate(description=if_else(type=="ECOFF", "ECOFF DISK", paste(guideline, method, "S"))) %>%
    select(ab, description, breakpoint_S) %>% 
    pivot_wider(names_from=description, values_from=breakpoint_S, values_fn=~min(.x, na.rm=T))

breakpoints_simple_R_MIC <- breakpoints %>% filter(method=="MIC" & type !="ECOFF") %>%
    mutate(description=paste(guideline, method, "R")) %>%
    select(ab, description, breakpoint_R) %>% 
    pivot_wider(names_from=description, values_from=breakpoint_R, values_fn=~max(.x, na.rm=T))
  
breakpoints_simple_R_disk <- breakpoints %>% filter(method=="DISK" & type !="ECOFF") %>%
    mutate(description=paste(guideline, method, "R")) %>%
    select(ab, description, breakpoint_R) %>% 
    pivot_wider(names_from=description, values_from=breakpoint_R, values_fn=~min(.x, na.rm=T))

data_per_drug_summary <- data_per_drug %>% 
  left_join(breakpoints_simple_S_MIC, join_by(drug_agent==ab)) %>%
  left_join(breakpoints_simple_R_MIC, join_by(drug_agent==ab)) %>%
  left_join(breakpoints_simple_S_disk, join_by(drug_agent==ab)) %>%
  left_join(breakpoints_simple_R_disk, join_by(drug_agent==ab)) %>%
  mutate(multiS=if_else(drug_agent %in% multi_S$ab, TRUE, FALSE)) %>%
  mutate(multiR=if_else(drug_agent %in% multi_R$ab, TRUE, FALSE)) %>%
  mutate(drug_name=ab_name(drug_agent)) %>%
  relocate(drug_name, .after=drug_agent)

data_per_drug_full <- left_join(data_per_drug_summary, breakpoints_wide, join_by(drug_agent==ab))

write_tsv(data_per_drug_full, file="data_per_drug_full.tsv")
```
### Generate a table of run parameters
Generate a table summarising the drugs and associated parameters to generate rules with, for all drugs that have EUCAST breakpoints and at least 1000 samples with assay measures, or drugs that we want to 'rescue' and include even though they have less data, because they have EUCAST breakpoints and we think there is enough data to at least explore.

``` {r}
# antibiotic, class, bp_site, use_mic, mic_S, mic_R, use_disk, disk_S, disk_R
# breakpoint site, mic_S, mic_R can be left blank (NA) if there is a single EUCAST breakpoint

# drugs to include even though they have too few samples, because they are clinically relevant (have EUCAST breakpoints)
rescue_list <- as.ab(c("TCC","TEM","CEC","CFM","CPT","CZA","IMR", "MEV","NOR","MFX","OFX"))

run_params_eucast <- data_per_drug_full %>% 
  filter(!(drug_agent %in% c("NAL", "FOX"))) %>% # exclude these screening drugs
  filter(n_measure>1000 | drug_agent %in% rescue_list) %>% # only considering drugs with at least 1000 MIC or disk measurements
  mutate(use_mic=if_else(n_mic>=20 & (!is.na(`EUCAST 2025 MIC S`) | !is.na(`ECOFF MIC`)), TRUE, FALSE)) %>% # if we have MIC data and ecoff or breakpoint
  mutate(mic_S=if_else(is.na(`EUCAST 2025 MIC S`) & !is.na(`ECOFF MIC`), `ECOFF MIC`, NA)) %>% # use ECOFF if no breakpoint
  mutate(mic_R=if_else(is.na(`EUCAST 2025 MIC R`) & !is.na(`ECOFF MIC`), `ECOFF MIC`, NA)) %>% # use ECOFF if no breakpoint
  mutate(use_disk=if_else(n_disk>=20 & (!is.na(`EUCAST 2025 DISK S`) | !is.na(`ECOFF DISK`)), TRUE, FALSE)) %>% # if we have disk data and ecoff or breakpoint
  mutate(disk_S=if_else(is.na(`EUCAST 2025 DISK S`) & !is.na(`ECOFF DISK`), `ECOFF DISK`, NA)) %>% # use ECOFF if no breakpoint
  mutate(disk_R=if_else(is.na(`EUCAST 2025 DISK R`) & !is.na(`ECOFF DISK`), `ECOFF DISK`, NA)) %>%
  filter(use_disk | use_mic) # exclude those with no data or no breakpoint

  
# for those with multiple (condition-specific) breakpoints, explore which we should use
breakpoints %>% filter(ab %in% multi_S$ab & guideline=="EUCAST 2025") %>% select(ab, method, site, breakpoint_S, breakpoint_R)

sites <- tibble(drug_agent=as.ab("CIP"), bp_site=c("Non-meningitis")) %>% # cipro, meningitis breakpoint is set to identify presence of resistance mechanisms, which we can achieve directly and through comparison to ECOFF
  bind_rows(tibble(drug_agent=as.ab("AMX"), bp_site=c("Intravenous"))) %>% # amoxicillin, other breakpoints are the same or equate to expected 'I'
  bind_rows(tibble(drug_agent=as.ab("AMC"), bp_site=c("Oral, Uncomplicated urinary tract infection"))) %>% # amoxi-clav, don't need iv breakpoint because same as ECOFF, don't need other oral breakpoints as they equate to expected 'I'
  bind_rows(tibble(drug_agent=as.ab("CXM"), bp_site=c("Uncomplicated urinary tract infection"))) # cefuroxime, don't need iv breakpoint as it amounts to expected 'I'

# make table of drugs to analyse, with their run parameters
run_params_eucast <- run_params_eucast %>% 
  select(drug_agent, use_mic, mic_S, mic_R, use_disk, disk_S, disk_R) %>% 
  left_join(sites) %>%
  mutate(mic_S=if_else(!is.na(bp_site), NA, mic_S)) %>% # don't specify the breakpoint if there are multiple defined
  mutate(mic_R=if_else(!is.na(bp_site), NA, mic_R)) %>%
  mutate(disk_S=if_else(!is.na(bp_site), NA, disk_S)) %>%
  mutate(disk_R=if_else(!is.na(bp_site), NA, disk_R))

> run_params_eucast
# A tibble: 49 × 8
# Groups:   drug_agent [49]
   drug_agent use_mic mic_S mic_R use_disk disk_S disk_R bp_site                                    
   <ab>       <lgl>   <dbl> <dbl> <lgl>     <dbl>  <dbl> <chr>                                      
 1 AMC        TRUE    NA    NA    TRUE         NA     NA Oral, Uncomplicated urinary tract infection
 2 AMK        TRUE    NA    NA    TRUE         NA     NA NA                                         
 3 AMP        TRUE    NA    NA    TRUE         NA     NA NA                                         
 4 AMX        TRUE    NA    NA    FALSE        NA     NA Intravenous                                
 5 ATM        TRUE    NA    NA    TRUE         NA     NA NA                                         
 6 AZM        TRUE    16    16    FALSE        NA     NA NA                                         
 7 BPR        TRUE     0.25  0.25 FALSE        NA     NA NA                                         
 8 CAZ        TRUE    NA    NA    TRUE         NA     NA NA                                         
 9 CCV        TRUE     0.5   0.5  FALSE        NA     NA NA                                         
10 CFM        TRUE    NA    NA    TRUE         NA     NA NA                                         
# ℹ 39 more rows
# ℹ Use `print(n = ...)` to see more rows
```

### Add the list of genotype classes to consider for each drug
Now we need to make a list of which AMRfinderplus classes to consider, for each drug. This has to be done manually currently, until we sort out a proper mapping of AMR package terms to AMRfp classes, and the fact that some drugs map to multiple classes (e.g. to analyse ampicillin phenotypes, you neeed to include markers in 3 classes: `Beta-lactams/penicillins`, `Carbapenems`, `Cephalosporins`. 

So for now it's best to do this manually and thoughtfully! You can use this table in the AMRgen package to help, but it may not be complete: `AMRgen::amrfp_drugs_table`

What you need to do is make a table that looks like this, for each drug in your `run_params_eucast` table:

```
> drug_to_amrfp_class
# A tibble: 74 × 5
   drug_agent drug_name                   AMRfp_class              class2      class3        
   <ab>       <chr>                       <chr>                    <chr>       <chr>         
 1 TOB        Tobramycin                  Aminoglycosides          NA          NA            
 2 GEN        Gentamicin                  Aminoglycosides          NA          NA            
 3 AMK        Amikacin                    Aminoglycosides          NA          NA            
 4 STR        Streptoduocin               Aminoglycosides          NA          NA            
 5 STR1       Streptomycin                Aminoglycosides          NA          NA            
 6 KAN        Kanamycin                   Aminoglycosides          NA          NA            
 7 CHL        Chloramphenicol             Amphenicols              NA          NA            
 8 AMP        Ampicillin                  Beta-lactams/penicillins Carbapenems Cephalosporins
 9 AMC        Amoxicillin/clavulanic acid Beta-lactams/penicillins Carbapenems Cephalosporins
10 TZP        Piperacillin/tazobactam     Beta-lactams/penicillins Carbapenems Cephalosporins
```

Then use this to make a list of classes, and add this list to you run_params_eucast table
```
drug_to_amrfp_class <- drug_to_amrfp_class %>% 
  mutate(classes=if_else(!is.na(class2), paste(AMRfp_class,class2,class3,sep=","), AMRfp_class))

run_params_eucast <- run_params_eucast %>% 
  left_join(drug_to_amrfp_class %>% select(drug_agent, classes))

```

The `makerules` function checks whether each drug has expected resistance, by checking the `AMR::intrinsic_resistant` table. If a drug is noted in this table, then any calls other than 'R' will be flagged for review in the generated rules table. To override this, you need to pass `expected_R=FALSE` to the `makerules` function, so record this parameter in the run params table:

```
# manually overrule AMR::intrinsic_resistant that Enterobacterales are expected R to azithromycin, so we can call S/R from data
run_params_eucast <- run_params_eucast %>% mutate(expected_R=if_else(drug_agent=="AZM", FALSE, NA))
```

### Set up fixed run parameters, and analyse all our drugs in one go
So now we have a table of all the drug-specific run parameters. 

Next let's set up all our fixed parameters, that will apply to all drugs:
```
# for geno/pheno analyses
species <- "Escherichia coli"
sir_col <- "pheno_eucast"
ecoff_col <- "ecoff"
geno_sample_col <- "Name"
pheno_sample_col <- "id"
minPPV <- 1
mafLogReg <- 5
mafUpset <- 1
marker_col <- "marker"

info_obj <- ast %>% select(id, source, method)

# for rule definitions
minObs <- 1
low_threshold <-20

# where should we start when numbering the new rules? make sure this is higher than the number of rules you currently have defined.
max_rule_id <- c(999)

# helper function to convert NA values from our run_params table to `NULL` values for our function call
na2null <- function(x) {
  if(is.na(x)) {return(NULL)}
  else{return(x)}
  }
```

Run analysis on all our rules. Note the code below uses the default settings for amrrules_save() which includes using the drafted rules to predict phenotypes. This step can take a long time if you have a lot of markers, if you don't need it then consider turning it off by setting `testRules=FALSE` in the amrrules_save() call.
```
for (i in 1:nrow(run_params_eucast)) {
  
  # run the analysis, including blaEC
  analysis <- amrrules_analysis(afp, ast, antibiotic=ab_name(run_params_eucast$drug_agent[i]),
                                drug_class_list=unlist(strsplit(run_params_eucast$classes[i], ",")),
                                species=species, sir_col=sir_col, ecoff_col=ecoff_col, 
                                geno_sample_col=geno_sample_col, 
                                pheno_sample_col=pheno_sample_col, 
                                minPPV=1, mafLogReg=5, mafUpset=1, info=info_obj,
                                use_mic=run_params_eucast$use_mic[i],
                                mic_S=na2null(run_params_eucast$mic_S[i]),
                                mic_R=na2null(run_params_eucast$mic_R[i]),
                                use_disk=run_params_eucast$use_disk[i],
                                disk_S=na2null(run_params_eucast$disk_S[i]),
                                disk_R=na2null(run_params_eucast$disk_R[i]),
                                call_manual=TRUE # recall SIR based on provided breakpoints
                                )
  
  # store the results and generate rules
  if (!is.na(run_params_eucast$bp_site[i])) {
    file_prefix <- paste0(gsub("/","_",ab_name(run_params_eucast$drug_agent[i])),"_",gsub(",","",run_params_eucast$bp_site[i]))
    }
  else{file_prefix=NULL}
  
  rules <- amrrules_save(analysis, dir_path="amrrulesEcoliFeb4", 
                         ruleID_start=(max_rule_id+1), 
                         minObs=minObs, low_threshold=low_threshold,
                         use_mic=run_params_eucast$use_mic[i],
                         mic_S=na2null(run_params_eucast$mic_S[i]),
                         mic_R=na2null(run_params_eucast$mic_R[i]),
                         use_disk=run_params_eucast$use_disk[i],
                         disk_S=na2null(run_params_eucast$disk_S[i]),
                         disk_R=na2null(run_params_eucast$disk_R[i]),
                         bp_site=na2null(run_params_eucast$bp_site[i]),
                         expected_R=na2null(run_params_eucast$expected_R[i]),
                         file_prefix=file_prefix)
  
  # update rule ID for next drug
  next_rule_id <- max(rules$rules$ruleID) %>% gsub(pattern="ECO", replacement="") %>% as.numeric()
  check_next_rule_id <- as.integer(next_rule_id)
  if(is.na(check_next_rule_id)) {max_rule_id <- max_rule_id+1000}
  else{max_rule_id = next_rule_id}
 
}
```

## Other considerations

### Core genes with WT S
If there are core genes in AMRfinderplus that definitely have no impact on phenotype, these will interfere with our ability to consider solo PPV for other markers in the same class, as those other markers will hardly ever be found 'solo' without the irrelevant core gene. E.g. in E. coli, blaEC is called in nearly every genome, but has no effect on ampicillin.

In this case, it makes sense to exclude the irrelevant core gene. Note this should only be done after you have confirmed the gene has no effect, even in combination with others, by running the analysis with it included first. But, if you are confident the gene is irrelevant, filter it out of your genotype file and then run the analysis on this filtered file:
```
afp_exclBlaEC <- afp %>% filter(!grepl("blaEC", gene))
```

### Impact of allelic variation amongst gene detection hits
For some genes where we are interested in presence/absence (i.e. detection), the AMRfinderplus database has multiple protein sequences that map to the same element symbol. E.g. refgene has 3 entries for [cmlA1](https://www.ncbi.nlm.nih.gov/pathogens/refgene/#cmlA1) which each have different reference protein sequences but map to the same node [cmlA1](https://www.ncbi.nlm.nih.gov/pathogens/genehierarchy/#node_id:cmlA1) in the gene hierarchy and are reported in AMRfinderplus results with the same element symbol `cmlA1`.
This means that where we see `cmlA1` in the `element symbol` (v4+) or `gene symbol` (pre-v4) field in AMRfinderplus results, (which via `import_amrfp` we represent as `node`='cmlA1', variation type='Gene presence detected', `marker.label`='cmlA1'), this may actually describe different allelic variants. Also, by default `import_amrfp` does not consider coverage and identity of hits, it just relies on AMRfinderplus to tell us if a hit is 'partial' (defined as <90% coverage) in which case the variation type will be called as 'Inactivating mutation detected' and the marker.label will be 'cmlA1:-' to indicate a partial hit.
By default, AMRgen and AMRrulemaker functions will use the 'marker.label' field as the unit for analysis, which is potentially grouping together several different allelic variants and comparing them as a group to the phenotype results (although excluding partial/broken hits). This will often be sensible, but it may be problematic if the different alleles actually have different functions.
To explore this, we can define a new 'marker.label' based on the closest reference sequence, and whether the match is identical or not to this, and then re-run the quantitative analysis and rule definitions using this alternative unit of analysis.

```
# define new marker labels
afp_accession <- afp %>% 
    mutate(exact_match=if_else(`% Coverage of reference sequence`==100 & `% Identity to reference sequence`==100.00, "", "x_")) %>%
    mutate(node_hit = paste0(node, "__", exact_match, `Accession of closest sequence`)) %>%
    mutate(marker.label = case_when(`variation type`=="Inactivating mutation detected" ~ paste0(node_hit,":-"),
                                    !is.na(mutation) ~ paste0(node_hit, ":", mutation),
                                    TRUE ~ node_hit
      
    ))

```
In this example in E. coli, we now have 5 different marker labels for hits to cmlA1, most are identical matches to the same reference protein (cmlA1__WP_000095725.1), but there are also inexact matches to this reference (cmlA1__x_WP_000095725.1), exact matches to an alternative allele (cmlA1__WP_001256776.1), and broken hits to both of these references (which would in the initial analysis have been grouped as 'cmlA1:-', as distinct from 'cmlA1'):
```
afp_accession %>% filter(grepl("cmlA1",marker.label)) %>% count(marker.label)
# A tibble: 5 × 2
  marker.label                  n
  <chr>                     <int>
1 cmlA1__WP_000095725.1       350
2 cmlA1__WP_001256776.1         8
3 cmlA1__x_WP_000095725.1      46
4 cmlA1__x_WP_000095725.1:-    12
5 cmlA1__x_WP_001256776.1:-    11
```

Running the chloramphenicol analysis on this new version of the genotyping table, we get separate PPVs for each of these types of hits, and see that while the evidence is strong for exact matches to the WP_000095725.1 reference sequence (positive predictive value 88.7% for R and 92.2% for NWT), neither of the 2 genomes that had exact matches to WP_001256776 and callable SIR phenotypes were called as R. Of the strains with inexact matches to WP_000095725.1 (cmlA1__x_WP_000095725.1), 3/5 were R and 4/6 were NWT. None of the genomes with broken hits were R, which is what we would expect.
```
> chloramphenicol %>% filter(grepl("cmlA1", marker))
# A tibble: 8 × 8
  marker                    category     x     n   ppv     se ci.lower ci.upper
  <chr>                     <chr>    <dbl> <dbl> <dbl>  <dbl>    <dbl>    <dbl>
1 cmlA1__WP_000095725.1     R           86    97 0.887 0.0322    0.823    0.950
2 cmlA1__WP_001256776.1     R            0     2 0     0         0        0    
3 cmlA1__x_WP_000095725.1   R            3     5 0.6   0.219     0.171    1    
4 cmlA1__x_WP_000095725.1:- R            0     1 0     0         0        0    
5 cmlA1__x_WP_001256776.1:- R            0     3 0     0         0        0    
6 cmlA1__WP_000095725.1     NWT         94   102 0.922 0.0266    0.869    0.974
7 cmlA1__WP_001256776.1     NWT          1     3 0.333 0.272     0        0.867
8 cmlA1__x_WP_000095725.1   NWT          4     6 0.667 0.192     0.289    1 
```
