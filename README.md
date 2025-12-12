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

# example data from AMRgen package: E. coli MIC data from NCBI, matching AMRfinderplus data
ecoli_ast
ecoli_geno <- import_amrfp(ecoli_geno_raw, "Name")

# run quantitative analyses
cip_analysis <- amrrules_analysis(ecoli_geno, ecoli_ast, antibiotic="Cipro", drug_class_list=c("Quinolones"), species="E. coli")

cip_analysis$ppv_plot
cip_analysis$upset_mic_plot
cip_analysis$upset_disk_plot

# save tables and plots and generate rules
cip_rules <- amrrules_save(cip_analysis, bp_site="Non-meningitis", dir_path="amrrules", file_prefix="Cipro", use_disk=F, guide="CLSI 2025")

# alternatively, call makerules directly on the analysis object
cip_rules <- makerules(cip_analysis, bp_site="Non-meningitis", use_disk=F, guide="CLSI 2025")

# view the proposed rules, in AMRrules specification format, with quantitative fields added
view(cip_rules$rules)

```

# Work in progress - suggested protocol for developing AMRrules using this package

## Collate and format phenotype data
For use with the AMRgen & AMRrulemakeR packages, you need to get the AST data into long format (one row per sample/drug result), with some key fields.
Data in NCBI or EBI antibiogram format can be automatically imported to the right dataframe using the `import_ast()` in the AMRgen package.
Or you can format your data manually. 

Example data frame included in the AMRgen package:

```
> ecoli_ast

# A tibble: 4,170 × 10
   id           drug_agent     mic  disk pheno_clsi ecoff guideline method pheno_provided spp_pheno   
   <chr>        <ab>         <mic> <dsk> <sir>      <sir> <chr>     <chr>  <sir>          <mo>        
 1 SAMN36015110 CIP        <128.00    NA   R          R   CLSI      NA       R            B_ESCHR_COLI
 2 SAMN11638310 CIP         256.00    NA   R          R   CLSI      NA       R            B_ESCHR_COLI
 3 SAMN05729964 CIP          64.00    NA   R          R   CLSI      Etest    R            B_ESCHR_COLI
 4 SAMN10620111 CIP         >=4.00    NA   R          R   CLSI      NA       R            B_ESCHR_COLI
 5 SAMN10620168 CIP         >=4.00    NA   R          R   CLSI      NA       R            B_ESCHR_COLI
 6 SAMN10620104 CIP         <=0.25    NA   S          R   CLSI      NA       S            B_ESCHR_COLI

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
For use with the AMRrulemakeR package, you need to process the AMRfinderplus data using `import_afp()`, which generates consistent marker labels that will be the units of analysis, and which can ultimately be represented in the AMRrules variant specification format.

Example data frame included in the AMRgen package:

```
ecoli_geno <- import_amrfp(ecoli_geno_raw, "Name")

ecoli_geno

# A tibble: 50,720 × 36
   Name         gene        mutation  node      `variation type` marker marker.label drug_agent drug_class
   <chr>        <chr>       <chr>     <chr>     <chr>            <chr>  <chr>        <ab>       <chr>     
 1 SAMN03177615 blaEC       NA        blaEC     Gene presence d… blaEC  blaEC        NA         Beta-lact…
 2 SAMN03177615 acrF        NA        acrF      Gene presence d… acrF   acrF         NA         Efflux    
 3 SAMN03177615 glpT        Glu448Lys glpT      Protein variant… glpT_… glpT:Glu448… FOS        Other ant…
 4 SAMN03177615 floR        NA        floR      Gene presence d… floR   floR         CHL        Amphenico…
 5 SAMN03177615 floR        NA        floR      Gene presence d… floR   floR         FLR        Other ant…
 6 SAMN03177615 mdtM        NA        mdtM      Gene presence d… mdtM   mdtM         NA         Efflux    
 7 SAMN03177615 blaTEM-1    NA        blaTEM-1  Gene presence d… blaTE… blaTEM-1     NA         Beta-lact…
 8 SAMN03177615 sul2        NA        sul2      Gene presence d… sul2   sul2         SSS        Other ant…
 9 SAMN03177615 aph(3'')-Ib NA        aph(3'')… Gene presence d… aph(3… aph(3'')-Ib  STR1       Aminoglyc…
10 SAMN03177615 aph(6)-Id   NA        aph(6)-Id Gene presence d… aph(6… aph(6)-Id    STR1       Aminoglyc…

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
# process AMRfp genotype data from Allthebacteria

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
```

## Run analyses need to define rules
The function `amrrules_analysis()` takes our phenotype table and extracts the data for a specified drug; then takes our genotype table and extracts the data for the relevant markers; and compares these geno/pheno data using several different `AMRgen` functions.

The function `amrrules_save()` takes these analysis results, and saves key tables and figures. It then uses these results to define rules in AMRrules format via the `makerules()` function, and writes these to an output TSV file. Finally, the rules are applied to the input genotypes to predict S/I/R and wildtype/nonwildtype for each sample using the `test_rules_amrfp()` function, and these results are compared to the input phenotypes and summarised in terms of positive predictive value (including stratified by assay method) to help assess the validity of the rules.


Example command
```
# import phenotype data in the right format
afp <-read_tsv("ATB_AFP_Ecoli_AMR.tsv.gz")

# import AMRfinderplus data in the right format
ast <- import_ebi_ast("EBI_CABBAGE_EcoliShigella.csv.gz")

# extract the information fields we want to consider in the analyses (source, method)
info_obj <- ast %>% select(id, source, method)

# run the required analyses to compare phenotypes for a specific drug (e.g.ciprofloxacin)
# with genetic markers associated with the corresponding drug class (e.g. quinolones)
analysis <- amrrules_analysis(geno_table=afp, pheno_table=ast,
                              antibiotic="Ciprofloxacin",
                              drug_class_list=c("Quinolones"),
                              species="Escherichia coli",
                              sir_col="pheno_eucast", ecoff_col="ecoff", 
                              minPPV=1, mafLogReg=5, mafUpset=1,
                              info=info_obj)

# use the results of these analyses to define rules, then apply the rules back to the data to predict phenotypes
rules <- amrrules_save(analysis, dir_path="amrrules", 
                         bp_site="Non-meningitis",
                         ruleID_start=1,
                         use_mic=TRUE, use_disk=TRUE,
                         geno_table=afp, pheno_table=ast,
                         file_prefix=file_prefix)
```

Notes:
1. Currently we use the supplied pheno_table to get the assay method to plot predictions vs method, but we could do this from the info_obj.
2. Add examples with multiple entries in drug_class_list (e.g. trim-sulfa, beta-lactams)
3. Add examples with manually supplied breakpoints (e.g. using azithromycin ECOFF as breakpoint)
