ecoli_pheno_ebi <- download_ebi(genus="Escherichia", reformat=T, interpret_clsi = TRUE, interpret_eucast = TRUE, interpret_ecoff = TRUE)

ecoli_pheno_ebi <- ecoli_pheno_ebi %>%
  filter(drug %in% c("AMK", "GEN", "TMP", "CIP", "SUL", "AZM", "AMP", "MEM"))

# don't use genotypes downloaded direct from EBI as they lack a of of necessary fields from the AMRfinderplus output
# ecoli_geno_ebi <- download_ebi(data="genotype", genus="Escherichia", reformat=T)

# use full AMRfp genotype results sourced previously from Allthebacteria
ecoli_geno_ebi <- import_amrfp("data-raw/atb_ecoli.tsv.gz") %>% filter(id %in% ecoli_pheno_ebi$id)

# filter phenotypes to only those with ATB genotypes
ecoli_pheno_ebi <- ecoli_pheno_ebi %>% filter(id %in% ecoli_geno_ebi$id)

usethis::use_data(ecoli_pheno_ebi, internal = FALSE, overwrite = TRUE)
usethis::use_data(ecoli_geno_ebi, internal = FALSE, overwrite = TRUE)
