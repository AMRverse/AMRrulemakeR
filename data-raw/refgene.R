refgene <- read_tsv("ReferenceGeneCatalog_v3.12_2024-01-31.1.txt.gz")
hierarchy <- read_tsv("ReferenceGeneHierarchy_v3.12_2024-01-31.1.txt.gz")

refgene_pubmed <- read_tsv("ReferenceGeneCatalog_v3.12_2024-01-31.1.txt.gz",
                     col_types = cols(allele=col_character(), gene_family=col_character(),
                                      whitelisted_taxa=col_character(), pubmed_reference=col_character(),
                                      .default = col_skip()))

usethis::use_data(refgene, internal = FALSE, overwrite = TRUE)
usethis::use_data(refgene_pubmed, internal = FALSE, overwrite = TRUE)
usethis::use_data(hierarchy, internal = FALSE, overwrite = TRUE)
