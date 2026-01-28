# download https://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip
# extract file names.dmp, name columns as per readme.txt
nodes <- read_delim("nodes.dmp.gz", delim="\t|\t", col_names=F)
bacteria <- nodes %>% filter(X5==0) %>% select(X1,X3) %>% rename(tax_id=X1, rank=X3)

taxid <- read_delim("names.dmp.gz", delim="\t|\t", col_names=F)
colnames(taxid) <- c("tax_id", "name_txt", "unique name", "name class")
taxid <- taxid %>% mutate(`name class`= sub("\t\\|", "", `name class`))

taxid_bacteria <- taxid %>%
  filter(tax_id %in% bacteria$tax_id) %>%
  filter(`name class`=="scientific name") %>%
  select(tax_id, name_txt) %>%
  left_join(bacteria) %>%
  filter(rank %in% c("family", "genus", "species", "species group", "species subgroup", "subspecies", "subgenus", "subfamily"))

usethis::use_data(taxid_bacteria, internal = FALSE, overwrite = TRUE)
