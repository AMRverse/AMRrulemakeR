ecoli_ast_ebi <- import_ebi_ast("EBI_AST_ecoli.csv.gz", interpret_clsi = TRUE, interpret_eucast = TRUE, interpret_ecoff = TRUE) %>% filter(id %in% afp$Name)

ecoli_afp_atb <- afp %>% filter(Name %in% ecoli_ast_ebi$id) %>% filter(drug_class %in% c("Aminoglycosides","Quinolones","Trimethoprims","Sulfonamides", "Macrolides", "Macrolides/lincosamides", "Beta-lactams/penicillins", "Carbapenems", "Cephalosporins"))

usethis::use_data(ecoli_ast_ebi, internal = FALSE, overwrite = TRUE)
usethis::use_data(ecoli_afp_atb, internal = FALSE, overwrite = TRUE)


