# interpret AMRfinderplus genotypes using a set of rules generated from quantitative data
# output = phenotype call and category call per strain
# note this function is designed  for internal validation/exploration of rule sets, 
#   NOT annotating each line of a genotype report with an interpretation
#   OR for generating a detailed AMR report

test_rules_amrfp <- function(geno_table, rules, species) {
  
  # convert mutations to AMRrules format, store concatenation as label
  # NOTE this will not work for manually defined rules which use a hierarchy node, it is just for checking rules defined from quantitative geno-pheno analysis of AMRfp data
  geno_table <- geno_table %>%
    separate(`Gene symbol`, into = c("gene", "mutation"), sep = "_", remove=F, fill="right") %>%
    mutate(mutation=if_else(`Element subtype`=="POINT", mutation, NA)) %>%
    mutate(mutation=if_else(!is.na(mutation), convert_mutation(mutation), "-")) %>%
    mutate(label=if_else(`Element subtype`=="POINT" & !is.na(mutation), paste(gene, mutation), paste(`Gene symbol`, "-")))
                                 
  rules <- rules %>% mutate(label=paste(nodeID, mutation)) %>%
    filter(organism==species) %>%
    mutate(`clinical category`=factor(`clinical category`, levels=c("S", "I", "R"), ordered=T)) %>%
    mutate(phenotype=factor(phenotype, levels=c("wildtype", "nonwildtype"), ordered=T))
  
  calls <- geno_table %>% group_by(Name) %>% 
    summarise(label_list=paste(label, collapse=",")) %>%
    #ungroup() %>%
    mutate(call_info = map(label_list, ~ getCall(.x, rules = rules))) %>%
    unnest(call_info, keep_empty = TRUE)

  return(calls)
}


getCall <- function(label_list, rules) {
  # get individual rule IDs
  label_list <- unlist(str_split(label_list, ","))
  matching_rules <- rules %>% filter(label %in% label_list)
  
  # if there's more than one, check for combination rules
  if (length(matching_rules)>1) { 
    ruleIDs <- matching_rules %>% pull(ruleID)
    matching_rules <- rules %>% filter(grepl("&", gene)) %>%
      rowwise() %>%
      filter(compare2sets(gene, ruleIDs)) %>% 
      bind_rows(matching_rules)
  }
  
  if (nrow(matching_rules)>0) {
    call <- matching_rules %>% 
      mutate(label=sub(" -","",label)) %>% 
      group_by(drug) %>% 
      summarise(phenotype=max(phenotype, na.rm=T), category=max(`clinical category`, na.rm=T),
              rules=paste(ruleID, collapse=",")) %>%
      select(drug, phenotype, category, rules)
  }
  else {call <- tibble(drug = NA_character_,
                       phenotype = NA_character_,
                       category = NA_character_,
                       rules = NA_character_)}
  return(call)
}
