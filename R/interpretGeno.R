#' Interpret AMRFinderPlus Genotypes Using Quantitative Rules
#'
#' This function interprets AMRFinderPlus genotype data using a set of rules generated from quantitative data,
#' producing phenotype and clinical category calls for each strain. It is intended for internal validation and 
#' exploration of rule sets derived from genotype-phenotype correlation analyses (i.e. using `makerules`). 
#' It is **not** suitable for line-by-line genotype annotation or for generating complete AMR reports.
#'
#' @param geno_table A data frame containing AMRFinderPlus genotype output. Must include columns 
#'   `Gene symbol`, `Element subtype`, and `Name`.
#' @param rules A data frame of interpretation rules, typically derived from quantitative genotype-phenotype 
#'   associations. Must include columns `nodeID`, `mutation`, `organism`, `clinical category`, and `phenotype`.
#' @param species A character string specifying the species of interest. Only rules for this species will be applied.
#'
#' @return A data frame with phenotype and clinical category calls for each strain (`Name`). 
#'   The result includes concatenated genotype labels and the interpretation results from applying the rules.
#'
#' @details 
#' Mutations are reformatted to match the expected format used in rule definitions. 
#' This function assumes that rules are flat (i.e., not hierarchical). It is not compatible with manually 
#' defined rule sets that use hierarchical `nodeID`s.
#'
#' @seealso \code{\link{convert_mutation}}, \code{\link{getCall}}
#'
#' @examples
#' \dontrun{
#' # apply the rules to generate S/I/R and wildtype/nonwildtype calls
#' calls <- test_rules_amrfp(geno_table = my_genotypes, rules = my_rules, species = "Escherichia coli")
#' 
#' # compare these calls to the AST data phenotypes in a separate dataframe, `pheno_table` with SIR phenotypes in `pheno`
#' calls_vs_pheno <- calls %>% left_join(pheno_table, join_by(Name==id))
#' calls_vs_pheno %>% group_by(pheno, category) %>% count() %>% filter(pheno %in% c("S", "I", "R")) 
#' }
#'
#' @export
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
