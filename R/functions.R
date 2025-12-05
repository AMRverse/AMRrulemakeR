#' Summarise Genotypic and Phenotypic Data for a given Antibiotic
#'
#' This function creates a summary of genotype and phenotype information for a specified antibiotic
#' and species using provided genotypic and phenotypic data tables. It reports sample counts with
#' relevant genetic markers and susceptibility test results, and lists EUCAST breakpoints and ECOFFs.
#'
#' @param geno_table A data frame of genotypic data, including a column `drug_class` and sample IDs. Can be generated from AMRfinderplus output using `AMRgen::import_amrfp`
#' @param pheno_table A data frame of phenotypic data, including susceptibility results and sample IDs. Required fields are `drug_agent`, `mic`, `disk`.
#' @param antibiotic A string naming the antibiotic of interest (eg: \code{"Ciprofloxacin"}), will be parsed using AMR::as.ab().
#' @param drug_class_list A character vector of drug classes to include from the genotypic data (eg: \code{c("Quinolones")}).
#' @param geno_sample_col The name of the column in \code{geno_table} that contains sample IDs (default: "Name").
#' @param pheno_sample_col The name of the column in \code{pheno_table} that contains sample IDs (default: "id").
#' @param species The species for which EUCAST breakpoints should be retrieved (default: "E. coli"), will be parsed using AMR::as.mo().
#' @param guide The guidelines to retrieve breakpoints from using the AMR package (default: "EUCAST 2025").
#'
#' @return A tibble summarizing the number of samples with genotype and phenotype data, and listing EUCAST breakpoints and ECOFFs.
#'
#' @details
#' - Uses `getBreakpoints()` to pull breakpoint and ECOFF data for the specified species, antibiotic, and guide.
#' - Separately summarizes MIC and disk diffusion data if available.
#'
#' @examples
#' \dontrun{
#' summarise_data(geno_table, pheno_table, antibiotic = "Ciprofloxacin", drug_class_list=c("Quinolones"), species="E. coli")
#' }
#'
#' @export
summarise_data <- function(geno_table, pheno_table, antibiotic, drug_class_list,
                           geno_sample_col="Name", pheno_sample_col="id", species, guide="EUCAST 2025") {
  geno_samples <- geno_table %>% pull(get(geno_sample_col)) %>% unique()

  pheno_rows <- pheno_table %>% filter(drug_agent==as.ab(antibiotic))
  pheno_rows_mic <- pheno_rows %>% filter(!is.na(mic))
  pheno_rows_disk <- pheno_rows %>% filter(!is.na(disk))
  pheno_samples <- pheno_rows %>% pull(get(pheno_sample_col)) %>% unique()
  pheno_samples_mic <- pheno_rows_mic %>% pull(get(pheno_sample_col)) %>% unique()
  pheno_samples_disk <- pheno_rows_disk %>% pull(get(pheno_sample_col)) %>% unique()

  breakpoints_eucast <- getBreakpoints(species, guide, antibiotic, type_filter="human", guide="EUCAST")
  breakpoints_clsi <- getBreakpoints(species, guide, antibiotic, type_filter="human", guide="CLSI")
  ecoff <- getBreakpoints(species, guide, antibiotic, type_filter="ECOFF")

  summary <- rbind(c(paste("Samples with", antibiotic,"phenotypes:"), length(pheno_samples)),
                   c(" - MIC", length(pheno_samples_mic)),
                   c(" - DISK", length(pheno_samples_disk)),
                   c(paste("Samples with genotypes and phenotypes:"), sum(geno_samples %in% pheno_samples)),
                   c(" - MIC", sum(geno_samples %in% pheno_samples_mic)),
                   c(" - DISK", sum(geno_samples %in% pheno_samples_disk)),
                   c("EUCAST breakpoint sites:", nrow(breakpoints_eucast))
  )
  if (nrow(breakpoints_eucast)>0) {
    for (i in 1:nrow(breakpoints_eucast)) {
      summary <- rbind(summary,
                       c(paste0(" - ", breakpoints_eucast$method[i], " / ",breakpoints_eucast$site[i]),
                         paste0(breakpoints_eucast$breakpoint_S[i], ", ", breakpoints_eucast$breakpoint_R[i]))
      )
    }
  }

  if (nrow(ecoff)>0) {
    summary <- rbind(summary,c("EUCAST ECOFFs:", nrow(ecoff)))
    for (i in 1:nrow(ecoff)) {
      summary <- rbind(summary,
                       c(paste0(" - ", ecoff$method[i], " / ",ecoff$site[i]),
                         paste0(ecoff$breakpoint_S[i], ", ", ecoff$breakpoint_R[i]))
      )
    }
  }

  return(summary %>% as_tibble())
}




#' Get Clinical Breakpoints for an Antibiotic
#'
#' This function retrieves the clinical breakpoints for a given species, antibiotic, and guideline, from the AMR package.
#' It attempts to find the breakpoints at various taxonomic levels (species, genus, family, order) if no direct match is found.
#'
#' @param species A character string representing the species of interest (e.g., "Escherichia coli").
#' @param guide A character string indicating the guideline for breakpoints (default, "EUCAST 2024").
#' @param antibiotic A character string indicating the antibiotic for which to retrieve breakpoints (e.g., "Ciprofloxacin").
#' @param type_filter A character string indicating the type of breakpoints to retrieve (e.g., "human"). Default is "human" which returns human clinical breakpoints, change to "ECOFF" to get the epidemiological cutoff.
#'
#' @return A data frame containing the clinical breakpoints for the specified species, antibiotic, and guideline.
#' If no exact match is found, the function attempts to retrieve breakpoints at the genus, family, and order levels.
#'
#' @examples
#' getBreakpoints("Escherichia coli", "EUCAST 2024", "Ciprofloxacin")
#'
#' @export
getBreakpoints <- function(species, guide="EUCAST 2024", antibiotic, type_filter="human") {
  bp <- AMR::clinical_breakpoints %>% filter(guideline==guide & mo==AMR::as.mo(species) & ab==AMR::as.ab(antibiotic)) %>% filter(type==type_filter)
  if(nrow(bp)==0) {
    sp_mo <- AMR::as.mo(species)
    this_mo <- AMR::microorganisms %>% filter(mo==sp_mo)
    # try genus
    bp <- AMR::clinical_breakpoints %>% filter(guideline==guide & mo==AMR::as.mo(this_mo$genus) & ab==AMR::as.ab(antibiotic)) %>% filter(type==type_filter)
    if (nrow(bp)==0) {
      # try family
      bp <- AMR::clinical_breakpoints %>% filter(guideline==guide & mo==AMR::as.mo(this_mo$family) & ab==AMR::as.ab(antibiotic)) %>% filter(type==type_filter)
      if (nrow(bp)==0) {
        # try order
        bp <- AMR::clinical_breakpoints %>% filter(guideline==guide & mo==AMR::as.mo(this_mo$order) & ab==AMR::as.ab(antibiotic)) %>% filter(type==type_filter)
      }
    }
  }
  return(bp)
}

#' Check and Retrieve Breakpoints for an Antibiotic
#'
#' This function checks the clinical breakpoints for a specified antibiotic and species using the `getBreakpoints` function.
#' It handles cases where multiple breakpoint sites exist and uses the specified site or the one with the highest susceptibility
#' breakpoint.
#'
#' @param species A character string representing the species of interest (e.g., "Escherichia coli").
#' @param guide A character string indicating the guideline for breakpoints (default, "EUCAST 2024").
#' @param antibiotic A character string indicating the antibiotic for which to check breakpoints (e.g., "Ciprofloxacin").
#' @param bp_site A character string specifying the breakpoint site to use (optional). If provided, the function uses this site; otherwise, if different breakpoints are specified for different sites it selects the one with the highest susceptibility breakpoint.
#' @param assay A character string specifying the assay type (either "MIC" or "Disk"). Default is "MIC".
#'
#' @return A list containing:
#' - \code{breakpoint_S}: The susceptibility breakpoint (e.g., MIC or disk size).
#' - \code{breakpoint_R}: The resistance breakpoint (e.g., MIC or disk size).
#' - \code{bp_standard}: The breakpoint site used, if multiple breakpoints are set in the guidelines.
#'
#' @details
#' The function first attempts to retrieve breakpoints using the `getBreakpoints` function. If multiple breakpoint sites are found, it handles the situation by:
#' - Using the specified site if it exists.
#' - Selecting the breakpoint with the highest susceptibility value if the specified site is not found.
#' - Returning a message about the selected site and breakpoint values.
#'
#' @examples
#' checkBreakpoints(species="Escherichia coli", guide="EUCAST 2024", antibiotic="Ciprofloxacin", assay="MIC")
#'
#' @export
checkBreakpoints <- function(species, guide="EUCAST 2024", antibiotic, assay="MIC", bp_site=NULL) {
  breakpoints <- getBreakpoints(species, guide, antibiotic) %>% filter(method==assay)
  if (nrow(breakpoints)==0) {stop(paste("Could not determine",assay,"breakpoints using AMR package, please provide your own breakpoints"))}
  else{
    breakpoint_sites <- unique(breakpoints$site)
    breakpoint_message_multibp <- NA
    bp_standard="-"
    # handle multiple breakpoints (e.g. for different conditions)
    if (length(breakpoint_sites)>1) {
      breakpoint_message_multibp <- paste("NOTE: Multiple breakpoint entries, for different sites:", paste(breakpoint_sites, collapse="; "))
      if (length(unique(breakpoints$breakpoint_R))==1 & length(unique(breakpoints$breakpoint_S))==1) {
        breakpoints <- breakpoints %>% arrange(-breakpoint_S) %>% first()
        breakpoint_message_multibp <- paste0(breakpoint_message_multibp, ". However S and R breakpoints are the same.")
        bp_standard<-breakpoints$site
      }
      else if (is.null(bp_site)) {
        breakpoints <- breakpoints %>% arrange(-breakpoint_S) %>% first()
        breakpoint_message_multibp <- paste0(breakpoint_message_multibp, ". Using the one with the highest S breakpoint (", breakpoints$site,").")
        bp_standard<-breakpoints$site
      }
      else {
        if (bp_site %in% breakpoint_sites) {
          breakpoints <- breakpoints %>% filter(site==bp_site)
          breakpoint_message_multibp <- paste0(breakpoint_message_multibp, ". Using the specified site (", bp_site,").")
          bp_standard<-breakpoints$site
        }
        else {
          breakpoints <- breakpoints %>% arrange(-breakpoint_S) %>% first()
          breakpoint_message_multibp <- paste0(breakpoint_message_multibp, ". Could not find the specified site (", bp_site,"), so using the one with the highest S breakpoint (", breakpoints$site,").")
          bp_standard <- breakpoints$site
        }
      }
    }
    breakpoint_S <- breakpoints$breakpoint_S
    breakpoint_R <- breakpoints$breakpoint_R
    if (is.na(breakpoint_S) | is.na(breakpoint_R)) {stop(paste("Could not determine",assay,"breakpoints using AMR package, please provide your own breakpoints"))}
    if (assay=="MIC") { breakpoint_message <- paste("MIC breakpoints determined using AMR package: S <=", breakpoint_S,"and R >", breakpoint_R) }
    else { breakpoint_message <- paste("Disk diffusion breakpoints determined using AMR package: S >=", breakpoint_S,"and R <", breakpoint_R) }
    cat(paste0("  ",breakpoint_message, "\n"))
    if(!is.na(breakpoint_message_multibp)) {cat(paste0("  ",breakpoint_message_multibp, "\n"))}
  }
  return(list(breakpoint_S=breakpoint_S,breakpoint_R=breakpoint_R, bp_standard=bp_standard))
}




# Helper functions
safe_execute <- function(expr) {
  tryCatch({
    expr
  }, error = function(e) {
    message("Error in executing command: ", e$message)
    return(NULL)
  })
}

na0 <- function(x) {
  ifelse(is.na(x), 0, x)
}

add_missing_cols <- function(df, cols) {
  for (col in cols) {
    if (!col %in% names(df)) {
      df <- df %>% mutate(!!sym(col) := NA)
    }
  }
  df
}

compareSets <- function(v, strings, split_delim) {
  matches <- NULL
  for (s in strings) {
    v2 <- unlist(str_split(s, split_delim))
    if(all(v2 %in% v)) {
      matches <- c(matches,s)
    }
  }
  return(matches)
}

# gene = string defining a combination rule ID
# ruleIDs = vector of ruleIDs
compare2sets <- function(gene, ruleIDs){
  component_ruleIDs <- unlist(str_split(gene, " & "))
  return(all(component_ruleIDs %in% ruleIDs))
}

getGenes <- function(data, combo) {
  marker_names <- unlist(str_split(combo, ", "))
  if (length(marker_names)>1) {
    ids <- data %>% filter(marker %in% marker_names) %>% pull(ruleID)
    return(paste(ids, collapse=" & "))
  }
  else{ return(combo)}
}

getSources <- function(amr_binary, combo, assay, exclude_range_values=FALSE) {

  # assume combo list is in proper syntax (with ':') but amr_binary has colnames with '..' in place of ':'
  marker_names <- unlist(str_split(gsub(":", "..", combo), ", "))

  dat <- amr_binary %>%
    filter(!is.na(get(assay))) 
    
  if (exclude_range_values) {
    dat <- dat %>% filter(!grepl("<|>", as.character(mic)))
  }
  
  total_sources <- dat %>%
    #select(source, all_of(marker_names)) %>%
    #filter(if_all(all_of(marker_names), ~ . == 1)) %>%
    select(source, any_of(marker_names)) %>%
    filter(if_all(any_of(marker_names), ~ . == 1)) %>%
    distinct() %>% nrow()

  sources_per_pheno <- dat %>%
    #select(source, pheno, all_of(marker_names)) %>%
    #filter(if_all(all_of(marker_names), ~ . == 1)) %>%
    select(source, pheno, any_of(marker_names)) %>%
    filter(if_all(any_of(marker_names), ~ . == 1)) %>%
    distinct() %>%
    group_by(pheno) %>% count(name="sources") %>%
    pivot_wider(names_from=pheno, names_prefix="sources.", values_from=sources) %>%
    mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))

  sources <- as_tibble(c(sources=total_sources, sources_per_pheno))

  source_names <- c("sources", "sources.S", "sources.I", "sources.R")
  
  if (exclude_range_values) {
    prefix <- paste0(assay,".x.")
  } else {prefix <- paste0(assay,".")}
  
  sources <- add_missing_cols(sources, source_names) %>%
    rename_with(.fn = ~ paste0(prefix, .x))
  
  return(sources)
}

enumerate_source_info <- function(data, info, solo_binary, amr_binary, column="source", use_mic=NULL, use_disk=NULL) {
  if (column %in% colnames(info)) {
    solo_sources_pheno <- solo_binary %>%
      left_join(info, by="id") %>%
      filter(!is.na(pheno) & !is.na(marker) & value==1) %>%
      select(all_of(column), marker) %>% distinct() %>%
      group_by(marker) %>% count(name=paste0("solo.",column,"s_pheno"))
    
    solo_sources_ecoff <- solo_binary %>%
      left_join(info, by="id") %>%
      filter(!is.na(ecoff) & !is.na(marker) & value==1) %>%
      select(all_of(column), marker) %>% distinct() %>%
      group_by(marker) %>% count(name=paste0("solo.",column,"s_ecoff"))
    
    # in case we have markers where only solo observations have an ecoff call or pheno call but not both
    solo_sources <- full_join(solo_sources_pheno, solo_sources_ecoff, by="marker") %>% 
      ungroup() %>% rowwise() %>%
      mutate(solo.sources = max(c_across(c(2, 3)), na.rm = TRUE)) %>%
      ungroup() %>%
      select(marker, solo.sources)
    
    # add sources per pheno value
    solo_sources_per_pheno <- solo_binary %>%
      left_join(info, by="id") %>%
      filter(!is.na(pheno) & !is.na(marker) & value==1) %>%
      select(all_of(column), marker, pheno) %>% distinct() %>%
      group_by(marker, pheno) %>% count() %>%
      pivot_wider(names_from=pheno, names_prefix=paste0("pheno_solo.",column,"s."), values_from=n, values_fill=0)
    
    solo_sources_per_ecoff <- solo_binary %>%
      left_join(info, by="id") %>%
      filter(!is.na(ecoff) & !is.na(marker) & value==1) %>%
      select(all_of(column), marker, ecoff) %>% distinct() %>%
      group_by(marker, ecoff) %>% count() %>%
      pivot_wider(names_from=ecoff, names_prefix=paste0("ecoff_solo.",column,"s."), values_from=n, values_fill=0) 
    
    solo_sources_per <- full_join(solo_sources_per_pheno, solo_sources_per_ecoff, by="marker") %>% 
      ungroup() %>% rowwise() %>%
      mutate(solo.sources.R=suppressWarnings(max(c_across(matches("solo\\.sources.R")), na.rm = TRUE))) %>% 
      mutate(solo.sources.R=if_else(is.infinite(solo.sources.R),0,solo.sources.R)) %>%
      mutate(solo.sources.I=suppressWarnings(max(c_across(matches("solo\\.sources.I")), na.rm = TRUE))) %>% 
      mutate(solo.sources.I=if_else(is.infinite(solo.sources.I),0,solo.sources.I)) %>%
      mutate(solo.sources.S=suppressWarnings(max(c_across(matches("solo\\.sources.S")), na.rm = TRUE))) %>% 
      mutate(solo.sources.S=if_else(is.infinite(solo.sources.S),0,solo.sources.S)) %>%
      ungroup() %>%
      select(any_of(c("marker", "solo.sources.R", "solo.sources.I", "solo.sources.S")))
      
    solo_sources <- full_join(solo_sources, solo_sources_per, by="marker")
    
    data <- data %>% left_join(solo_sources, by="marker")
    
    amr_binary <- amr_binary %>% left_join(info, by="id")
    
    if (use_mic) {
      data <- data %>% rowwise() %>% 
        mutate(source_info=getSources(amr_binary, marker, "mic")) %>% 
        unnest_wider(source_info) %>% ungroup()
      
      data <- data %>% rowwise() %>% 
        mutate(source_info=getSources(amr_binary, marker, "mic", exclude_range_values=TRUE)) %>% 
        unnest_wider(source_info) %>% ungroup()
    }
    
    if (use_disk) {
      data <- data %>% rowwise() %>% 
        mutate(source_info=getSources(amr_binary, marker, "disk")) %>% 
        unnest_wider(source_info) %>% ungroup()
    }
    
    source_names <- c(paste0(column,"s"), paste0(column,"s.S"), paste0(column,"s.I"), paste0(column,"s.R"))
    
    data <- add_missing_cols(data, c(paste0("solo.",source_names), paste0("disk.",source_names), paste0("mic.",source_names))) %>%
      mutate(solo.sources.SIR = if_else(!is.na(solo.sources),
                                        paste0(na0(solo.sources.S), " S, ", na0(solo.sources.I), " I, ", na0(solo.sources.R), " R"),
                                        "-")) %>%
      mutate(mic.sources.SIR = if_else(!is.na(mic.sources) & use_mic,
                                       paste0(na0(mic.sources.S), " S, ", na0(mic.sources.I), " I, ", na0(mic.sources.R), " R"),
                                       "-")) %>%
      mutate(mic.x.sources.SIR = if_else(!is.na(mic.x.sources) & use_mic,
                                       paste0(na0(mic.x.sources.S), " S, ", na0(mic.x.sources.I), " I, ", na0(mic.x.sources.R), " R"),
                                       "-")) %>%
      mutate(disk.sources.SIR = if_else(!is.na(disk.sources) & use_disk,
                                        paste0(na0(disk.sources.S), " S, ", na0(disk.sources.I), " I, ", na0(disk.sources.R), " R"),
                                        "-"))
  }
  return(data)
}


categorise_from_PPV <- function(data, suffix="", use_mic=F, mic_S=NULL, mic_R=NULL, use_disk=F, disk_S=NULL, disk_R=NULL) {
  nwt_ppv=grep("^NWT\\..*\\.ppv$", colnames(data), value = TRUE)
  if(length(nwt_ppv)==0) {nwt_ppv=""}
  r_ppv=grep("^R\\..*\\.ppv$", colnames(data), value = TRUE)
  if(length(r_ppv)==0) {r_ppv=""}
  
  if (nwt_ppv %in% colnames(data)) {
    # call pheno from NWT PPV (i.e. based on ecoff)
    data <- data %>% mutate(phenotype=if_else(get(nwt_ppv)>0.5, "NWT", "WT"))
    if (r_ppv %in% colnames(data)) {
      # call category using R & NWT PPV (update to include I PPV if available)
      data <- data %>% 
        mutate(category = case_when(get(r_ppv)>=0.8 ~ "R",
                                    get(nwt_ppv)<=0.5 ~ "S",
                                    get(nwt_ppv)>0.5 & (is.na(get(r_ppv)) | get(r_ppv)<0.8) ~ "I",
                                    is.na(get(nwt_ppv)) & get(r_ppv)<=0.5 ~ "S",
                                    is.na(get(nwt_ppv)) & get(r_ppv)>0.5 ~ "I",
                                    TRUE ~ NA))
      if (!is.null(mic_S) & !is.null(mic_R) & use_mic) { 
        if (mic_S==mic_R) {# no I category for this drug
          data <- data %>% mutate(category=if_else(category=="I", "R", category))
        }
      }
      else if (!is.null(disk_S) & !is.null(disk_R) & use_disk) { # no I category for this drug
        if (disk_S==disk_R) {# no I category for this drug
          data <- data %>% mutate(if_else(category=="I", "R", category))
        }
      }
    } else { # NWT PPV but no R PPV, infer category from NWT PPV only (ie based on ecoff)
      data <- data %>% mutate(category = case_when(get(nwt_ppv)<=0.5 ~ "S",
                                                           get(nwt_ppv)>0.5 ~ "R",
                                                           is.na(get(nwt_ppv)) ~ NA,
                                                           TRUE ~ NA))
    }
  } else if (r_ppv %in% colnames(data)) { 
    # R PPV but no NWT PPV; call pheno from R PPV (update to I PPV if we have one)
    data <- data %>% mutate(phenotype=if_else(get(r_ppv)>0.5, "NWT", "WT"))
  }
  
  # always return with category_suffix and phenotype_suffix columns added with notes
  pheno_colname <- paste0("phenotype",suffix)
  if ("phenotype" %in% colnames(data)) {
    data <- data %>% rename(!!pheno_colname := phenotype)
  } else {data <- data %>% mutate(!!pheno_colname := NA)}

  category_colname <- paste0("category",suffix)
  if ("category" %in% colnames(data)) {
    data <- data %>% rename(!!category_colname := category)
  } else {data <- data %>% mutate(!!category_colname := NA)}
  
  return(data)
}


compare_categories <- function(category_soloPPV, category_soloExtPPV, category_micPPV, category_diskPPV, category_logReg, category_logRegExt, R.solo.n, R.soloExt.n, MIC.n, Disk.n, minObs, low_threshold) {
  
  category <- NA
  grade <- ""
  limitations <- list()
  note <- list()
  breakpoints_source <- list()
  
  if (!is.na(category_soloPPV)) {
    category <- category_soloPPV
    
    if(R.solo.n < low_threshold) {
      grade <- "low"
      limitations <- list("Limited samples.")
    }
    else {grade <- "moderate"}
    
    # ignore any missing category calls and record MIC or disk breakpoint if these measures contributed 
    if (is.na(category_soloExtPPV)) {category_soloExtPPV <- category}
    if (is.na(category_micPPV)) {category_micPPV <- category
    } else if (category_micPPV==category_soloPPV) {breakpoints_source <- append(breakpoints_source, "mic")}
    if (is.na(category_diskPPV)) {category_diskPPV <- category
    } else if (category_diskPPV==category_soloPPV) {breakpoints_source <- append(breakpoints_source, "disk")}
    
    if (category_soloPPV!=category_soloExtPPV) {
      if (R.solo.n<minObs & R.soloExt.n>=minObs) { # grade already set to low
        category <- category_soloExtPPV
        limitations <- append(limitations, list("Limited assay data."))
        note <- append(note, list("Call is based on solo PPV data from extended data set including SIR calls without assay measurements."))
      } else if (R.solo.n>=minObs) {
        limitations <- append(limitations, list("Conflicting evidence."))
        note <- append(note, list("Call is based on solo PPV data from combined SIR from MIC/disk measures, but conflicts with extended data set including SIR calls without assay measurements."))
      }
    } else if (category_soloPPV!=category_micPPV) {
      if (MIC.n>=minObs) {
        limitations <- append(limitations, list("Conflicting evidence."))
        note <- append(note, list("Call is based on solo PPV data based on combined SIR from MIC/disk measures, but conflicts with call from MIC measurements alone."))
      }
    } else if (category_soloPPV!=category_diskPPV) {
      if (Disk.n>=minObs) {
        limitations <- append(limitations, list("Conflicting evidence."))
        note <- append(note, list("Call is based on solo PPV data based on combined SIR from MIC/disk measures, but conflicts with call from disk measurements alone."))
      }
    }
  }
  else { # no solo PPV, should mean this is a combination (which only has MIC.ppv and/or Disk.ppv data) or marker/s found only in extended dataset without assay measures
    if (!is.na(category_micPPV) & !is.na(category_diskPPV)) {
      if (category_micPPV==category_diskPPV) { 
        category <- category_micPPV
        breakpoints_source <- append(breakpoints_source, "mic")
        breakpoints_source <- append(breakpoints_source, "disk")
        if ((MIC.n + Disk.n) < low_threshold) {
          grade <- "low"
          limitations <- append(limitations, list("Limited samples."))
        } else {grade <- "moderate"}
      } else if(MIC.n > Disk.n) { # more mic data than disk, use mic call
        category <- category_micPPV
        breakpoints_source <- append(breakpoints_source, "mic")
        if (MIC.n <- low_threshold) {
          grade <- "low"
          limitations <- append(limitations, list("Limited samples."))
        }  else {grade <- "moderate"}
        limitations <- append(limitations, list("Conflicting evidence."))
        note <- append(note, list("Using MIC call, disk data disagrees."))
      } else if(MIC.n < Disk.n) { # more mic data than disk, use mic call
        category <- category_diskPPV
        breakpoints_source <- append(breakpoints_source, "disk")
        if (Disk.n <- low_threshold) {
          grade <- "low"
          limitations <- append(limitations, list("Limited samples."))
        } else {grade <- "moderate"}
        limitations <- append(limitations, list("Conflicting evidence."))
        note <- append(note, list("Using disk call, MIC data disagrees."))
      }
    } else if (!is.na(category_micPPV)) {
      category <- category_micPPV
      breakpoints_source <- append(breakpoints_source, "mic")
      if (MIC.n <- low_threshold) {
        grade <- "low"
        limitations <- append(limitations, list("Limited samples."))
      } else {grade <- "moderate"}
    } else if (!is.na(category_diskPPV)) {
      category <- category_diskPPV
      breakpoints_source <- append(breakpoints_source, "disk")
      if (Disk.n <- low_threshold) {
        grade <- "low"
        limitations <- append(limitations, list("Limited samples."))
      } else {grade <- "moderate"}
    } else if (!is.na(category_soloExtPPV)) {
      category <- category_soloExtPPV
      grade <- "low"
      limitations <- append(limitations, list("Limited assay data."))
      if (R.soloExt.n <- low_threshold) {
        limitations <- append(limitations, list("Limited samples."))
      }
    } else if (!is.na(category_logReg)) {
      category <- category_logReg
      grade <- "low"
      limitations <- append(limitations, list("No solo data."))
      note <- append(note, list("Call based only on negative result in multiple logistic regression, based on combined SIR from MIC/disk measures."))
    } else if (!is.na(category_logRegExt)) {
      category <- category_logRegExt
      grade <- "low"
      limitations <- append(limitations, list("No solo data."))
      note <- append(note, list("Call based only on negative result in multiple logistic regression, using extended data set including SIR calls without assay measurements."))
    }
  }
  return(list(category=category, 
              limitations_list=limitations, 
              grade=grade, 
              breakpoints_source=breakpoints_source, 
              note=note))
}

breakpoints <- function(breakpoint_mic, breakpoint_disk, bp_source_list) {
  bp_list <- list()
  if ("mic" %in% bp_source_list) {
    bp_list <- append(bp_list, breakpoint_mic)
  }
  if ("disk" %in% bp_source_list) {
    bp_list <- append(bp_list, breakpoint_disk)
  }
  bp_string <- paste0(unlist(bp_list), collapse=", ")
  return(bp_string)
}

