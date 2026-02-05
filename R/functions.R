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
      df <- df %>% mutate(!!sym(col) := 0)
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

enumerate_source_info <- function(data, info=NULL, solo_binary, amr_binary, column="source", use_mic=NULL, use_disk=NULL) {
  if (!(column %in% colnames(solo_binary))) {
    if (!is.null(info)) {
      if (column %in% colnames(info)) {
        solo_binary <- solo_binary %>%
          left_join(info %>% select(id, !!sym(column)), by="id")
      }
    }
  }

  if (!(column %in% colnames(amr_binary))) {
    if (!is.null(info)) {
      if (column %in% colnames(info)) {
        amr_binary <- amr_binary %>%
          left_join(info %>% select(id, !!sym(column)), by="id")
      }
    }
  }

  if (column %in% colnames(solo_binary)) {
    solo_sources_pheno <- solo_binary %>%
      filter(!is.na(pheno) & !is.na(marker) & value==1) %>%
      select(all_of(column), marker) %>% distinct() %>%
      group_by(marker) %>% count(name=paste0("solo.",column,"s_pheno"))

    solo_sources_ecoff <- solo_binary %>%
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
      filter(!is.na(pheno) & !is.na(marker) & value==1) %>%
      select(all_of(column), marker, pheno) %>% distinct() %>%
      group_by(marker, pheno) %>% count() %>%
      pivot_wider(names_from=pheno, names_prefix=paste0("pheno_solo.",column,"s."), values_from=n, values_fill=0)

    solo_sources_per_ecoff <- solo_binary %>%
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
  }

  if (column %in% colnames(amr_binary)) {
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
                                        "-"))
    if ("mic.sources" %in% colnames(data) & use_mic) {
      data <- data %>% mutate(mic.sources.SIR = paste0(na0(mic.sources.S), " S, ", na0(mic.sources.I), " I, ", na0(mic.sources.R)))
    } else {data$mic.sources.SIR <- "-"}

    if ("mic.x.sources" %in% colnames(data) & use_mic) {
      data <- data %>% mutate(mic.x.sources.SIR = paste0(na0(mic.x.sources.S), " S, ", na0(mic.x.sources.I), " I, ", na0(mic.x.sources.R), " R"))
    } else {data$mic.x.sources.SIR <- "-"}

    if ("disk.sources" %in% colnames(data) & use_disk) {
      data <- data %>% mutate(disk.sources.SIR = paste0(na0(disk.sources.S), " S, ", na0(disk.sources.I), " I, ", na0(disk.sources.R), " R"))
    } else {data$disk.sources.SIR <- "-"}
  }

  return(data)
}

#' Classify markers as S/I/R and WT/NWT based on solo PPV values
#'
#' This function takes as input a dataframe with columns indicating PPV values for R, I and/or NWT, and returns the data frame with new fields indicating the SIR and NWT/WT calls made from the PPV estimates.
#'
#' @param data A data frame with column names including PPV estimates, named in the form of 'R.x.ppv', 'I.x.ppv', 'NWT.x.ppv'. This is intended as a helper function, called within the makerules() function.
#' @param suffix A character string indicating the suffix to append to the field names for the generated category (SIR) and phenotype (NWT/WT) fields.
#'
#' @return A copy of the input data frame with new fields indicating the SIR category and NWT/WT phenotype calls made from the PPV estimates.
#'
#' @examples
#' \dontrun{
#' data <- solo_stats %>%
#'  pivot_wider(id_cols=marker, names_from=category, values_from=c(ppv, x, n, ci.lower, ci.upper),
#'            names_glue= "{category}.solo.{.value}") %>%
#'  categorise_from_PPV(suffix="_soloPPV")
#' }
#'
#' @export
categorise_from_PPV <- function(data, suffix="") {

  nwt_ppv=grep("^NWT\\..*\\.ppv$", colnames(data), value = TRUE)
  if(length(nwt_ppv)==0) {nwt_ppv=""}

  r_ppv=grep("^R\\..*\\.ppv$", colnames(data), value = TRUE)
  if(length(r_ppv)==0) {r_ppv=""}

  i_ppv=grep("^I\\..*\\.ppv$", colnames(data), value = TRUE)
  if(length(i_ppv)==0) {i_ppv=""}

  if (nwt_ppv %in% colnames(data)) {
    # call pheno from NWT PPV (i.e. based on ecoff)
    data <- data %>% mutate(phenotype=if_else(get(nwt_ppv)>0.5, "NWT", "WT"))
  } else if (i_ppv %in% colnames(data)) {
    # I PPV but no NWT PPV; call pheno from I PPV
    data <- data %>% mutate(phenotype=if_else(get(i_ppv)>0.5, "NWT", "WT"))
  } else if (r_ppv %in% colnames(data)) {
    # R PPV but no NWT PPV; call pheno from R PPV
    data <- data %>% mutate(phenotype=if_else(get(r_ppv)>0.5, "NWT", "WT"))
  }

  if (r_ppv %in% colnames(data)) {
    # call category using R
    data <- data %>% mutate(category = if_else(get(r_ppv)>0.5, "R", "S"))
    if (i_ppv %in% colnames(data)) {
      data <- data %>% mutate(category = if_else(get(i_ppv)>0.5 & get(r_ppv)<0.8, "I", category))
    }
  } else if (nwt_ppv %in% colnames(data)) { # NWT PPV but no R PPV, infer category from NWT PPV only (ie based on ecoff)
    data <- data %>% mutate(category = case_when(get(nwt_ppv)<=0.5 ~ "S",
                                                 get(nwt_ppv)>0.5 ~ "R",
                                                 is.na(get(nwt_ppv)) ~ NA,
                                                 TRUE ~ NA))
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

#' Compare category calls and sample sizes for different sources of PPV (solo MIC+disk, extended solo, MIC, disk, logistic regression) and adjudicate a final call, recording evidence grade, limitations and notes for AMRrules. Helper function used by makerules().
#'
#' @param category_soloPPV Call from solo PPV based on combined interpretation of MIC+disk assay data.
#' @param category_soloExtPPV Call from solo PPV based on all available SIR calls including those without raw assay data.
#' @param category_micPPV Call from solo PPV based on MIC data.
#' @param category_diskPPV Call from solo PPV based on disk data.
#' @param category_logReg Call from multivariable logistic regression including all markers exceeding a minimum frequency, for samples with assay measures.
#' @param category_logRegExt Call from multivariable logistic regression including all markers exceeding a minimum frequency, for all samples with SIR calls including those without raw assay data.
#' @param R.solo.n Number of samples contributing to category_soloPPV call.
#' @param R.soloExt.n Number of samples contributing to category_soloExtPPV call.
#' @param MIC.n Number of samples contributing to category_micPPV call.
#' @param Disk.n Number of samples contributing to category_diskPPV call.
#' @param minObs Minimum number of samples to make a call.
#' @param low_threshold Minimum number of observations required for `evidence grade` to be assigned as "moderate" rather than "low". Default is 20.
#'
#' @return A list of values
#'
#' @examples
#' \dontrun{
#' data <- data %>%
#'   rowwise() %>%
#'   mutate(results=list(compare_categories(category_soloPPV=category_soloPPV,
#'              category_soloExtPPV=category_soloExtPPV,
#'              category_micPPV=category_micPPV,
#'              category_diskPPV=category_diskPPV,
#'              category_logReg=category_logReg,
#'              category_logRegExt=category_logRegExt,
#'              R.solo.n=R.solo.n, R.soloExt.n=R.soloExt.n,
#'              MIC.n=MIC.n, Disk.n=Disk.n))) %>%
#'  ungroup() %>%
#'  unnest_wider(results)
#' }
#'
#' @export
compare_categories <- function(category_soloPPV, category_soloExtPPV, category_micPPV,
                               category_diskPPV, category_logReg, category_logRegExt,
                               R.solo.n, R.soloExt.n, MIC.n, Disk.n,
                               minObs=1, low_threshold=20) {

  category <- NA
  grade <- ""
  limitations <- list()
  note <- list("")
  breakpoints_source <- list()

  if (!is.na(category_soloPPV)) {
    category <- category_soloPPV

    if (R.solo.n < low_threshold) {
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
      if (!is.na(MIC.n)) {
        if (MIC.n>=minObs) {
          limitations <- append(limitations, list("Conflicting evidence."))
          note <- append(note, list("Call is based on solo PPV data based on combined SIR from MIC/disk measures, but conflicts with call from MIC measurements alone."))
        }
      }
    } else if (category_soloPPV!=category_diskPPV) {
      if (!is.na(Disk.n)) {
        if (Disk.n>=minObs) {
          limitations <- append(limitations, list("Conflicting evidence."))
          note <- append(note, list("Call is based on solo PPV data based on combined SIR from MIC/disk measures, but conflicts with call from disk measurements alone."))
        }
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
      breakpoints_source <- append(breakpoints_source, "logreg")
    } else if (!is.na(category_logRegExt)) {
      category <- category_logRegExt
      grade <- "low"
      limitations <- append(limitations, list("No solo data."))
      note <- append(note, list("Call based only on negative result in multiple logistic regression, using extended data set including SIR calls without assay measurements."))
      breakpoints_source <- append(breakpoints_source, "logreg")
    }
  }
  return(list(category=category,
              limitations_list=limitations,
              grade=grade,
              breakpoints_source=breakpoints_source,
              note=note))
}

breakpoints <- function(breakpoint_mic, breakpoint_disk, bp_source_list, mic.sources, disk.sources) {
  bp_list <- list()
  if ("mic" %in% bp_source_list) {
    bp_list <- append(bp_list, breakpoint_mic)
  }
  if ("disk" %in% bp_source_list) {
    bp_list <- append(bp_list, breakpoint_disk)
  }
  if ("logreg" %in% bp_source_list) {
    if (mic.sources >0) {bp_list <- append(bp_list, breakpoint_mic)}
    if (disk.sources >0) {bp_list <- append(bp_list, breakpoint_disk)}
  }
  bp_string <- paste0(unlist(bp_list), collapse=", ")
  return(bp_string)
}

# check MICs expressed as ranges, how these are interpreted and if it impacts PPV estimates
# plot MIC distributions for isolates with no resistance markers, stratified by platform/method to help check the impact of range values in different methods
checkMICranges <- function(geno_table, pheno_table, antibiotic, drug_class_list, sir_col="pheno_eucast", ecoff_col="ecoff", minPPV=1, marker_col="marker.label", icat=TRUE, species, bp_site=NULL, excludeRanges="NWT") {

  pheno_table_micdisk <- pheno_table %>% filter(!is.na(mic) | !is.na(disk))

  soloPPV_micdisk <- safe_execute(AMRgen::solo_ppv_analysis(geno_table=geno_table, pheno_table=pheno_table_micdisk, antibiotic=antibiotic, drug_class_list=drug_class_list, sir_col=sir_col, ecoff_col=ecoff_col, min=minPPV, marker_col=marker_col, icat=icat, excludeRanges = NULL))

  soloPPV_micdisk_norange <- safe_execute(AMRgen::solo_ppv_analysis(geno_table=geno_table, pheno_table=pheno_table_micdisk, antibiotic=antibiotic, drug_class_list=drug_class_list, sir_col=sir_col, ecoff_col=ecoff_col, min=minPPV, marker_col=marker_col, icat=icat, excludeRanges = c("NWT", "R", "I")))

  compare_solo <- soloPPV_micdisk_norange$solo_stats %>% full_join(soloPPV_micdisk$solo_stats, by=c("marker", "category"), suffix=c("", ".ranges"))

  compare_solo_plot <- compare_solo %>%
    ggplot(aes(x=ppv, y=ppv.ranges, col=category)) +
    geom_abline(intercept=0, slope=1, col="grey", linetype=2) +
    geom_point() +
    labs(x="Excluding MICs expressed as >x", y="Including all MIC values")

  if (excludeRanges=="all") {
    cat("Returning solo stats ignoring MICs expressed as ranges\n")
    soloPPV_micdisk <- soloPPV_micdisk_norange # return analysis ignoring range values for R, I and NWT PPVs
  } else if (excludeRanges=="NWT") { # redo analysis ignoring range values for NWT only
    cat("Returning solo stats that ignore MICs expressed as ranges for NWT only (I/R stats still include range values)\n")
    soloPPV_micdisk <- safe_execute(AMRgen::solo_ppv_analysis(geno_table=geno_table, pheno_table=pheno_table_micdisk, antibiotic=antibiotic, drug_class_list=drug_class_list,
                                                              sir_col=sir_col, ecoff_col=ecoff_col, min=minPPV, marker_col=marker_col, icat=icat, excludeRanges = c("NWT")))

  }

  marker_count <- soloPPV_micdisk$amr_binary %>% select(-any_of(c("id", "pheno", "ecoff", "mic", "disk", "R", "I", "NWT"))) %>% rowSums()

  mic_plot_all <- assay_by_var(pheno_table=pheno_table_micdisk, antibiotic=antibiotic, measure="mic", facet_var="method",
                           species=species, bp_site=bp_site)

  marker_free_strains <- soloPPV_micdisk$amr_binary$id[marker_count==0]
  mic_plot_nomarkers <- assay_by_var(pheno_table_micdisk %>% filter(id %in% marker_free_strains),
                                     antibiotic=antibiotic, measure="mic", facet_var="method",
                                     species=species, bp_site=bp_site)

  return(list(mic_plot_nomarkers=mic_plot_nomarkers,
              mic_plot_all=mic_plot_all,
              plot=compare_solo_plot,
              table=compare_solo,
              soloPPV_micdisk=soloPPV_micdisk))
}

compare_interpretations <- function(pred, obs, antibiotic, sir_col="pheno_eucast", ecoff_col="ecoff", var="method",
                                    mic_S=NULL, mic_R=NULL, mic_ecoff=NULL, disk_S=NULL, disk_R=NULL, disk_ecoff=NULL) {
  antibiotic <- as.ab(antibiotic)
  obs <- obs %>% filter(!is.na(get(sir_col)) | !is.na(get(ecoff_col)))

  if (antibiotic %in% as.ab(pred$drug)) {
    pred <- pred %>% mutate(drug=as.ab(drug)) %>% filter(drug == antibiotic)
  } else { stop(paste("Antibiotic", antibiotic, "not found in predictions input"))}
  if (antibiotic %in% as.ab(obs$drug_agent)) {
    obs <- obs %>% filter(drug_agent == antibiotic)
  } else { stop(paste("Antibiotic", antibiotic, "not found in observations input"))}

  true_vs_predict <- left_join(pred, obs, join_by("Name"=="id"))

  ppv <- safe_execute(predictive_value(true_vs_predict, sir_col=sir_col, ecoff_col=ecoff_col))

  true_vs_predict <- true_vs_predict %>%
    mutate(type=case_when(!is.na(mic) ~ "mic", !is.na(disk) ~ "disk", TRUE ~ "none"))
  if (var %in% colnames(true_vs_predict)) {
    true_vs_predict <- true_vs_predict %>% mutate(type=paste(type, !!sym(var)))
  }

  ppv_bymethod <- safe_execute(predictive_value_by_var(true_vs_predict, sir_col=sir_col, ecoff_col=ecoff_col, var="type"))

  dist_mic_bypred <- safe_execute(dist_by_pred(true_vs_predict %>% filter(!is.na(mic)), antibiotic, assay="mic", var=NULL,
                                               breakpoint_ecoff=mic_ecoff, breakpoint_S=mic_S, breakpoint_R = mic_R))
  dist_mic_bypred_bymethod <- safe_execute(dist_by_pred(true_vs_predict %>% filter(!is.na(mic)), antibiotic, assay="mic", var=var,
                                              breakpoint_ecoff=mic_ecoff, breakpoint_S=mic_S, breakpoint_R = mic_R))

  dist_disk_bypred <- safe_execute(dist_by_pred(true_vs_predict %>% filter(!is.na(disk)), antibiotic, assay="disk", var=NULL,
                                                breakpoint_ecoff=disk_ecoff, breakpoint_S=disk_S, breakpoint_R = disk_R))
  dist_disk_bypred_bymethod <- safe_execute(dist_by_pred(true_vs_predict %>% filter(!is.na(disk)), antibiotic, assay="disk", var=var,
                                                breakpoint_ecoff=disk_ecoff, breakpoint_S=disk_S, breakpoint_R = disk_R))

  return(list(true_vs_predict=true_vs_predict,
         pred_ppv=ppv,
         pred_ppv_bymethod=ppv_bymethod,
         dist_mic_bypred=dist_mic_bypred,
         dist_disk_bypred=dist_disk_bypred,
         dist_mic_bypred_bymethod=dist_mic_bypred_bymethod,
         dist_disk_bypred_bymethod=dist_disk_bypred_bymethod))
}

predictive_value <- function(true_vs_predict, sir_col="pheno_eucast", ecoff_col="ecoff") {
    sir_ppv <- true_vs_predict %>%
      filter(!!sym(sir_col) %in% c("S", "I", "R")) %>%
      count(!!sym(sir_col), category) %>%
      group_by(category) %>% mutate(pct= prop.table(n) * 100) %>%
      mutate(category=as.sir(category)) %>% arrange(category) %>%
      group_by(category) %>% mutate(total=sum(n)) %>%
      mutate(pred=paste0(category, " (n=",total,")"))
    sir_ppv_plot <- sir_ppv %>%
      ggplot(aes(fill=!!sym(sir_col), x=pred, y=n)) +
      geom_col(position="fill") + scale_fill_sir() + coord_flip() +
      labs(fill="Observed", x="Predicted", y="Proportion", subtitle="Proportional") +
      geom_text(aes(label=paste0(sprintf("%1.1f", pct),"%")), position=position_fill(vjust=0.5))
    sir_ppv_n_plot <- sir_ppv %>%
      ggplot(aes(fill=!!sym(sir_col), x=pred, y=n)) +
      geom_col() + scale_fill_sir() + coord_flip() +
      labs(fill="Observed", x="Predicted", y="Count", subtitle="Counts")
    test_sir_ppv_plot <- sir_ppv_n_plot + sir_ppv_plot +
      patchwork::plot_layout(guides="collect", axes="collect") +
      patchwork::plot_annotation("Predicted vs observed (S/I/R calls)")

    ecoff_ppv <- true_vs_predict %>%
      filter(!!sym(ecoff_col) %in% c("S", "R")) %>%
      count(!!sym(ecoff_col), phenotype) %>%
      group_by(phenotype) %>% mutate(pct= prop.table(n) * 100) %>%
      arrange(phenotype) %>%
      group_by(phenotype) %>% mutate(total=sum(n)) %>%
      mutate(pred=paste0(phenotype, " (n=",total,")"))
    ecoff_ppv_plot <- ecoff_ppv %>%
      ggplot(aes(fill=!!sym(ecoff_col), x=pred, y=n)) +
      geom_col(position="fill") + scale_fill_sir() + coord_flip() +labs(fill="Observed", x="Predicted", y="Proportion", subtitle="Proportional") +
      geom_text(aes(label=paste0(sprintf("%1.1f", pct),"%")), position=position_fill(vjust=0.5))
    ecoff_ppv_n_plot <- ecoff_ppv %>%
      ggplot(aes(fill=!!sym(ecoff_col), x=pred, y=n)) +
      geom_col() + scale_fill_sir() + coord_flip() +
      labs(fill="Observed", x="Predicted", y="Count", subtitle="Counts")
    test_ecoff_ppv_plot <- ecoff_ppv_n_plot + ecoff_ppv_plot +
      patchwork::plot_layout(guides="collect", axes="collect") +
      patchwork::plot_annotation("Predicted vs observed (ECOFF calls)")

    return(list(table_sir=sir_ppv, plot_sir=test_sir_ppv_plot,
                table_ecoff=ecoff_ppv, plot_ecoff=test_ecoff_ppv_plot))
}

predictive_value_by_var <- function(true_vs_predict, sir_col="pheno_eucast", ecoff_col="ecoff", var="method") {
  sir_ppv <- true_vs_predict %>%
    filter(!!sym(sir_col) %in% c("S", "I", "R")) %>%
    count(!!sym(sir_col), category, !!sym(var)) %>%
    group_by(category, !!sym(var)) %>% mutate(pct= prop.table(n) * 100) %>%
    mutate(category=as.sir(category)) %>% arrange(category)
  totals <- sir_ppv %>%
    group_by(!!sym(var), category) %>%
    summarise(n_total = sum(n, na.rm = TRUE), .groups = 'drop')

  sir_ppv_plot <- sir_ppv %>%
    ggplot(aes(fill=!!sym(sir_col), x=!!sym(var), y=n)) +
    geom_col(position="fill") + scale_fill_sir() + coord_flip() +
    labs(fill="Observed", x="Assay method", y="Proportion", subtitle="Proportional") +
    geom_text(aes(label=paste0(sprintf("%1.1f", pct),"%")), position=position_fill(vjust=0.5)) +
    facet_wrap(~category, ncol=1)
  sir_ppv_n_plot <- sir_ppv %>%
    ggplot(aes(fill=!!sym(sir_col), x=!!sym(var), y=n)) +
    geom_col() + scale_fill_sir() + coord_flip() +
    labs(fill="Observed", x="Assay method", y="Count", subtitle="Counts") +
    geom_text(data=totals, aes(x = !!sym(var), y=n_total, label=n_total),
      inherit.aes = FALSE, hjust = -0.2,size = 3.5) +
    facet_wrap(~category, ncol=1)
  test_sir_ppv_plot <- sir_ppv_n_plot + sir_ppv_plot +
    patchwork::plot_layout(guides="collect", axes="collect") +
    patchwork::plot_annotation("Predicted vs observed (S/I/R calls), by method")

  ecoff_ppv <- true_vs_predict %>%
    filter(!!sym(ecoff_col) %in% c("S", "R")) %>%
    count(!!sym(ecoff_col), phenotype, !!sym(var)) %>%
    group_by(phenotype, !!sym(var)) %>% mutate(pct= prop.table(n) * 100) %>%
    arrange(phenotype)
  totals_ecoff <- ecoff_ppv %>%
    group_by(!!sym(var), phenotype) %>%
    summarise(n_total = sum(n, na.rm = TRUE), .groups = 'drop')

  ecoff_ppv_plot <- ecoff_ppv %>%
    ggplot(aes(fill=!!sym(ecoff_col), x=!!sym(var), y=n)) +
    geom_col(position="fill") + scale_fill_sir() + coord_flip() +
    labs(fill="Observed", x="Assay method", y="Proportion", subtitle="Proportional") +
    geom_text(aes(label=paste0(sprintf("%1.1f", pct),"%")), position=position_fill(vjust=0.5)) +
    facet_wrap(~phenotype, ncol=1)
  ecoff_ppv_n_plot <- ecoff_ppv %>%
    ggplot(aes(fill=!!sym(ecoff_col), x=!!sym(var), y=n)) +
    geom_col() + scale_fill_sir() + coord_flip() +
    labs(fill="Observed", x="Assay method", y="Count", subtitle="Counts") +
    geom_text(data=totals_ecoff, aes(x = !!sym(var), y=n_total, label=n_total),
              inherit.aes = FALSE, hjust = -0.2,size = 3.5) +
    facet_wrap(~phenotype, ncol=1)
  test_ecoff_ppv_plot <- ecoff_ppv_n_plot + ecoff_ppv_plot +
    patchwork::plot_layout(guides="collect", axes="collect") +
    patchwork::plot_annotation("Predicted vs observed (ECOFF calls), by method")

  return(list(table_sir=sir_ppv, plot_sir=test_sir_ppv_plot,
              table_ecoff=ecoff_ppv, plot_ecoff=test_ecoff_ppv_plot))
}

dist_by_pred <- function(true_vs_predict, antibiotic, assay="mic", var="method",
                              breakpoint_S=NULL, breakpoint_R=NULL, breakpoint_ecoff=NULL) {

  if (assay %in% assay & assay %in% colnames(true_vs_predict)) {
    true_vs_predict <- true_vs_predict %>%
      filter(!is.na(!!sym(assay))) %>%
      mutate(category=as.sir(category))

    pred_data <- true_vs_predict %>% count(factor(!!sym(assay)), category)
    colnames(pred_data)[1] <- assay

    if(!is.null(breakpoint_ecoff)) {ecoff_subtitle <- paste0("ECOFF: ", breakpoint_ecoff)}
    else {ecoff_subtitle=""}

    if (!is.null(breakpoint_S) | !is.null(breakpoint_R) | !is.null(breakpoint_ecoff)) {
      subtitle <- list()
      if(!is.null(breakpoint_S)) {subtitle <- append(subtitle, paste0("S: ", breakpoint_S))}
      if(!is.null(breakpoint_R)) {subtitle <- append(subtitle, paste0("R: ", breakpoint_R))}
      subtitle <- paste(subtitle, collapse=",")

      # get x coordinates for breakpoints after converting to factor
      assay_order <- true_vs_predict %>% count(factor(!!sym(assay)))
      colnames(assay_order)[1] <- assay
      breakpoint_S <- c(1:nrow(assay_order))[assay_order[,1]==breakpoint_S]
      breakpoint_R <- c(1:nrow(assay_order))[assay_order[,1]==breakpoint_R]
      breakpoint_ecoff <- c(1:nrow(assay_order))[assay_order[,1]==breakpoint_ecoff]
    } else {subtitle=""}

    pred <- true_vs_predict %>%
      ggplot(aes(x=factor(!!sym(assay)), fill=category)) +
      geom_bar() + scale_fill_sir() +
      labs(x="Measure", y="count", fill="Predicted",
           title=paste(ab_name(as.ab(antibiotic)), "assay distributions for all samples"),
           subtitle=paste("coloured by predicted S/I/R", subtitle)) +
      theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))

    pred_ecoff <- true_vs_predict %>%
      ggplot(aes(x=factor(!!sym(assay)), fill=phenotype)) +
      geom_bar() +
      labs(x="Measure", y="count", fill="Predicted (ECOFF)",
           title=paste(ab_name(as.ab(antibiotic)), "assay distributions for all samples"),
           subtitle=paste("coloured by predicted WT/NWT", ecoff_subtitle)) +
      theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1)) +
      scale_fill_manual(values=c(wildtype="#3CAEA3", nonwildtype="#ED553B", `NA`="grey"))

    # add breakpoints
    if (!is.null(breakpoint_S)) {
      pred <- pred + geom_vline(xintercept=breakpoint_S, color="#3CAEA3", linetype=2) }
    if (!is.null(breakpoint_R)) {
      pred <- pred + geom_vline(xintercept=breakpoint_R, color="#ED553B", linetype=2) }
    if (!is.null(breakpoint_ecoff)) {
      pred_ecoff <- pred_ecoff + geom_vline(xintercept=breakpoint_ecoff, color="#ED553B", linetype=2) }

    # add facets
    if (!is.null(var)) {
      if (var %in% colnames(true_vs_predict)) {
        if (true_vs_predict %>% filter(!is.na(get(var))) %>% nrow() > 0) {
          pred <- pred + facet_wrap(~get(var), ncol=1, scales="free_y")
          pred_ecoff <- pred_ecoff + facet_wrap(~get(var), ncol=1, scales="free_y")
        }}
    }
  } else {
    pred <- NULL
    pred_ecoff <- NULL
  }
  return(list(pred=pred, pred_ecoff=pred_ecoff))
}
