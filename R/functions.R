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
#' @param guide The guidelines to retrieve breakpoints from using the AMR package (default: "EUCAST 2024").
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
                           geno_sample_col="Name", pheno_sample_col="id", species, guide="EUCAST 2024") {
  geno_rows <- geno_table %>% filter(drug_class %in% drug_class_list)
  geno_samples <- geno_rows %>% pull(get(geno_sample_col)) %>% unique()

  pheno_rows <- pheno_table %>% filter(drug_agent==as.ab(antibiotic))
  pheno_rows_mic <- pheno_rows %>% filter(!is.na(mic))
  pheno_rows_disk <- pheno_rows %>% filter(!is.na(disk))
  pheno_samples <- pheno_rows %>% pull(get(pheno_sample_col)) %>% unique()
  pheno_samples_mic <- pheno_rows_mic %>% pull(get(pheno_sample_col)) %>% unique()
  pheno_samples_disk <- pheno_rows_disk %>% pull(get(pheno_sample_col)) %>% unique()

  breakpoints <- getBreakpoints(species, guide, antibiotic, type_filter="human")

  summary <- rbind(c(paste("Samples with", paste0(drug_class_list,collapse="/"),"genotypes:"), length(geno_samples)),
                   c(paste("Samples with", antibiotic,"phenotypes:"), length(pheno_samples)),
                   c(" - MIC", length(pheno_samples_mic)),
                   c(" - DISK", length(pheno_samples_disk)),
                   c(paste("Samples with genotypes and phenotypes:"), sum(geno_samples %in% pheno_samples)),
                   c(" - MIC", sum(geno_samples %in% pheno_samples_mic)),
                   c(" - DISK", sum(geno_samples %in% pheno_samples_disk)),
                   c("EUCAST breakpoint sites:", nrow(breakpoints))
  )
  if (nrow(breakpoints)>0) {
    for (i in 1:nrow(breakpoints)) {
      summary <- rbind(summary,
                       c(paste0(" - ", breakpoints$method[i], " / ",breakpoints$site[i]),
                         paste0(breakpoints$breakpoint_S[i], ", ", breakpoints$breakpoint_R[i]))
      )
    }
  }

  ecoff <- getBreakpoints(species, guide, antibiotic, type_filter="ECOFF")
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


#' AMR Genotype-Phenotype Analysis to support defining AMRrules
#'
#' This function runs a full AMR rule analysis pipeline, using the AMR and AMRgen packages, for a specified antibiotic and organism.
#' It includes EUCAST reference distribution plotting, phenotype/genotype overlap, solo PPV calculations, logistic regression,
#' upset plots for disk and MIC data, and summary statistics for AMR markers.
#'
#' It outputs data and plots used to assess AMR marker relevance, which can be used as input to the `makerules` function.
#'
#' @param geno_table A data frame of genotypic data with a `drug_class` column, sample IDs (specified using `geno_sample_col`), and AMR determinants (`marker`). Can be generated from AMRfinderplus output using `AMRgen::import_amrfp`.
#' @param pheno_table A data frame of phenotypic data, including columns for sample IDs (specified using `pheno_sample_col`), antibiotic (`drug_agent`), MIC (`mic`), disk diffusion (`disk`), and S/I/R classification (specified using `sir_col`).
#' @param antibiotic A string specifying the antibiotic to analyze (e.g., `"Ciprofloxacin"`).
#' @param drug_class_list A character vector of drug classes whose attributed markers should be included in the analysis.
#' @param species A string indicating the species name used for breakpoint and reference distribution lookups (e.g., `"E. coli"`).
#' @param sir_col The name of the phenotypic classification column to use (default: `"pheno"`).
#' @param geno_sample_col The column in \code{geno_table} with sample IDs (default: `"Name"`).
#' @param pheno_sample_col The column in \code{pheno_table} with sample IDs (default: `"id"`).
#' @param marker_col A character string specifying the column name in `geno_table` containing the marker identifiers. Defaults to `"marker.label"`.
#' @param minPPV Minimum number of samples required for a marker to be included in PPV analysis (default: 1).
#' @param mafLogReg Minor allele frequency threshold for logistic regression inclusion (default: 5).
#' @param mafUpset Minor allele frequency threshold for upset plot inclusion (default: 5).
#'
#' @return A list containing:
#' \item{reference_mic_plot}{EUCAST reference MIC distribution plot}
#' \item{reference_disk_plot}{EUCAST reference disk zone distribution plot}
#' \item{summary}{Output of \code{\link{summarise_data}} showing sample and breakpoint summaries}
#' \item{solo_stats}{PPV statistics for individual markers}
#' \item{solo_binary}{Binary matrix of individual marker presence by sample}
#' \item{amr_binary}{Binary matrix of AMR marker presence by sample}
#' \item{ppv_plot}{Bar plot summarizing PPV results}
#' \item{logistic_mat}{Binary matrix used in logistic regression}
#' \item{logistic_plot}{Plot of logistic regression estimates}
#' \item{ppv_logistic_plot}{Combined PPV/logistic plot}
#' \item{modelR}{Logistic model for predicting resistance}
#' \item{modelNWT}{Logistic model for predicting NWT (non-wild-type)}
#' \item{allstatsR}{Merged statistics for R category}
#' \item{allstatsNWT}{Merged statistics for NWT category}
#' \item{upset_mic_plot}{Upset plot for MIC data}
#' \item{upset_disk_plot}{Upset plot for disk data}
#' \item{upset_mic_summary}{MIC data summarised per marker or combination}
#' \item{upset_disk_summary}{Disk diffusion data summarised per marker or combination}
#' \item{combination_summary_values}{Summary of genotype combinations from both MIC and disk}
#' \item{afp_hits}{List of AMR markers detected}
#' \item{species}{The species used in the analysis}
#' \item{antibiotic}{The antibiotic used in the analysis}
#'
#' @details
#' This function is designed to assist in AMR rule derivation by integrating multiple evidence streams: EUCAST distributions,
#' solo PPV calculations, logistic regression models, and marker co-occurrence patterns.
#'
#' All steps are wrapped in `safe_execute()` to allow the pipeline to run partially in case of errors.
#'
#' @examples
#' \dontrun{
#' amrrules_analysis(geno_table, pheno_table, antibiotic = "Ciprofloxacin",
#'                   drug_class_list = c("Quinolones"), species = "E. coli")
#' }
#'
#' @export
amrrules_analysis <- function(geno_table, pheno_table, antibiotic, drug_class_list, species, sir_col="pheno",
                              geno_sample_col="Name", pheno_sample_col="id", marker_col="marker.label",
                              minPPV=1, mafLogReg=5, mafUpset=5) {

  # plot EUCAST reference distributions
  reference_mic <- safe_execute(AMRgen::get_eucast_mic_distribution(antibiotic, species))
  reference_mic <- safe_execute(rep(reference_mic$mic, reference_mic$count))
  reference_mic_plot <- safe_execute(ggplot2::autoplot(reference_mic, ab = antibiotic, mo = species, title = "EUCAST reference MIC distribtion"))

  reference_disk <- safe_execute(AMRgen::get_eucast_disk_distribution(antibiotic, species))
  reference_disk <- safe_execute(rep(reference_disk$disk_diffusion, reference_disk$count))
  reference_disk_plot <- safe_execute(ggplot2::autoplot(reference_disk, ab = antibiotic, mo = species, title = "EUCAST reference disk zone distribtion"))

  summary <- safe_execute(summarise_data(geno_table, pheno_table, antibiotic=antibiotic, drug_class_list=drug_class_list, geno_sample_col=geno_sample_col, pheno_sample_col=pheno_sample_col, species=species))

  cat("Running solo PPV analysis\n")
  soloPPV <- safe_execute(AMRgen::solo_ppv_analysis(geno_table=geno_table, pheno_table=pheno_table, antibiotic=antibiotic, drug_class_list=drug_class_list, sir_col=sir_col, min=minPPV, marker_col=marker_col))

  cat("Running logistic regression\n")
  logistic <- safe_execute(AMRgen::amr_logistic(geno_table=geno_table, pheno_table=pheno_table, antibiotic=antibiotic, drug_class_list=drug_class_list, sir_col=sir_col, maf=mafLogReg, geno_sample_col=geno_sample_col, pheno_sample_col=pheno_sample_col, marker_col=marker_col, ecoff_col=NULL))

  cat("Summarising stats\n")
  allstatsR <- safe_execute(AMRgen::merge_logreg_soloppv(logistic$modelR, soloPPV$solo_stats %>% filter(category=="R"),
                                                 title=paste0(paste0(drug_class_list, collapse="/")," markers vs ",antibiotic, " R")))

  allstatsNWT <- safe_execute(AMRgen::merge_logreg_soloppv(logistic$modelNWT, soloPPV$solo_stats %>% filter(category=="NWT"),
                                                   title=paste0(paste0(drug_class_list, collapse="/")," markers vs ",antibiotic, " NWT")))

  # combined PPV/logistic plot
  ppv_logistic_plot <- safe_execute(soloPPV$combined_plot + AMRgen::compare_estimates(logistic$modelNWT, logistic$modelR, title1="NWT", title2="R",
                                                                              colors=c(NWT="blue4", R="maroon"), y_title="",
                                                                              marker_order=levels(as.factor(soloPPV$solo_stats$marker))) +
                                      theme(legend.position="none") + ggtitle("Logistic regression", subtitle="for R and NWT"))

  cat("Generating upset plots\n")

  upset_mic <- safe_execute(AMRgen::amr_upset(soloPPV$amr_binary %>% filter(!is.na(pheno)), min_set_size=mafUpset, order="value", assay="mic"))

  upset_disk <- safe_execute(AMRgen::amr_upset(soloPPV$amr_binary %>% filter(!is.na(pheno)), min_set_size=mafUpset, order="value", assay="disk"))
  
  if (!is.null(upset_mic$summary) & !is.null(upset_disk$summary)) {
    combination_summary_values <- safe_execute(full_join(upset_mic$summary, upset_disk$summary, by=c("marker_list", "marker_count"), suffix=c(".mic", ".disk")))
  } else {
    combination_summary_values <- NULL
  }

  # return AMRFP info for all unique markers in this class with any pheno data
  overlap <- AMRgen::compare_geno_pheno_id(geno_table %>% filter(drug_class %in% drug_class_list),
                                   pheno_table %>% filter(drug_agent==as.ab(antibiotic)),
                                   geno_sample_col = geno_sample_col,
                                   pheno_sample_col = pheno_sample_col, rename_id_cols = T)

  afp_hits <- overlap$geno_matched %>%
    filter(`Element type`=="AMR") %>%
    select(any_of(c("marker", "gene", "mutation", "marker.label", "Gene symbol", "Element subtype", "Hierarchy node", "HMM id"))) %>%
    distinct()

  # add gene frequencies to help define core/accessory
  afp_hits <- overlap$geno_matched %>%
    group_by(!!sym(marker_col)) %>%
    count() %>%
    mutate(freq=n/length(unique(overlap$geno_matched$id))) %>%
    select(-n) %>% right_join(afp_hits, by=marker_col)

  # merge geno/pheno and count unique sources per gene + pheno/mic/disk
  soloPPV$amr_binary <- pheno_table %>% select(any_of(c("id", "source", "pheno", "mic", "disk"))) %>%
    left_join(soloPPV$amr_binary, by=c("id", "pheno", "mic", "disk"))

  soloPPV$solo_binary <- pheno_table %>% select(any_of(c("id", "source", "pheno", "mic", "disk"))) %>%
    left_join(soloPPV$solo_binary, by=c("id", "pheno", "mic", "disk"))

  return(list(reference_mic_plot=reference_mic_plot,
              reference_disk_plot=reference_disk_plot,
              summary=summary,
              solo_stats=soloPPV$solo_stats,
              solo_binary=soloPPV$solo_binary,
              amr_binary=soloPPV$amr_binary,
              ppv_plot=soloPPV$combined_plot,
              logistic_mat=logistic$bin_mat,
              logistic_plot=logistic$plot,
              ppv_logistic_plot=ppv_logistic_plot,
              modelR=logistic$modelR,
              modelNWT=logistic$modelNWT,
              allstatsR=allstatsR$combined,
              allstatsNWT=allstatsNWT$combined,
              upset_mic_plot=upset_mic$plot,
              upset_mic_summary=upset_mic$summary,
              upset_disk_plot=upset_disk$plot,
              upset_disk_summary=upset_disk$summary,
              combination_summary_values=combination_summary_values,
              afp_hits=afp_hits,
              species=species,
              antibiotic=antibiotic))
}

#' Save AMR Geno-Pheno Analysis Results and Plots to Files and Draft AMRrules
#'
#' This function saves the results generated using `amrrules_analysis`, including figures and summary tables, to a specified directory.
#' It also uses these results to draft AMRrules using the `makerules` function and saves these to file
#' The output includes plots of reference distributions, PPV analysis, logistic regression, upset plots, and more, along with data summaries and rule tables.
#'
#' @param amrrules A list containing the output of the \code{\link{amrrules_analysis}} function, including plots and statistical summaries.
#' @param width Numeric, the width of the output plot images (default: 9).
#' @param height Numeric, the height of the output plot images (default: 9).
#' @param dir_path The path to the directory where the output will be saved.
#' @param outdir_name The name of the subdirectory within \code{dir_path} to save the results in (default: the antibiotic name).
#' @param file_prefix The prefix for filenames of output files (default: the antibiotic name).
#' @param minObs The minimum number of observations required to include a rule (default: 3).
#' @param weak_threshold The threshold for identifying weak rules (default: 20).
#' @param bp_site The breakpoint site to filter on (Optional, use the \code{\link{summarise_data}} function to check whether there are different breakpoints for different sites and choose which one to specify here. By default, if multiple breakpoints are available the most conservative will be used.).
#' @param ruleID_start The starting numeric ID with which to number rules (integer, default: 1000).
#' @param mic_S The MIC breakpoint to define S (Optional, by default breakpoints are extracted from EUCAST using the AMR package, however if none are available or user wants to use a different one e.g. an ECOFF then it can be specified here).
#' @param mic_R The MIC breakpoint to define R (Optional, by default breakpoints are extracted from EUCAST using the AMR package, however if none are available or user wants to use a different one e.g. an ECOFF then it can be specified here).
#'
#' @return A data frame containing the generated AMR rules.
#'
#' @details
#' This function saves the plots and data tables generated by `amrrules_analysis`, including:
#' - Reference MIC and disk distribution plots (EUCAST)
#' - PPV analysis plots
#' - Logistic regression plots
#' - Summary statistics tables
#' It then uses these data to draft AMRrules, returns the data and rules, and writes them to a TSV file for review.
#'
#' The output is saved in a directory structure specified by \code{dir_path} and \code{outdir_name}.
#' If no output directory name is provided, the default is the antibiotic name.
#'
#' @examples
#' \dontrun{
#' amrrules_save(amrrules, dir_path = "output/", outdir_name = "Ciprofloxacin", file_prefix = "Ciprofloxacin")
#' }
#'
#' @export
amrrules_save <- function(amrrules, width=9, height=9, dir_path, outdir_name=NULL, file_prefix=NULL, minObs=3, weak_threshold=20, bp_site=NULL, ruleID_start=1000, mic_S=NULL, mic_R=NULL, makeRules=TRUE) {

  if (is.null(outdir_name)) { outdir_name <- amrrules$antibiotic }
  if (is.null(file_prefix)) { file_prefix <- amrrules$antibiotic }

  outdir_path <- file.path(dir_path, outdir_name)
  if (!dir.exists(outdir_path)) {
    safe_execute(dir.create(outdir_path, recursive=TRUE))
    cat(paste0("Directory '", outdir_path, "' created successfully.\n"))
  }

  outpath <- file.path(outdir_path, file_prefix)

  cat(paste0("\nWriting figs and tables to ", outpath,"_*\n"))

  if (!is.null(amrrules$reference_mic_plot)) {safe_execute(ggsave(amrrules$reference_mic_plot, filename=paste0(outpath,"_reference_mic_plot.pdf"), width=width, height=height))}
  if (!is.null(amrrules$reference_disk_plot)) {safe_execute(ggsave(amrrules$reference_disk_plot, filename=paste0(outpath,"_reference_disk_plot.pdf"), width=width, height=height))}
  if (!is.null(amrrules$ppv_plot)) {safe_execute(ggsave(amrrules$ppv_plot, filename=paste0(outpath,"_soloPPV_plot.pdf"), width=width, height=height))}
  if (!is.null(amrrules$logistic_plot)) {safe_execute(ggsave(amrrules$logistic_plot, filename=paste0(outpath,"_logistic_plot.pdf"), width=width, height=height))}
  if (!is.null(amrrules$ppv_logistic_plot)) {safe_execute(ggsave(amrrules$ppv_logistic_plot, filename=paste0(outpath,"_soloPPV_logistic_plot.pdf"), width=width, height=height))}
  if (!is.null(amrrules$upset_mic_plot)) {safe_execute(ggsave(amrrules$upset_mic_plot, filename=paste0(outpath,"_MIC_upset_plot.pdf"), width=width, height=height))}
  if (!is.null(amrrules$upset_disk_plot)) {safe_execute(ggsave(amrrules$upset_disk_plot, filename=paste0(outpath,"_disk_upset_plot.pdf"), width=width, height=height))}

  safe_execute(readr::write_tsv(as.data.frame(amrrules$summary), col_names=F, file=paste0(outpath,"_data_summary.tsv")))
  safe_execute(readr::write_tsv(amrrules$allstatsR, file=paste0(outpath,"_stats_R.tsv")))
  safe_execute(readr::write_tsv(amrrules$allstatsNWT, file=paste0(outpath,"_stats_NWT.tsv")))
  
  if (!is.null(amrrules$upset_mic_summary)) {safe_execute(readr::write_tsv(amrrules$upset_mic_summary, file=paste0(outpath,"_MIC_summary.tsv")))}
  else {cat ("  (No MIC data summary available to write)\n")}
  if (!is.null(amrrules$upset_disk_summary)) {safe_execute(readr::write_tsv(amrrules$upset_disk_summary, file=paste0(outpath,"_DD_summary.tsv")))}
  else {cat ("  (No disk diffusion data summary available to write)\n")}

  # make rules and write them out to same directory
  if (makeRules) {
    cat ("\n")
    rules <- safe_execute(makerules(amrrules, minObs=minObs, weak_threshold=weak_threshold, bp_site=bp_site, ruleID_start=ruleID_start, mic_S=mic_S, mic_R=mic_R))
    if (!is.null(rules)) {
      safe_execute(readr::write_tsv(rules$rules, file=paste0(outpath,"_AMRrules.tsv")))
      safe_execute(readr::write_tsv(rules$data, file=paste0(outpath,"_AMRrules_data.tsv")))
      cat(paste0("\nWriting rules to ",outpath,"_AMRrules.tsv\n"))
      }
    else{cat("Failed to make rules\n\n")}
  }
  else {rules <- NULL}

  cat ("\n")
  
  return(rules)
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
#' @param bp_standard A character string to return as the breakpoint standard (optional). This is updated with the selected breakpoint site.
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
#' checkBreakpoints(species="Escherichia coli", guide="EUCAST 2024"", antibiotic="Ciprofloxacin", assay="MIC")
#'
#' @export
checkBreakpoints <- function(species, guide="EUCAST 2024", antibiotic, assay="MIC", bp_site=NULL, bp_standard="") {
  breakpoints <- getBreakpoints(species, guide, antibiotic) %>% filter(method==assay)
  if (nrow(breakpoints)==0) {stop(paste("Could not determine",assay,"breakpoints using AMR package, please provide your own breakpoints"))}
  else{
    breakpoint_sites <- unique(breakpoints$site)
    breakpoint_message_multibp <- NA
    # handle multiple breakpoints (e.g. for different conditions)
    if (length(breakpoint_sites)>1) {
      breakpoint_message_multibp <- paste("NOTE: Multiple breakpoint entries, for different sites:", paste(breakpoint_sites, collapse="; "))
      if (length(unique(breakpoints$breakpoint_R))==1 & length(unique(breakpoints$breakpoint_S))==1) {
        breakpoints <- breakpoints %>% arrange(-breakpoint_S) %>% first()
        breakpoint_message_multibp <- paste0(breakpoint_message_multibp, ". However S and R breakpoints are the same.")
      }
      else if (is.null(bp_site)) {
        breakpoints <- breakpoints %>% arrange(-breakpoint_S) %>% first()
        breakpoint_message_multibp <- paste0(breakpoint_message_multibp, ". Using the one with the highest S breakpoint (", breakpoints$site,").")
        bp_standard <- paste0(" (", breakpoints$site, ")")
      }
      else {
        if (bp_site %in% breakpoint_sites) {
          breakpoints <- breakpoints %>% filter(site==bp_site)
          breakpoint_message_multibp <- paste0(breakpoint_message_multibp, ". Using the specified site (", bp_site,").")
          bp_standard <- paste0(" (", bp_site, ")")
        }
        else {
          breakpoints <- breakpoints %>% arrange(-breakpoint_S) %>% first()
          breakpoint_message_multibp <- paste0(breakpoint_message_multibp, ". Could not find the specified site (", bp_site,"), so using the one with the highest S breakpoint (", breakpoints$site,").")
          bp_standard <- paste0(" (", breakpoints$site, ")")
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

getSources <- function(amr_binary, combo, assay) {

  # assume combo list is in proper syntax (with ':') but amr_binary has colnames with '..' in place of ':'
  marker_names <- unlist(str_split(gsub(":", "..", combo), ", "))

  total_sources <- amr_binary %>%
    filter(!is.na(get(assay))) %>%
    select(source, all_of(marker_names)) %>%
    filter(if_all(all_of(marker_names), ~ . == 1)) %>%
    distinct() %>% nrow()

  sources_per_pheno <- amr_binary %>%
    filter(!is.na(get(assay))) %>%
    select(source, pheno, all_of(marker_names)) %>%
    filter(if_all(all_of(marker_names), ~ . == 1)) %>%
    distinct() %>%
    group_by(pheno) %>% count(name="sources") %>%
    pivot_wider(names_from=pheno, names_prefix="sources.", values_from=sources) %>%
    mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))

  sources <- as_tibble(c(sources=total_sources, sources_per_pheno))

  source_names <- c("sources", "sources.S", "sources.I", "sources.R")
  sources <- add_missing_cols(sources, source_names) %>%
    rename_with(.fn = ~ paste0(assay,".", .x))

  return(sources)
}


convert_mutation <- function(mut) {
  aa_mapping <- c(
    A = "Ala", R = "Arg", N = "Asn", D = "Asp", C = "Cys", E = "Glu",
    Q = "Gln", G = "Gly", H = "His", I = "Ile", L = "Leu", K = "Lys",
    M = "Met", F = "Phe", P = "Pro", S = "Ser", T = "Thr", W = "Trp",
    Y = "Tyr", V = "Val"
  )
  aa1 <- substr(mut, 1, 1)
  position <- substr(mut, 2, nchar(mut) - 1)
  aa2 <- substr(mut, nchar(mut), nchar(mut))
  return(paste0("p.", aa_mapping[aa1], position, aa_mapping[aa2]))
}
