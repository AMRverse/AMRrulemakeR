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
#' @param info Data frame indicating information about each sample, so we can summarise the source/s of data supporting each rule. First column should be sample id, matching that in the `pheno_table` object.
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
amrrules_analysis <- function(geno_table, pheno_table, antibiotic, drug_class_list, species, sir_col="pheno", ecoff_col="ecoff",
                              geno_sample_col="Name", pheno_sample_col="id", marker_col="marker.label",
                              minPPV=1, mafLogReg=5, mafUpset=5, info=NULL) {
  
  # plot EUCAST reference distributions
  reference_mic <- safe_execute(AMRgen::get_eucast_mic_distribution(antibiotic, species))
  reference_mic <- safe_execute(rep(reference_mic$mic, reference_mic$count))
  reference_mic_plot <- safe_execute(ggplot2::autoplot(reference_mic, ab = antibiotic, mo = species, title = "EUCAST reference MIC distribtion"))
  
  reference_disk <- safe_execute(AMRgen::get_eucast_disk_distribution(antibiotic, species))
  reference_disk <- safe_execute(rep(reference_disk$disk_diffusion, reference_disk$count))
  reference_disk_plot <- safe_execute(ggplot2::autoplot(reference_disk, ab = antibiotic, mo = species, title = "EUCAST reference disk zone distribtion"))
  
  summary <- safe_execute(summarise_data(geno_table, pheno_table, antibiotic=antibiotic, drug_class_list=drug_class_list, geno_sample_col=geno_sample_col, pheno_sample_col=pheno_sample_col, species=species))
  
  cat("Running solo PPV analysis\n")
  soloPPV <- safe_execute(AMRgen::solo_ppv_analysis(geno_table=geno_table, pheno_table=pheno_table, antibiotic=antibiotic, drug_class_list=drug_class_list, sir_col=sir_col, ecoff_col=ecoff_col, min=minPPV, marker_col=marker_col))
  
  cat("Running logistic regression\n")
  logistic <- safe_execute(AMRgen::amr_logistic(geno_table=geno_table, pheno_table=pheno_table, antibiotic=antibiotic, drug_class_list=drug_class_list, sir_col=sir_col, ecoff_col=ecoff_col, maf=mafLogReg, geno_sample_col=geno_sample_col, pheno_sample_col=pheno_sample_col, marker_col=marker_col))
  
  cat("Summarising stats\n")
  if (!is.null(logistic$modelR)) {
    allstatsR <- safe_execute(AMRgen::merge_logreg_soloppv(logistic$modelR, soloPPV$solo_stats %>% filter(category=="R"),
                                                           title=paste0(paste0(drug_class_list, collapse="/")," markers vs ",antibiotic, " R")))
  } else {allstatsR <- NULL}
  
  if (!is.null(logistic$modelNWT)) {
    allstatsNWT <- safe_execute(AMRgen::merge_logreg_soloppv(logistic$modelNWT, soloPPV$solo_stats %>% filter(category=="NWT"),
                                                             title=paste0(paste0(drug_class_list, collapse="/")," markers vs ",antibiotic, " NWT")))
  } else {allstatsNWT <- NULL}
  
  # combined PPV/logistic plot
  ppv_logistic_plot <- safe_execute(soloPPV$combined_plot + logistic$plot + scale_y_discrete(limits = names(soloPPV$plot_order)) + labs(y="") +
                                      theme(legend.position="none") + ggtitle("Logistic regression"))
  
  cat("Generating upset plots\n")
  
  cat("...MIC ")
  upset_mic <- safe_execute(AMRgen::amr_upset(soloPPV$amr_binary, min_set_size=mafUpset, order="value", assay="mic"))
  
  cat("...disk ")
  upset_disk <- safe_execute(AMRgen::amr_upset(soloPPV$amr_binary, min_set_size=mafUpset, order="value", assay="disk"))
  
  cat("\n")
  
  if (!is.null(upset_mic$summary) & !is.null(upset_disk$summary)) {
    combination_summary_values <- safe_execute(full_join(upset_mic$summary, upset_disk$summary, by=c("marker_list", "marker_count"), suffix=c(".mic", ".disk")))
  } else { combination_summary_values <- NULL }
  
  cat("Extracting gene info\n")
  
  # return AMRFP info for all unique markers in this class with any pheno data
  overlap <- AMRgen::compare_geno_pheno_id(geno_table %>% filter(drug_class %in% drug_class_list),
                                           pheno_table %>% filter(drug_agent==as.ab(antibiotic)),
                                           geno_sample_col = geno_sample_col,
                                           pheno_sample_col = pheno_sample_col, rename_id_cols = T)
  
  afp_hits <- overlap$geno_matched %>%
    filter(`Element type`=="AMR") %>%
    select(any_of(c("marker", "gene", "mutation", "node", "variation type", "marker.label", "Gene symbol", "Element subtype", "HMM id"))) %>%
    distinct()
  
  # add gene frequencies to help define core/accessory
  afp_hits <- overlap$geno_matched %>%
    group_by(!!sym(marker_col)) %>%
    summarise(freq_n=n()) %>%
    mutate(freq=freq_n/length(unique(overlap$geno_matched$id))) %>%
    right_join(afp_hits, by=marker_col)
  
  # extract info for relevant samples
  if (!is.null(info)) {
    samples_with_pheno <- pheno_table %>% filter(drug_agent==as.ab(antibiotic))
    colnames(info)[1] <- "id"
    info <- info %>% filter(id %in% samples_with_pheno$id)
  }
  
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
              antibiotic=antibiotic,
              info=info))
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
#' @param low_threshold The threshold below which rules are assigned to evidence grade 'low' (default: 20).
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
amrrules_save <- function(amrrules, width=9, height=9, dir_path, outdir_name=NULL, file_prefix=NULL, 
                          minObs=3, low_threshold=20, bp_site=NULL, ruleID_start=1000, 
                          mic_S=NULL, mic_R=NULL, disk_S=NULL, disk_R=NULL, 
                          use_mic=TRUE, use_disk=TRUE, makeRules=TRUE, guide="EUCAST 2025") {
  
  if (is.null(outdir_name)) { outdir_name <- gsub("/","_",amrrules$antibiotic) }
  if (is.null(file_prefix)) { file_prefix <- gsub("/","_",amrrules$antibiotic) }
  
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
  if (!is.null(amrrules$ppv_logistic_plot)) {safe_execute(ggsave(amrrules$ppv_logistic_plot, filename=paste0(outpath,"_soloPPV_logistic_plot.pdf"), width=width*1.5, height=height))}
  if (!is.null(amrrules$upset_mic_plot)) {safe_execute(ggsave(amrrules$upset_mic_plot, filename=paste0(outpath,"_MIC_upset_plot.pdf"), width=width*1.5, height=height))}
  if (!is.null(amrrules$upset_disk_plot)) {safe_execute(ggsave(amrrules$upset_disk_plot, filename=paste0(outpath,"_disk_upset_plot.pdf"), width=width*1.5, height=height))}
  
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
    rules <- safe_execute(makerules(amrrules, minObs=minObs, low_threshold=low_threshold, bp_site=bp_site, 
                                    ruleID_start=ruleID_start, mic_S=mic_S, mic_R=mic_R,
                                    disk_S=disk_S, disk_R=disk_R, use_disk=use_disk, guide=guide))
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
