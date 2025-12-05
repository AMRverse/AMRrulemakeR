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
#' amrrules_analysis2(geno_table, pheno_table, antibiotic = "Ciprofloxacin",
#'                   drug_class_list = c("Quinolones"), species = "E. coli")
#' }
#'
#' @export
amrrules_analysis2 <- function(geno_table, pheno_table, antibiotic, drug_class_list, species, 
                               sir_col="pheno_eucast", ecoff_col="ecoff", sir_provided_col="pheno_provided",
                              geno_sample_col="Name", pheno_sample_col="id", marker_col="marker.label",
                              minPPV=1, mafLogReg=5, mafUpset=5, info=NULL) {
  
  # plot EUCAST reference distributions
  reference_mic <- safe_execute(AMRgen::get_eucast_mic_distribution(antibiotic, species))
  reference_mic <- safe_execute(rep(reference_mic$mic, reference_mic$count))
  reference_mic_plot <- safe_execute(ggplot2::autoplot(reference_mic, ab = antibiotic, mo = species, title = "EUCAST reference MIC distribtion"))
  
  reference_disk <- safe_execute(AMRgen::get_eucast_disk_distribution(antibiotic, species))
  reference_disk <- safe_execute(rep(reference_disk$disk_diffusion, reference_disk$count))
  reference_disk_plot <- safe_execute(ggplot2::autoplot(reference_disk, ab = antibiotic, mo = species, title = "EUCAST reference disk zone distribtion"))
  
  ### TO DO: plot MIC and disk distributions for input data, stratified by method
  
  ### TO DO: improve this
  #summary <- safe_execute(summarise_data(geno_table, pheno_table, antibiotic=antibiotic, drug_class_list=drug_class_list, geno_sample_col=geno_sample_col, pheno_sample_col=pheno_sample_col, species=species))
  
  
  ### NEW: filter to pheno data with disk/MIC measures first
  
  cat("Running solo PPV analysis on samples with MIC or disk measures\n")
  pheno_table_micdisk <- pheno_table %>% filter(!is.na(mic) | !is.na(disk))
  soloPPV_micdisk <- safe_execute(AMRgen::solo_ppv_analysis(geno_table=geno_table, pheno_table=pheno_table_micdisk, antibiotic=antibiotic, drug_class_list=drug_class_list, sir_col=sir_col, ecoff_col=ecoff_col, min=minPPV, marker_col=marker_col))
  
  cat("Running logistic regression on samples with MIC or disk measures\n")
  logistic_micdisk <- safe_execute(AMRgen::amr_logistic(geno_table=geno_table, pheno_table=pheno_table_micdisk, antibiotic=antibiotic, drug_class_list=drug_class_list, sir_col=sir_col, ecoff_col=ecoff_col, maf=mafLogReg, geno_sample_col=geno_sample_col, pheno_sample_col=pheno_sample_col, marker_col=marker_col))
  
  cat("Combining PPV and regression stats for samples with MIC or disk measures\n")
  # combined PPV/logistic plot for samples with MIC or disk measures
  ppv_logistic_plot_micdisk <- safe_execute(soloPPV_micdisk$combined_plot + logistic_micdisk$plot + scale_y_discrete(limits = names(soloPPV_micdisk$plot_order)) + labs(y="") +
                                      theme(legend.position="none") + ggtitle("Logistic regression"))
  
  cat("Generating upset plots, with PPV and median measures for marker combinations\n")
  
  cat("...MIC ")
  upset_mic <- safe_execute(AMRgen::amr_upset(soloPPV_micdisk$amr_binary, min_set_size=mafUpset, order="value", assay="mic"))
  
  cat("...disk ")
  upset_disk <- safe_execute(AMRgen::amr_upset(soloPPV_micdisk$amr_binary, min_set_size=mafUpset, order="value", assay="disk"))
  
  cat("\n")
  
  
  ### NEW: merge pheno and ecoff data from MIC/disk values with those from SIR, to give tier-2 info
  
  pheno_table_sir <- pheno_table %>% 
    filter(is.na(mic) & is.na(disk))
  
  if (nrow(pheno_table_sir)>0) {
   pheno_table_sir <- pheno_table_sir %>% 
    rename(pheno=!!sir_provided_col) %>%
    mutate(ecoff=case_when(!is.na(!!ecoff_col) ~ get(ecoff_col), # define ecoff from SIR if not already done
                           pheno %in% c("I", "R") ~ as.sir("R"),
                           pheno %in% c("S") ~ as.sir("R"),
                           TRUE ~ NA)) %>%
    bind_rows(pheno_table_micdisk %>% rename(pheno=!!sir_col) %>% rename(ecoff=!!ecoff_col))
  
    cat("Running solo PPV analysis on all samples, including those with SIR calls only\n")
    soloPPV_sir <- safe_execute(AMRgen::solo_ppv_analysis(geno_table=geno_table, pheno_table=pheno_table_sir, antibiotic=antibiotic, drug_class_list=drug_class_list, sir_col="pheno", ecoff_col="ecoff", min=minPPV, marker_col=marker_col))
  
    cat("Running logistic regression on all samples, including those with SIR calls only\n")
    logistic_sir <- safe_execute(AMRgen::amr_logistic(geno_table=geno_table, pheno_table=pheno_table_sir, antibiotic=antibiotic, drug_class_list=drug_class_list, sir_col="pheno", ecoff_col="ecoff", maf=mafLogReg, geno_sample_col=geno_sample_col, pheno_sample_col=pheno_sample_col, marker_col=marker_col))
  
    cat("Combining PPV and regression stats for all samples\n")
    # combined PPV/logistic plot for samples with MIC or disk measures
    ppv_logistic_plot_sir <- safe_execute(soloPPV_sir$combined_plot + logistic_sir$plot + scale_y_discrete(limits = names(soloPPV_sir$plot_order)) + labs(y="") +
                                              theme(legend.position="none") + ggtitle("Logistic regression"))
  }
  
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
  
  # extract info for relevant samples (for source enumeration)
  if (!is.null(info)) {
    samples_with_pheno <- pheno_table %>% filter(drug_agent==as.ab(antibiotic))
    colnames(info)[1] <- "id"
    info <- info %>% filter(id %in% samples_with_pheno$id)
  }
  
  return(list(reference_mic_plot=reference_mic_plot,
              reference_disk_plot=reference_disk_plot,
              solo_stats=soloPPV_micdisk$solo_stats,
              solo_binary=soloPPV_micdisk$solo_binary,
              amr_binary=soloPPV_micdisk$amr_binary,
              ppv_plot=soloPPV_micdisk$combined_plot,
              logistic_mat=logistic_micdisk$bin_mat,
              logistic_plot=logistic_micdisk$plot,
              ppv_logistic_plot=ppv_logistic_plot_micdisk,
              modelR=logistic_micdisk$modelR,
              modelNWT=logistic_micdisk$modelNWT,
              upset_mic_plot=upset_mic$plot,
              upset_mic_summary=upset_mic$summary,
              upset_disk_plot=upset_disk$plot,
              upset_disk_summary=upset_disk$summary,
              solo_stats_all=soloPPV_sir$solo_stats,
              solo_binary_all=soloPPV_sir$solo_binary,
              amr_binary_all=soloPPV_sir$amr_binary,
              ppv_plot_all=soloPPV_sir$combined_plot,
              logistic_mat_all=logistic_sir$bin_mat,
              logistic_plot_all=logistic_sir$plot,
              ppv_logistic_plot_all=ppv_logistic_plot_sir,
              modelR_all=logistic_sir$modelR,
              modelNWT_all=logistic_sir$modelNWT,
              afp_hits=afp_hits,
              species=species,
              antibiotic=antibiotic,
              info=info))
}

#' Save AMR Geno-Pheno Analysis Results and Plots to Files and Draft AMRrules
#'
#' This function saves the results generated using `amrrules_analysis2`, including figures and summary tables, to a specified directory.
#' It also uses these results to draft AMRrules using the `makerules2` function and saves these to file
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
#' amrrules_save2(amrrules, dir_path = "output/", outdir_name = "Ciprofloxacin", file_prefix = "Ciprofloxacin")
#' }
#'
#' @export
amrrules_save2 <- function(amrrules, width=9, height=9, dir_path, outdir_name=NULL, file_prefix=NULL, 
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
  
  #safe_execute(readr::write_tsv(as.data.frame(amrrules$summary), col_names=F, file=paste0(outpath,"_data_summary.tsv")))
  safe_execute(readr::write_tsv(amrrules$modelR, file=paste0(outpath,"_model_R.tsv")))
  safe_execute(readr::write_tsv(amrrules$modelNWT, file=paste0(outpath,"_model_NWT.tsv")))
  safe_execute(readr::write_tsv(amrrules$solo_stats, file=paste0(outpath,"_soloPPV.tsv")))
  
  if (!is.null(amrrules$upset_mic_summary)) {safe_execute(readr::write_tsv(amrrules$upset_mic_summary, file=paste0(outpath,"_MIC_summary.tsv")))}
  else {cat ("  (No MIC data summary available to write)\n")}
  if (!is.null(amrrules$upset_disk_summary)) {safe_execute(readr::write_tsv(amrrules$upset_disk_summary, file=paste0(outpath,"_DD_summary.tsv")))}
  else {cat ("  (No disk diffusion data summary available to write)\n")}
  
  safe_execute(readr::write_tsv(amrrules$modelR_all, file=paste0(outpath,"_model_R_all.tsv")))
  safe_execute(readr::write_tsv(amrrules$modelNWT_all, file=paste0(outpath,"_model_NWT_all.tsv")))
  safe_execute(readr::write_tsv(amrrules$solo_stats_all, file=paste0(outpath,"_soloPPV_all.tsv")))
  
  if (!is.null(amrrules$ppv_plot_all)) {safe_execute(ggsave(amrrules$ppv_plot_all, filename=paste0(outpath,"_soloPPV_plot_all.pdf"), width=width, height=height))}
  if (!is.null(amrrules$logistic_plot_all)) {safe_execute(ggsave(amrrules$logistic_plot_all, filename=paste0(outpath,"_logistic_plot_all.pdf"), width=width, height=height))}
  if (!is.null(amrrules$ppv_logistic_plot_all)) {safe_execute(ggsave(amrrules$ppv_logistic_plot_all, filename=paste0(outpath,"_soloPPV_logistic_plot_all.pdf"), width=width*1.5, height=height))}
  
  # make rules and write them out to same directory
  if (makeRules) {
    cat ("\n")
    rules <- safe_execute(makerules2(amrrules, minObs=minObs, low_threshold=low_threshold, bp_site=bp_site, 
                                    ruleID_start=ruleID_start, mic_S=mic_S, mic_R=mic_R,
                                    disk_S=disk_S, disk_R=disk_R, use_disk=use_disk, guide=guide))
    if (!is.null(rules)) {
      safe_execute(readr::write_tsv(rules$rules, file=paste0(outpath,"_AMRrules.tsv")))
      safe_execute(readr::write_tsv(rules$data, file=paste0(outpath,"_AMRrules_data.tsv")))
      cat(paste0("\nWriting rules to ",outpath,"_AMRrules.tsv\n"))
      
      cat ("  Comparing primary calls (solo PPV from MIC/disk data) with other evidence\n")
      rule_compare_R <- safe_execute(compareRulesData(rules$data, antibiotic=amrrules$antibiotic, drug_class_list=amrrules$drug_class_list, type="R"))
      rule_compare_NWT <- safe_execute(compareRulesData(rules$data, antibiotic=amrrules$antibiotic, drug_class_list=amrrules$drug_class_list, type="NWT"))  
      
      if (!is.null(rule_compare_R$plot)) {safe_execute(ggsave(rule_compare_R$plot, filename=paste0(outpath,"_soloPPV_R_vsOther.pdf"), width=width, height=height/2))}
      if (!is.null(rule_compare_NWT$plot)) {safe_execute(ggsave(rule_compare_NWT$plot, filename=paste0(outpath,"_soloPPV_NWT_vsOther.pdf"), width=width, height=height/2))}
      safe_execute(readr::write_tsv(rule_compare_R$diffs, file=paste0(outpath,"_soloPPV_R_mismatchCalls.tsv")))
      safe_execute(readr::write_tsv(rule_compare_NWT$diffs, file=paste0(outpath,"_soloPPV_NWT_mismatchCalls.tsv")))
    }
    else{cat("Failed to make rules\n\n")}
  }
  else {
    cat("Not making rules (rerun with makeRules=TRUE to generate rules)\n")
    rules <- NULL}
  
  cat ("\n")
  
  return(rules)
}

compareRulesData <- function(rules_dat, antibiotic, drug_class_list, type="R") {
  
  if (type=="R") {exclude<-"NWT"}
  else if (type=="NWT") {exclude<-"R"}
  
  # compare using rules data
  compare_r_ppv <- rules_dat %>%
    select(marker, ends_with("ppv")) %>% 
    select(marker, starts_with(type)) %>% 
    pivot_longer(cols=-marker)
  
  solo_ppv_colname=paste0(type,".solo.ppv")
  solo_ppv_sym <- sym(solo_ppv_colname)
  
  compare_r_ppv <- compare_r_ppv %>% filter(name==solo_ppv_colname) %>% 
    select(marker, value) %>% rename(!!solo_ppv_sym:=value) %>% 
    left_join(compare_r_ppv %>% filter(name!=solo_ppv_colname), by="marker") %>%
    mutate(match=case_when(!!solo_ppv_sym>0.5 & value>0.5 ~ TRUE,
                           !!solo_ppv_sym<=0.5 & value<=0.5 ~ TRUE,
                           is.na(!!solo_ppv_sym) | is.na(value) ~ NA,
                           TRUE ~ FALSE)) 
  
  mic_colname=paste0(type,".MIC.n")
  disk_colname=paste0(type,".Disk.n")
  
  compare_r_n <- rules_dat %>%
    select(marker, ends_with(".n")) %>% 
    select(-starts_with(exclude)) 
    
  if (mic_colname %in% colnames(compare_r_n)){
    compare_r_n <- compare_r_n %>% rename(!!mic_colname:=MIC.n)
  }
  if (disk_colname %in% colnames(compare_r_n)) {
    compare_r_n <- compare_r_n %>% rename(!!disk_colname:=Disk.n)
  }
  
  compare_r_n <- compare_r_n %>% 
    pivot_longer(cols=-marker)
  
  plot <- ggplot(data=compare_r_ppv, aes(x=.data[[solo_ppv_colname]], 
                                         y=value, col=match)) + 
    geom_abline(slope=1,intercept=0) + 
    geom_point() + 
    facet_wrap(~name) +
    labs(x=paste0("Solo PPV for ",type,", using data with MIC/disk assay measures"),
         y="Alternative solo PPV estimates",
         title=paste("Solo markers for class:", 
                     paste0(drug_class_list, collapse = ", ")),
         subtitle=paste("vs phenotype for drug:", ab_name(as.ab(antibiotic)))
    )
  
  solo_n_colname=paste0(type,".solo.n")
  
  compare_r_n <- compare_r_n %>% filter(name==solo_n_colname) %>% 
    select(marker, value) %>% rename(!!solo_n_colname:=value) %>% 
    left_join(compare_r_n %>% filter(name!=solo_n_colname), by="marker") %>%
    mutate(name=gsub(".n", ".ppv", name)) %>% 
    rename(alt.n=value)
  
  diffs <- compare_r_ppv %>% filter(!match) %>% left_join(compare_r_n, by = join_by(marker, name))
  
  plot 
  
  return(list(diffs=diffs, plot=plot))
}