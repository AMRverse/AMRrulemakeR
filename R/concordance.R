# ===================================================================== #
#  Licensed as GPL-v3.0.                                                #
#                                                                       #
#  Developed as part of the AMRverse (https://github.com/AMRverse):     #
#  https://github.com/AMRverse/AMRgen                                   #
#                                                                       #
#  We created this package for both routine data analysis and academic  #
#  research and it was publicly released in the hope that it will be    #
#  useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
#                                                                       #
#  This R package is free software; you can freely use and distribute   #
#  it for both personal and commercial purposes under the terms of the  #
#  GNU General Public License version 3.0 (GNU GPL-3), as published by  #
#  the Free Software Foundation.                                        #
# ===================================================================== #

#' Calculate genotype-phenotype concordance from rules-based predictions
#'
#' Predictions are
#' compared to the observed phenotypes using standard classification metrics
#' (via the `yardstick` pkg) and AMR-specific error rates (major error, ME
#' and very major error, VME) per ISO 20776-2 (and see
#' [FDA definitions](https://www.fda.gov/medical-devices/guidance-documents-medical-devices-and-radiation-emitting-products/antimicrobial-susceptibility-test-ast-systems-class-ii-special-controls-guidance-industry-and-fda)).
#'
#' @param data A data frame containing observed and predicted phenotypes.
#' @param pheno_col String indicating the name of the column containing the true phenotypes (default `"pheno_eucast"`).
#' @param pred_col String indicating the name of the column containing the predicted phenotypes (default `"category"`).
#'
#' @details
#' ...
#'
#' Standard metrics (sensitivity, specificity, PPV, NPV, accuracy, kappa,
#' F-measure) are calculated using pkg `yardstick`. AMR-specific error rates are
#' computed internally:
#' - **VME** (Very Major Error): FN / (TP + FN) = 1 - sensitivity. Proportion of
#'   truly resistant isolates not predicted as such from genotype.
#' - **ME** (Major Error): FP / (TN + FP) = 1 - specificity. Proportion of
#'   truly susceptible isolates incorrectly predicted resistant from genotype.
#'
#' @return A list containing:
#' - `metrics`: A tibble with columns `metric`, `estimate`.
#' - `matrix`: Confusion matrix comparing observed vs predicted.
#'
#' @importFrom dplyr bind_rows mutate select
#' @importFrom tibble tibble
#' @importFrom yardstick conf_mat sens spec ppv npv accuracy kap f_meas
#' @seealso [yardstick]
#' @export
#' @examples
#' \dontrun{
#' metrics_SIR <- rules_concordance(true_vs_predict))
#' metrics_NWT <- rules_concordance(true_vs_predict, "ecoff", "phenotype"))
#' }
rules_concordance <- function(df, pheno_col="pheno_eucast", pred_col="category") {

  df <- df %>%
    mutate(truth_value=case_when(get(pheno_col) %in% c("R", "NWT", "nonwildtype") ~ 1,
                                 get(pheno_col) %in% c("S", "I", "WT", "wildtype") ~ 0,
                           TRUE ~ NA)
           ) %>%
    mutate(est_value=case_when(get(pred_col) %in% c("R", "NWT", "nonwildtype") ~ 1,
                               get(pred_col) %in% c("S", "I", "WT", "wildtype") ~ 0,
                         TRUE ~ NA))

  df$truth_value <- as.factor(df$truth_value)
  df$est_value <- as.factor(df$est_value)

  # --- compute confusion matrix ---
  cm <- yardstick::conf_mat(df, truth = truth_value, estimate = est_value)

  # --- compute yardstick metrics ---
  ys_metrics <- dplyr::bind_rows(
      yardstick::sens(df, truth = truth_value, estimate = est_value),
      yardstick::spec(df, truth = truth_value, estimate = est_value),
      yardstick::ppv(df, truth = truth_value, estimate = est_value),
      yardstick::npv(df, truth = truth_value, estimate = est_value),
      yardstick::accuracy(df, truth = truth_value, estimate = est_value),
      yardstick::kap(df, truth = truth_value, estimate = est_value),
      yardstick::f_meas(df, truth = truth_value, estimate = est_value)
    )

  # --- compute AMR-specific metrics (VME and ME) ---
  sensitivity <- ys_metrics$.estimate[ys_metrics$.metric == "sens"]
  specificity <- ys_metrics$.estimate[ys_metrics$.metric == "spec"]

  vme <- 1 - sensitivity # FN / (TP + FN)
  me <- 1 - specificity # FP / (TN + FP)

  amr_metrics <- tibble(
    .metric = c("VME", "ME"),
    .estimator = c("binary", "binary"),
    .estimate = c(vme, me)
  )

  # --- clean metrics format ---
  outcome_metrics <- dplyr::bind_rows(ys_metrics, amr_metrics) %>%
    select(metric = .metric, estimate = .estimate)

  return(list(metrics=outcome_metrics, matrix=cm))
}
