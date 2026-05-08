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

AMRrulemakeR_env <- new.env()

utils::globalVariables(c(
	".estimate",
	".metric",
	":=",
	"Disk.I.ppv",
	"Disk.NWT.ppv",
	"Disk.R.ppv",
	"Disk.marker_list",
	"Disk.median",
	"Disk.n",
	"Disk.q25",
	"Disk.q75",
	"Gene symbol",
	"MIC.I.ppv",
	"MIC.NWT.ppv",
	"MIC.R.ppv",
	"MIC.marker_list",
	"MIC.median",
	"MIC.median_excludeRangeValues",
	"MIC.median_ignoreRanges",
	"MIC.n",
	"MIC.n_excludeRangeValues",
	"MIC.q25",
	"MIC.q25_excludeRangeValues",
	"MIC.q25_ignoreRanges",
	"MIC.q75",
	"MIC.q75_excludeRangeValues",
	"MIC.q75_ignoreRanges",
	"NWT.logreg.ci.lower",
	"NWT.logreg.ci.upper",
	"NWT.logreg.marker",
	"NWT.logregExt.ci.lower",
	"NWT.logregExt.ci.upper",
	"NWT.logregExt.marker",
	"NWT.solo.n",
	"NWT.soloExt.n",
	"PMID",
	"Prefix",
	"R.logreg.ci.lower",
	"R.logreg.ci.upper",
	"R.logreg.marker",
	"R.logregExt.ci.lower",
	"R.logregExt.ci.upper",
	"R.logregExt.marker",
	"R.solo.n",
	"R.soloExt.n",
	"ab",
	"ab_name",
	"allele",
	"as.sir",
	"assay_by_var",
	"breakpoint",
	"breakpoint condition",
	"breakpoint standard",
	"breakpoint_S",
	"breakpoint_disk",
	"breakpoint_mic",
	"breakpoints_source",
	"call_info",
	"category",
	"category_breakpoints_source",
	"category_diskPPV",
	"category_limitations_list",
	"category_logReg",
	"category_logRegExt",
	"category_micPPV",
	"category_note",
	"category_soloExtPPV",
	"category_soloPPV",
	"checkBreakpoints",
	"ci.lower",
	"ci.upper",
	"clinical category",
	"dat_note",
	"date_stamp",
	"disk",
	"disk.sources",
	"disk.sources.I",
	"disk.sources.R",
	"disk.sources.S",
	"drug",
	"drug class",
	"drug_class",
	"ecoff",
	"ecoff standard",
	"ecoff_disk",
	"ecoff_mic",
	"est_value",
	"evidence code",
	"evidence grade",
	"evidence limitations",
	"freq",
	"freq_n",
	"gene",
	"gene context",
	"gene_family",
	"getBreakpoints",
	"grade",
	"id",
	"label_list",
	"limitations_list",
	"map",
	"marker",
	"marker.label",
	"marker_count",
	"marker_list",
	"method",
	"mic",
	"mic.sources",
	"mic.sources.I",
	"mic.sources.R",
	"mic.sources.S",
	"mic.x.sources.I",
	"mic.x.sources.R",
	"mic.x.sources.S",
	"mic_note",
	"mo",
	"mutation",
	"n_total",
	"na.omit",
	"name",
	"name_txt",
	"node",
	"note",
	"organism",
	"organism_codes",
	"pct",
	"pheno",
	"phenotype",
	"phenotype_breakpoints_source",
	"phenotype_diskPPV",
	"phenotype_limitations_list",
	"phenotype_logReg",
	"phenotype_logRegExt",
	"phenotype_micPPV",
	"phenotype_note",
	"phenotype_soloExtPPV",
	"phenotype_soloPPV",
	"pivot_longer",
	"pivot_wider",
	"ppv.ranges",
	"pred",
	"pubmed_reference",
	"refgene_pubmed",
	"results",
	"rule curation note",
	"ruleID",
	"scale_fill_sir",
	"solo.sources",
	"solo.sources.I",
	"solo.sources.R",
	"solo.sources.S",
	"source_info",
	"str_count",
	"str_split",
	"tax_id",
	"taxid_bacteria",
	"total",
	"truth_value",
	"txid",
	"type",
	"unnest",
	"unnest_wider",
	"value",
	"variation type",
	"whitelisted_taxa",
	"x"
))

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