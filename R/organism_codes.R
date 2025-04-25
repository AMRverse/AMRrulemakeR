#' Organism Code Reference Table
#'
#' A dataset linking bacterial species names and AMR microorganisms (mo) to AMRrules
#' organism codes to use as prefixes for rules.
#'
#'
#' @format A tibble with 45 rows and 3 columns:
#' \describe{
#'   \item{Species}{The full species name, e.g. \code{"E. coli"}}
#'   \item{Prefix}{A 3-letter code representing the organism, e.g. \code{"ECO"}}
#'   \item{mo}{Microorganism ID as returned by \code{AMR::as.mo()}}
#' }
#'
#' @source Manually curated for this package.
#'
#' @examples
#' data(organism_codes)
#' head(organism_codes)
#' dplyr::filter(organism_codes, Prefix == "ECO")
#'
#' @usage data(organism_codes)
"organism_codes"
