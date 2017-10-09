#' Proteomics metadata
#'
#' Metadata describing a proteomics experiment in two conditions with 7 replicates each.
#'
#' @docType data
#' @name xpmeta
#' @usage data(xpmeta)
#' @format A \code{data.frame} with 14 rows and 3 columns (sample, condition and replicate).
NULL

#' Proteomics metadata
#'
#' Metadata describing a proteomics experiment in two conditions with 7 replicates each. One "bad" replicate was removed.
#'
#' @docType data
#' @name xpmeta.clean
#' @usage data(xpmeta.clean)
#' @format A \code{data.frame} with 13 rows and 3 columns (sample, condition and replicate).
NULL

#' Evidence data
#'
#' Evidence data from a proteomics experiment in two conditions, 7 replicates each.
#'
#' @docType data
#' @name xpevi
#' @usage data(xpevi)
#' @format A data frame with 470318 rows and 10 columns
NULL

#' Peptide data
#'
#' Peptide data aggregated from evidence data \code{xpevi}. Contains all replicates.
#'
#' @docType data
#' @name xppepdat
#' @usage data(xppepdat)
#' @format A \code{proteusData} object
NULL

#' Clean peptide data
#'
#' Peptide data aggregated from evidence data \code{xpevi}. One "bad" replicate was removed.
#'
#' @docType data
#' @name xppepdat.clean
#' @usage data(xppepdat.clean)
#' @format A \code{proteusData} object
NULL


#' Proteins data
#'
#' Protein data aggregated from evidence data \code{xpevi}.
#'
#' @docType data
#' @name xpprodat
#' @usage data(xpprodat)
#' @format A \code{proteusData} object
NULL
