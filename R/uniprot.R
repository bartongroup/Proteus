#' Allowed UniProt columns
#'
#' A character vector containing possible columns to query with
#' \code{\link{fetchFromUniProt}}. Any selection of these names can be used as
#' \code{columns} parameter.
#'
#' @examples
#' allowedUniProtColumns
#'
#' @export
allowedUniProtColumns <- c("citation", "clusters", "comments", "domains", "domain", "ec",
                           "id", "entry name", "existence", "families", "features", "genes",
                           "go", "go-id", "interactor", "keywords", "last-modified", "length",
                           "organism", "organism-id", "pathway", "protein names", "reviewed",
                           "sequence", "3d", "version", "virus hosts")


#' Internal function to query UniProt
#'
#' @param url Full query URL
#' @param maxtry Maximum number of tries
#'
#' @return Query result
tryQuery <- function(url, maxtry=5) {
  for(i in 1:maxtry) {
    res <- tryCatch(read.delim(url, stringsAsFactors=FALSE), error=function(err) NULL)
    if(!is.null(res)) return(res)
    Sys.sleep(3)
  }
  stop(paste("UniProt not responding."))
}


#' Fetch annotations from UniProt
#'
#' @details
#'
#' For a given list of UniProt identifiers this function will bring back
#' annotations from UniProt servers. What information is downloaded is
#' controlled by the \code{columns} parameter. By default it fetches gene names
#' and protein name/description. The full list of available columns is in a
#' vector \code{allowedUniProtColumns}.
#'
#' The column names in the returned data frame are the same as in \code{columns}
#' parameter, unless alternative names are provided in parameter
#' \code{col.names}. The \code{id} column is added by default.
#'
#' @param unis A character vector with UniProt identifiers
#' @param columns Data columns requested (see \code{\link{allowedUniProtColumns}})
#' @param col.names How to name data columns in the returned data frame
#' @param batchsize Size of batch of proteins in a single query
#' @param verbose Logical, if true, query progress will be displayed
#'
#' @return A data frame with protein annotations.
#' @export
#'
#' @examples
#' library(proteusLabelFree)
#' data(proteusLabelFree)
#'
#' # Extract UniProt identifiers from protein IDs
#' unis <- sapply(as.character(prodat$proteins), function(prot) {
#'  s <- unlist(strsplit(prot, "|", fixed=TRUE))
#'  s[2]
#' })
#'
#' # Fetch first 100 annotations (for a quick example)
#' anno <- fetchFromUniProt(unis[1:100])
fetchFromUniProt <- function(unis, columns=c("genes", "protein names"),
                             col.names=NULL, batchsize=400, verbose=FALSE) {

  good.cols <- columns %in% allowedUniProtColumns
  if(sum(!good.cols) > 0) {
    bad.cols <- paste(columns[!good.cols], collapse=",")
    stop(paste("Columns", bad.cols, "are not allowed in UniProt queries. Here is a list of allowed columns:",
               paste(allowedUniProtColumns, collapse=", "), "."))
  }

  if(is.null(col.names)) {
    col.names <- columns
  } else {
    if(length(col.names) != length(columns)) stop("col.names has to be the same length as columns.")
  }

  columns <- c("id", columns)
  col.names <- c("id", col.names)
  colstr <- paste(columns, collapse=",")

  unis <- as.character(na.omit(unis))
  # test if IDs conform to UniProt accession format
  # regular expression from https://www.uniprot.org/help/accession_numbers
  acc.test <- grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$", unis, perl=TRUE)
  unis <- unis[acc.test]
  if(sum(!acc.test) > 0) {
    warning("Some identifiers do not conform to UniProt accession number format. Skipping.")
  }
  if(length(unis) == 0) {
    stop("No valid UniProt accession numbers found.")
  }

  url <- "http://legacy.uniprot.org/uniprot/"

  # split id list into batches of maximum size of batchsize
  batches <- split(unis, ceiling(seq_along(unis) / batchsize))
  nbatch <- length(batches)
  res <- lapply(seq_along(batches), function(i) {
    if(verbose) cat(paste("Batch", i, "out of", nbatch, "\n"))
    ids <- batches[[i]]
    qry <- paste(paste0("id:", ids), collapse="+or+")
    qurl <- URLencode(paste0(url, "?query=", qry, "&format=tab&columns=", colstr))
    tryQuery(qurl)
  })
  df <- do.call(rbind, res)
  colnames(df) <- col.names
  df
}


