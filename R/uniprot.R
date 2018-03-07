#' A character vector containing possible columns to query using \code{fetchFromUniProt}
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
#' annotation from UniProt servers. What information is downloaded is controlled
#' by the \code{columns} parameter. By default it fetches gene names and protein
#' name/description. The full list of available columns is in a vector
#' \code{allowedUniProtColumns}.
#'
#' The column names in the returned data frame are the same as in \code{columns}
#' unless alterntive names are provided in parameter \code{col.names}.
#'
#' @param unis A character vector with UniProt identifiers
#' @param columns Data columns requested
#' @param col.names How to name data columns in the returned data frame
#' @param batchsize Size of batch of proteins in a single query
#' @param verbose Logical, if true, query progress will be displayed
#'
#' @return A data frame with protein annotations.
#' @export
#'
#' @examples
#' library(proteusUnlabelled)
#' data(proteusUnlabelled)
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
  good <- columns %in% allowedUniProtColumns
  if(sum(!good) > 0) {
    bad <- paste(columns[!good], collapse=",")
    stop(paste("Columns", bad, "are not allowed in UniProt queries. Here is a list of allowed columns:", paste(allowedUniProtColumns, collapse=", "), "."))
  }

  if(is.null(col.names)) {
    col.names <- columns
  } else {
    if(length(col.names) != length(columns)) stop("col.names has to be the same length as columns.")
  }
  columns <- c("id", columns)
  col.names <- c("id", col.names)

  unis.good <- as.character(na.omit(unis))
  cols <- paste(columns, collapse=",")

  url <- "http://www.uniprot.org/uniprot/"
  batches <- split(unis.good, ceiling(seq_along(unis.good) / batchsize))
  nbatch <- length(batches)
  res <- lapply(seq_along(batches), function(i) {
    if(verbose) cat(paste("Batch", i, "out of", nbatch, "\n"))
    ids <- batches[[i]]
    qry <- paste(paste0("id:", ids), collapse="+or+")
    qurl <- URLencode(paste0(url, "?query=", qry, "&format=tab&columns=", cols))
    tryQuery(qurl)
  })
  df <- do.call(rbind, res)
  colnames(df) <- col.names
  df
}


