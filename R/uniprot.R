

tryQuery <- function(url, maxtry=5) {
  for(i in 1:maxtry) {
    res <- tryCatch(read.delim(url, stringsAsFactors=FALSE), error=function(err) NULL)
    if(!is.null(res)) return(res)
    Sys.sleep(3)
  }
  stop(paste("UniProt not responding."))
}


fetchFromUniProt <- function(unis, columns=c("genes", "protein names", "go", "go-id"), batchsize=400) {
  allowed.columns <- c("citation", "clusters", "comments", "domains", "domain", "ec",
                       "id", "entry name", "existence", "families", "features", "genes",
                       "go", "go-id", "interactor", "keywords", "last-modified", "length",
                       "organism", "organism-id", "pathway", "protein names", "reviewed",
                       "sequence", "3d", "version", "virus hosts")
  good <- columns %in% allowed.columns
  if(sum(!good) > 0) {
    bad <- paste(columns[!good], collapse=",")
    stop(paste("Columns", bad, "are not allowed in UniProt queries. Here is a list of allowed columns:", paste(allowed.columns, collapse=", "), "."))
  }

  url <- "http://www.uniprot.org/uniprot/"
  cols <- paste(c("id", columns), collapse=",")
  batches <- split(unis, ceiling(seq_along(unis) / batchsize))
  res <- lapply(batches, function(ids) {
    qry <- paste(ids, collapse="+or+")
    qurl <- URLencode(paste0(url, "?query=", qry, "&format=tab&columns=", cols))
    tryQuery(qurl)
  })
  d <- do.call(rbind, res)
}


