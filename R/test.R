


makeProteins <- function(evi, peptab, norm="median") {
  meta <- attr(peptab, "metadata")
  normfac <- normalizingFactors(peptab, method=norm)
  peptab <- normalizeTable(peptab, normfac)

  pep2prot <- select(evi, sequence, protein)
  pep2prot <- unique(pep2prot)

  # make sure order of peptides is as in intensity tables
  peptides <- data.frame(sequence=rownames(peptab))
  peptides <- merge(peptides, pep2prot, by="sequence")
  proteins <- levels(as.factor(peptides$protein))

  intlist <- splitConditions(peptab)
  protlist <- list()
  for(condition in names(intlist)) {
    w <- intlist[[condition]]
    samples <- colnames(w)
    protint <- NULL
    for(prot in proteins) {
      sel <- which(peptides$protein == prot)
      if(length(sel) > 0)
      {
        wp <- w[sel,, drop=FALSE]  #ARGHHH! Took me forever to get it right
        medPep <- apply(wp, 1, function(x) {median(x, na.rm=TRUE)})  # median across replicates
        itop <- sort.int(medPep, decreasing=TRUE, index.return=TRUE)$ix[1:hifly] #index of top peptides
        row <- data.frame(protein=prot, t(colMeans(wp[itop,]))) # mean of the top peptides
        protint <- rbind(protint, row)
      }
    }
    colnames(protint) <- c("protein", samples)
    protlist[[condition]] <- protint
  }

  # dplyr join all tables
  protlist %>% Reduce(function(df1, df2) full_join(df1,df2, by="protein"), .) -> protab
  rownames(protab) <- protab$protein
  protab <- protab[,2:ncol(protab)]
  attr(protab, "metadata") <- meta
  attr(protab, "pep2prot") <- pep2prot
  attr(protab, "norm") <- norm
  attr(protab, "logflaf") <- FALSE


  # split again to calculate statistics
  intlist <- splitConditions(protab)
  stats <- list()
  for(condition in conditions) {
    w <- intlist[[condition]]
    m <- rowMeans(w, na.rm=TRUE)
    v <- apply(w, 1, function(v) sd(v, na.rm=TRUE)^2)
    stats[[condition]] <- data.frame(mean=m, variance=v)
  }
  attr(protab, "stats") <- stats

  return(protab)
}

