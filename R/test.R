


makeProteins <- function(evi, peptab, norm="median", min.peptides=1) {
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
      if(length(sel) >= min.peptides)
      {
        # ARGHHH! Took me forever to get it right
        # without drop=FALSE it will drop a dimension for one-row selection
        # and all downstream analysis goes to hell
        wp <- w[sel,, drop=FALSE]
        medPep <- apply(wp, 1, function(x) {median(x, na.rm=TRUE)})  # median across replicates
        nmed <- length(medPep)
        top <- min(c(nmed, hifly))
        itop <- sort.int(medPep, decreasing=TRUE, index.return=TRUE)$ix[1:top] #index of top peptides
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
  attr(protab, "logflag") <- FALSE
  attr(protab, "stats") <- intensityStats(protab)

  return(protab)
}

