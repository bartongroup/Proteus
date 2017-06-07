


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

plotMV <- function(intab, with.loess=FALSE, xmin=5, xmax=10, ymin=7, ymax=20) {
  meta <- attr(intab, "metadata")
  if(is.null(meta)) stop("No metadata found.")
  conditions <- meta$condition

  stats <- attr(intab, "stats")
  if(is.null(stats)) stats <- intensityStats(intab)

  P <- list()
  for(condition in conditions) {
    st <- stats[[condition]]
    df <- data.frame(mean=log10(st$mean), variance=log10(st$variance))
    P[[condition]] <- ggplot(df, aes(mean, variance)) + geom_point(alpha=0.3) +
      labs(x="log mean", y="log variance", title=condition) +
      xlim(xmin, xmax) +
      ylim(ymin, ymax)
    if(with.loess) {
      ls <- loess(variance ~ mean, data=df)
      x <- seq(from=min(na.omit(df$mean)), to=max(na.omit(df$mean)), by=0.05)
      pr <- predict(ls, x)
      pf <- data.frame(x=x, y=pr)
      P[[condition]] <- P[[condition]] + geom_line(data=pf, aes(x,y), color='yellow')
    }
  }
  grid.arrange(grobs=P, ncol=length(P))
}
