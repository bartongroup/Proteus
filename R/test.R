


makeProteins <- function(evi, int) {
  pep2prot <- select(evi, sequence, protein)
  pep2prot <- unique(pep2prot)

  # make sure order of peptides is as in intensity tables
  peptides <- data.frame(sequence=rownames(int[[1]]))
  peptides <- merge(peptides, pep2prot, by="sequence")
  proteins <- levels(as.factor(peptides$protein))

  for(condition in names(int)) {
    w <- int[[condition]]
    protint <- NULL
    for(prot in proteins) {
      sel <- which(peptides$protein == prot)
      if(length(sel) > 0)
      {
        wp <- w[sel,, drop=FALSE]  #ARGHHH! Took me forever to get it right
        medPep <- apply(wp, 1, function(x) {median(x, na.rm=TRUE)})  # median across replicates
        itop <- sort.int(medPep, decreasing=TRUE, index.return=TRUE)$ix[1:hifly] #index of top peptides
        row <- data.frame(protein=prot, t(colMeans(wp[itop,])))
        protint <- rbind(protint, row)
      }
    }


  }
}


makeOneProtein <- function(wp, hifly=3) {
  medPep <- apply(wp, 1, function(x) {median(x, na.rm=TRUE)})

}
