library(proteus)
library(proteusUnlabelled)
library(UniProt.ws)

# Annotation proteins

### Extract UniProt IDs

# Create a named list Protein ID -> UniProt ID
unis <- sapply(as.character(prodat$proteins), function(prot) {
  s <- unlist(strsplit(prot, "|", fixed=TRUE))
  s[2]
})
head(unis)

#### Extract protein info from UniProt

# Find the right taxon ID
availableUniprotSpecies(pattern="cerevisiae")

# We have strain ATCC 204508 / S288c
# UniProt object
up <- UniProt.ws(559292)

# all UniProt IDs for the strain
ks <- keys(up, "UNIPROTKB")

# we can have a look at all available columns
cols <- columns(up)
cols

# we need, as a miniumum, gene names and protein names (descriptions)
# UNIPROTKB is the main key, corresponding to UniProt IDs
res <- select(up, unis, c("GENES", "PROTEIN-NAMES"), "UNIPROTKB")
rownames(res) <- names(unis)


