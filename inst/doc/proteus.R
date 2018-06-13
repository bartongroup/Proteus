## ----setup, include=FALSE-----------------------------------------------------
library(proteus)
library(knitr)
library(ggplot2)
options(width = 80)
knitr::opts_chunk$set(echo = TRUE)

## ----load_data, echo=FALSE----------------------------------------------------
library(proteusLabelFree, warn.conflicts=FALSE)
data(proteusLabelFree)

## ----load_data_dummy, eval=FALSE----------------------------------------------
#  library(proteusLabelFree)
#  data(proteusLabelFree)

## ----quick_start, eval=FALSE--------------------------------------------------
#  evidenceFile <- system.file("extdata", "evidence.txt.gz", package="proteusLabelFree")
#  metadataFile <- system.file("extdata", "metadata.txt", package="proteusLabelFree")
#  
#  evi <- readEvidenceFile(evidenceFile)
#  meta <- read.delim(metadataFile, header=TRUE, sep="\t")
#  pepdat <- makePeptideTable(evi, meta)
#  prodat <- makeProteinTable(pepdat)
#  prodat.med <- normalizeData(prodat)
#  res <- limmaDE(prodat.med)
#  plotVolcano_live(prodat.med, res)

## ----read_evidence, eval=FALSE------------------------------------------------
#  evidenceFile <- system.file("extdata", "evidence.txt.gz", package="proteusLabelFree")
#  evi <- readEvidenceFile(evidenceFile)

## ----head_evidence------------------------------------------------------------
head(evi)

## ----size_evidence------------------------------------------------------------
dim(evi)
format(object.size(evi), units="Mb")

## ----str_measure_columns------------------------------------------------------
str(measureColumns)

## ----str_evidence_columns-----------------------------------------------------
str(evidenceColumns)

## ----my_columns, eval=FALSE---------------------------------------------------
#  myColumns <- c(evidenceColumns, mz="m/z")

## ----rename_column, eval=FALSE------------------------------------------------
#  myColumns$protein <- "Leading Razor Protein"

## ----my_columns_evidence, eval=FALSE------------------------------------------
#  evi_mz <- readEvidenceFile(evidenceFile, columns=myColumns)

## ----evidence_column_names----------------------------------------------------
evidenceFile <- system.file("extdata", "evidence.txt.gz", package="proteusLabelFree")
evidence.columns <- readColumnNames(evidenceFile)
evidence.columns

## ----metadata-----------------------------------------------------------------
metadataFile <- system.file("extdata", "metadata.txt", package="proteusLabelFree")
meta <- read.delim(metadataFile, header=TRUE, sep="\t")

## ----show_metadata------------------------------------------------------------
meta

## ----metadata_samples---------------------------------------------------------
unique(evi$experiment)

## ----make_peptides, eval=FALSE------------------------------------------------
#  pepdat <- makePeptideTable(evi, meta)

## ----show_peptides------------------------------------------------------------
pepdat$tab[1:5, 1:5]

## ----summary_peptides---------------------------------------------------------
summary(pepdat)

## ----plot_peptide_count, fig.width=5, fig.height=4----------------------------
plotCount(pepdat)

## ----plot_detection_similarity, fig.width=5, fig.height=4---------------------
plotDetectionSimilarity(pepdat, bin.size = 0.02)

## ----plot_correlation_matrix, fig.width=6, fig.height=5-----------------------
plotDistanceMatrix(pepdat)

## ----plot_clusterin, fig.width=6, fig.height=5--------------------------------
plotClustering(pepdat)

## ----remove_bad_replicate-----------------------------------------------------
meta.clean <- meta[which(meta$sample  != 'BMO-7'),]

## ----make_clean_peptides, eval=FALSE------------------------------------------
#  pepdat.clean <- makePeptideTable(evi, meta.clean)

## ----clean_peptides_samples---------------------------------------------------
as.character(meta.clean$sample)

## ---- fig.width=6, fig.height=5-----------------------------------------------
plotClustering(pepdat.clean)

## ----make_proteins, eval=FALSE------------------------------------------------
#  prodat <- makeProteinTable(pepdat.clean)

## ----summary_proteins---------------------------------------------------------
summary(prodat)

## ----normalize_proteins-------------------------------------------------------
prodat.med <- normalizeData(prodat)

## ----normalize_proteins_quantiles---------------------------------------------
prodat.quant <- normalizeData(prodat, norm.fun=limma::normalizeQuantiles)

## ----intensity_distributions_1, fig.width=5, fig.height=4---------------------
plotSampleDistributions(prodat, title="Not normalized", fill="condition", method="violin")

## ----intensity_distributions_2, fig.width=5, fig.height=4---------------------
plotSampleDistributions(prodat.med, title="Median normalization", fill="condition", method="violin")

## ----intensity_distributions_3, fig.width=5, fig.height=4---------------------
plotSampleDistributions(prodat.quant, title="Quantile normalization", fill="condition", method="violin")

## ----plot_mv, fig.width=6, fig.height=4, warning=FALSE------------------------
plotMV(prodat.med, with.loess=TRUE)

## ----plot_clustering_proteins, fig.width=6, fig.height=5----------------------
plotClustering(prodat.med)

## ----protein_annotations------------------------------------------------------
luni <- lapply(as.character(prodat$proteins), function(prot) {
 if(grepl("sp\\|", prot)) {
   uniprot <- unlist(strsplit(prot, "|", fixed=TRUE))[2]
   c(prot, uniprot)
 }
})
ids <- as.data.frame(do.call(rbind, luni))
names(ids) <- c("protein", "uniprot")

## ----ids_head-----------------------------------------------------------------
head(ids)

## ----fetch_from_uniprot, eval=FALSE-------------------------------------------
#  annotations <- fetchFromUniProt(ids$uniprot, verbose=TRUE)

## ----annotations_head---------------------------------------------------------
head(annotations)

## ----merge_annotations--------------------------------------------------------
annotations.id <- merge(ids, annotations, by.x="uniprot", by.y="id")
annotations.id <- unique(annotations.id)

## ----annotate_proteins--------------------------------------------------------
prodat <- annotateProteins(prodat, annotations.id)

## ----evi_column_names---------------------------------------------------------
names(evi)

## ----peptide_table_modified_sequence, eval=FALSE------------------------------
#  pepdat.mod <- makePeptideTable(evi, meta, sequence.col = "modified_sequence")
#  prodat.mod <- makeProteinTable(pepdat.mod)

## ----protein_groups, eval=FALSE-----------------------------------------------
#  pepdat.group <- makePeptideTable(evi, meta, protein.col="protein_group")
#  prodat.group <- makeProteinTable(pepdat.group)

## ----peptide_aggregate_matrix-------------------------------------------------
evitab.example

## ----peptide_aggregate_default------------------------------------------------
aggregateSum(evitab.example)

## ----peptide_aggregate_maximum_function---------------------------------------
aggregateMax <- function(wp) {
  s <- apply(wp, 2, function(x) max(x, na.rm=TRUE))
  return(as.vector(s))
}

## ----peptide_aggregate_maximum_example----------------------------------------
aggregateMax(evitab.example)

## ----peptide_aggregate_maximum_create, eval=FALSE-----------------------------
#  pepdat.max <- makePeptideTable(evi, meta, aggregate.fun=aggregateMax)

## ----limma, warning=FALSE-----------------------------------------------------
res <- limmaDE(prodat.med, sig.level=0.05)

## ----show_limma_res-----------------------------------------------------------
head(res)

## ----show_significant_proteins------------------------------------------------
res[which(res$significant), c("protein", "logFC", "adj.P.Val")]

## ----limma_conditions, eval=FALSE---------------------------------------------
#  res <- limmaDE(prodat.med, conditions=c("1112", "BMO"))

## ----protein_detection--------------------------------------------------------
head(prodat$detect)

## ----missing_condition--------------------------------------------------------
only.1112 <- which(prodat$detect$`1112` & !prodat$detect$BMO)
only.BMO <- which(!prodat$detect$`1112` & prodat$detect$BMO)

## ----missing_condition_proteins-----------------------------------------------
as.character(prodat$proteins[only.1112])

## ----plot_fid, fig.width=4, fig.height=4, warning=FALSE-----------------------
plotFID(prodat.med)

## ----plot_volcano, fig.width=4, fig.height=4, warning=FALSE-------------------
plotVolcano(res)

## ----plot_pdist, fig.width=4, fig.height=4, warning=FALSE---------------------
plotPdist(res)

## ----plot_proteins, fig.width=4, fig.height=4, warning=FALSE------------------
plotIntensities(prodat.med, id='sp|P26263|PDC6_YEAST', log=TRUE)

## ----plot_prot_peptides, fig.width=7, fig.height=6, warning=FALSE-------------
plotProtPeptides(pepdat.clean, 'sp|P26263|PDC6_YEAST', prodat.med)

## ----shiny_volcano, eval=FALSE------------------------------------------------
#  plotVolcano_live(prodat.med, res)

## ----shiny_fid, eval=FALSE----------------------------------------------------
#  plotFID_live(prodat.med, res)

