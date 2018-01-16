## ----setup, include=FALSE-----------------------------------------------------
library(proteus)
library(knitr)
library(dplyr)
library(ggplot2)
library(grid)
require(gridExtra)
library(dendextend)
options(width = 80)
knitr::opts_chunk$set(echo = TRUE)

## ----load_data----------------------------------------------------------------
library(proteusUnlabelled)
data(proteusUnlabelled)

## ----read_evidence, eval=FALSE------------------------------------------------
#  evidenceFile <- system.file("extdata", "evidence.txt.gz", package="proteusUnlabelled")
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

## ----my_columns_evidence, eval=FALSE------------------------------------------
#  evi_mz <- readEvidenceFile(evidenceFile, columns=myColumns)

## ----evidence_column_names----------------------------------------------------
evidenceFile <- system.file("extdata", "evidence.txt.gz", package="proteusUnlabelled")
evicols <- read.delim(evidenceFile, header=TRUE, sep="\t", check.names=FALSE, nrows=1)
names(evicols)

## ----metadata-----------------------------------------------------------------
metadataFile <- system.file("extdata", "metadata.txt", package="proteusUnlabelled")
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

