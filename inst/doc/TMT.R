## ----setup, include=FALSE-----------------------------------------------------
library(proteus)
library(knitr)
options(width = 80)
knitr::opts_chunk$set(echo = TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  vignette("unlabelled", package="proteus")

## ----load_data----------------------------------------------------------------
library(proteusTMT)
data(proteusTMT)

## ----measure_columns, cache=FALSE---------------------------------------------
measCols <- paste0("Reporter intensity ", 0:9)
names(measCols) <- paste0("reporter_", 0:9)

## ----measure_columns_object---------------------------------------------------
str(as.list(measCols))

## ----read_evidence, eval=FALSE------------------------------------------------
#  evidenceFile <- system.file("extdata", "evidence.txt.gz", package="proteusTMT")
#  evi <- readEvidenceFile(evidenceFile, measure.cols=measCols)

## ----head_evidence------------------------------------------------------------
head(evi)

## ----metadata, eval=FALSE-----------------------------------------------------
#  metadataFile <- system.file("extdata", "metadata.txt", package="proteusTMT")
#  meta <- read.delim(metadataFile, header=TRUE, sep="\t")

## ----show_metadata------------------------------------------------------------
meta

## ----make_peptides, eval=FALSE------------------------------------------------
#  pepdat <- makePeptideTable(evi, meta, measure.cols=measCols, fun.aggregate=median, experiment.type="TMT")

## ----summary_peptides---------------------------------------------------------
summary(pepdat)

## ----plot_peptide_count, fig.width=5, fig.height=4----------------------------
plotCount(pepdat)

## ----make_proteins, eval=FALSE------------------------------------------------
#  prodat <- makeProteinTable(pepdat, method="hifly", hifly=3)

## ----summary_proteins---------------------------------------------------------
summary(prodat)

## ----normalize_proteins-------------------------------------------------------
prodat.norm <- normalizeTMT(prodat)

## ---- sample_dist, fig.width=7, fig.height=4----------------------------------
plotSampleDistributions(prodat, fill="replicate") + labs(title="Before")
plotSampleDistributions(prodat.norm, log.scale=FALSE, fill="replicate") + labs(title="After")

## ----plot_clustering_proteins, fig.width=6, fig.height=5----------------------
plotClustering(prodat.norm)

## ----DE_C0_C4-----------------------------------------------------------------
res <- limmaDE(prodat.norm, conditions=c("C0", "C5"))

## ----p_value_dist, fig.width=5, fig.height=4----------------------------------
plotPdist(res)

## ----de_table_top-------------------------------------------------------------
res <- res[order(res$adj.P.Val),]
head(res)

## ----DE_example1, fig.width=8, fig.height=3-----------------------------------
prot <- as.character(res$protein[1])
plotIntensities(prodat.norm, id=prot)

