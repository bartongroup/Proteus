## ----setup, include=FALSE-----------------------------------------------------
library(proteus)
library(ggplot2)
library(knitr)
options(width = 80)
knitr::opts_chunk$set(echo = TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  vignette("unlabelled", package="proteus")

## ----load_data, echo=FALSE----------------------------------------------------
library(proteusSILAC, warn.conflicts=FALSE)
data(proteusSILAC)

## ----load_data_dummy, eval=FALSE----------------------------------------------
#  library(proteusSILAC)
#  data(proteusSILAC)

## ----measure_columns, cache=FALSE---------------------------------------------
measCols <- list(
  ML = "ratio_ml"
)

## ----evidence_columns, cache=FALSE--------------------------------------------
eviCols <- list(
  sequence = 'pep_sequence',
  modseq = 'modified_sequence',
  modifications = 'modifications',
  proteins = 'proteins',
  protein = 'leading_razor_protein',
  experiment = 'experiment',
  charge = 'charge',
  reverse = 'reverse',
  contaminant = 'potential_contaminant'
)

## ----read_evidence, eval=FALSE------------------------------------------------
#  evi <- readEvidenceFile(evidenceFile, measure.cols=measCols, data.cols=eviCols, zeroes.are.missing=FALSE)

## ----head_evidence------------------------------------------------------------
head(evi)

## ----metadata, eval=FALSE-----------------------------------------------------
#  metadataFile <- system.file("extdata", "metadata.txt", package="proteusSILAC")
#  meta <- read.delim(metadataFile, header=TRUE, sep="\t")

## ----show_metadata------------------------------------------------------------
meta

## ----make_peptides, eval=FALSE------------------------------------------------
#  pepdat <- makePeptideTable(evi, meta, measure.cols=measCols, aggregate.fun=aggregateMedian, experiment.type="SILAC")

## ----summary_peptides---------------------------------------------------------
summary(pepdat)

## ----plot_peptide_count, fig.width=5, fig.height=4----------------------------
plotCount(pepdat)

## ----make_proteins, eval=FALSE------------------------------------------------
#  prodat <- makeProteinTable(pepdat, aggregate.fun=aggregateMedian)

## ----summary_proteins---------------------------------------------------------
summary(prodat)

## ----normalize_proteins-------------------------------------------------------
prodat.norm <- normalizeData(prodat)

## ---- sample_dist, fig.width=4, fig.height=3----------------------------------
plotSampleDistributions(prodat, fill="replicate") + labs(title="Before")
plotSampleDistributions(prodat.norm, fill="replicate") + labs(title="After")

## ----limma_one_sample---------------------------------------------------------
res <- limmaRatioDE(prodat.norm, condition="T48")
res <- res[order(res$P.Value),]
head(res)

## ----volcano_plot, fig.width=4, fig.height=3----------------------------------
plotVolcano(res, binhex=FALSE)

## ----limma_two_samples--------------------------------------------------------
res2 <- limmaDE(prodat.norm)
res2 <- res2[order(res2$P.Value),]
head(res2)

## ----protein_example, fig.width=4, fig.height=3-------------------------------
plotIntensities(prodat.norm, id="P03372-3")

