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

## ---- eval=FALSE--------------------------------------------------------------
#  vignette("unlabelled", package="proteus")

## ----load_data----------------------------------------------------------------
library(proteusSILAC)
data(proteusSILAC)

## ----measure_columns, cache=FALSE---------------------------------------------
measCols <- list(
  ML = "ratio_ml"
)

## ----evidence_columns, cache=FALSE--------------------------------------------
eviCols <- evidenceColumns
eviCols$protein <- "Leading Razor Protein"
eviCols$modseq <- "Modified Sequence"
eviCols$contaminant <- "Contaminant"


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
#  evidenceFile <- system.file("extdata", "evidence.txt.gz", package="proteusSILAC")
#  #evidenceFile <- "/home/mgierlinski/projects/turnprot/pilot_data/evidence_pilot2.txt"
#  #evidenceFile <- "/home/mgierlinski/projects/noprot/evidence_e3_Phn.txt"
#  #evidenceFile <- "/home/mgierlinski/projects/temp/evi3.txt"
#  evi <- readEvidenceFile(evidenceFile, measure.cols=measCols, data.cols=eviCols, zeroes.are.missing=FALSE)
#  #evi <- readRDS("/home/mgierlinski/projects/temp/SILAC_evi_turnprot_selection.rds")

## ----head_evidence------------------------------------------------------------
head(evi)

## ----metadata, eval=FALSE-----------------------------------------------------
#  metadataFile <- system.file("extdata", "metadata.txt", package="proteusSILAC")
#  meta <- read.delim(metadataFile, header=TRUE, sep="\t")

## ----show_metadata------------------------------------------------------------
meta

## ----make_peptides, eval=FALSE------------------------------------------------
#  pepdat <- makePeptideTable(evi, meta, measure.cols=measCols, fun.aggregate=median, experiment.type="SILAC")

## ----summary_peptides---------------------------------------------------------
summary(pepdat)

## ----plot_peptide_count, fig.width=5, fig.height=4----------------------------
plotCount(pepdat)

## ----make_proteins, eval=FALSE------------------------------------------------
#  prodat <- makeProteinTable(pepdat, method="median")

## ----summary_proteins---------------------------------------------------------
summary(prodat)

## ----normalize_proteins-------------------------------------------------------
prodat.norm <- normalizeData(prodat)

## ---- sample_dist, fig.width=7, fig.height=4----------------------------------
plotSampleDistributions(prodat, fill="replicate") + labs(title="Before")
plotSampleDistributions(prodat.norm, fill="replicate") + labs(title="After")

## ----limma_one_sample---------------------------------------------------------
res <- limmaRatioDE(prodat.norm, condition="T48")
res <- res[order(res$P.Value),]
head(res)

