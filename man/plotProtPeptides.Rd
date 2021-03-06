% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/func.R
\name{plotProtPeptides}
\alias{plotProtPeptides}
\title{Protein-peptide plot}
\usage{
plotProtPeptides(pepdat, protein, prodat = NULL, palette = cbPalette)
}
\arguments{
\item{pepdat}{Peptide \code{proteusData} object.}

\item{protein}{Protein name.}

\item{prodat}{(optional) protein \code{proteusData} object.}

\item{palette}{Palette of colours}
}
\value{
A \code{ggplot} object.
}
\description{
\code{plotProtPeptides} creates a plot consisting of two panels. The top
panel shows peptide log intensity. Each box represents one peptide, peptide
numbering follows alphabetical sequence order. The bottom panel shows sample
intensity. Each box represents one sample. White boxes show derived protein
intensities (if \code{prodat} is provided).
}
\examples{
library(proteusLabelFree)
data(proteusLabelFree)
prodat.med <- normalizeData(prodat)
plotProtPeptides(pepdat.clean, 'sp|P26263|PDC6_YEAST', prodat.med)

}
