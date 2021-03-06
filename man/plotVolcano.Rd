% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/func.R
\name{plotVolcano}
\alias{plotVolcano}
\title{Volcano plot}
\usage{
plotVolcano(
  res,
  bins = 80,
  xmax = NULL,
  ymax = NULL,
  marginal.histograms = FALSE,
  text.size = 12,
  show.legend = TRUE,
  plot.grid = TRUE,
  binhex = TRUE
)
}
\arguments{
\item{res}{Result table from  \code{\link{limmaDE}}.}

\item{bins}{Number of bins for binhex.}

\item{xmax}{Upper limit on x-axis. If used, the lower limit is -xmax.}

\item{ymax}{Upper limit on y-axis. If used, the lower limit is -ymax.}

\item{marginal.histograms}{A logical to add marginal histograms.}

\item{text.size}{Text size.}

\item{show.legend}{Logical to show legend (colour key).}

\item{plot.grid}{Logical to plot grid.}

\item{binhex}{Logical. If TRUE, a hexagonal density plot is made, otherwise it
is a simple point plot.}
}
\value{
A \code{ggplot} object.
}
\description{
\code{plotVolcano} makes a volcano plot from limma results. Uses
\code{\link{stat_binhex}} function from ggplot2 to make a hexagonal heatmap.
}
\examples{
library(proteusLabelFree)
data(proteusLabelFree)
prodat.med <- normalizeData(prodat)
res <- limmaDE(prodat.med)
plotVolcano(res)

}
