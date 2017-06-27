
#' Evidence columns
#'
#' \code{evidenceColumns} contains default columns to be read from the evidence file.
#' The names of list elements are used internally to reference evidence data.
#' @export
evidenceColumns <- list(
  sequence = 'Sequence',
  modseq = 'Modified sequence',
  modifications = 'Modifications',
  proteins = 'Proteins',
  protein = 'Leading razor protein',
  experiment = 'Experiment',
  charge = 'Charge',
  intensity = 'Intensity',
  reverse = 'Reverse',
  contaminant = 'Potential contaminant'
)

#' @import ggplot2
#' @import graphics
#' @import methods
#' @import stats
#' @import utils
simple_theme <- ggplot2::theme_bw() +
  ggplot2::theme(
    panel.border = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(colour = "black")
  )
simple_theme_grid <- ggplot2::theme_bw() +
  ggplot2::theme(
    panel.border = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(colour = "grey90"),
    panel.grid.minor = ggplot2::element_line(colour = "grey95"),
    axis.line = ggplot2::element_line(colour = "black")
  )


#' \code{proteusData} constructor
#'
#' @param tab A matrix with data, rows are peptides/proteins, columns are
#'   samples
#' @param metadata A data frame with metadata, needs at least two columns:
#'   "sample" and "condition"
#' @param content Either "peptides" or "proteins"
#' @param pep2prot A data frame with two columns: "sequence" and "protein". Can
#'   use "other" for non-proteus data.
#' @param peptides A character vector with peptide sequences
#' @param proteins A character vector with protein names
#' @param values A character vector with names of intensity and/or ratio columns
#' @param type Type of experiment: "unlabelled" or "SILAC"
#' @param npep A data frame with number of peptides per protein
#' @param pepseq Name of the sequence used, either "sequence" or "modseq"
#' @param hifly Number of high-flyers
#' @param min.peptides Integer, minimum number of peptides to combine into a
#'   protein
#' @param norm.fun Normalizing function
#' @param protein.method Name of the method with which proteins were quantified
#'
#' @return A \code{proteusData} object.
#' @export
proteusData <- function(tab, metadata, content, pep2prot, peptides, proteins, values,
                        npep=NULL, type="unlabelled", pepseq="sequence", hifly=NULL,
                        min.peptides=NULL, norm.fun=identity, protein.method=NULL) {
  stopifnot(
    ncol(tab) == nrow(metadata),
    is(tab, "matrix"),
    content %in% c("peptide", "protein", "other"),
    type %in% c("unlabelled", "SILAC"),
    pepseq %in% c("sequence", "modseq")
  )

  colnames(tab) <- metadata$sample
  if(content == "peptide") {
    stopifnot(nrow(tab) == length(peptides))
    rownames(tab) <- peptides
  } else if(content == "protein") {
    stopifnot(
      nrow(tab) == length(proteins),
      is.numeric(hifly),
      is.numeric(min.peptides),
      !is.null(protein.method)
    )
    rownames(tab) <- proteins
  }

  if(!is.null(npep)) {
    stopifnot(nrow(npep) == nrow(tab))
  }

  pdat <- list(
    tab = tab,
    metadata = metadata,
    content = content,
    conditions = levels(as.factor(metadata$condition)),
    type = type,
    values = values,
    pepseq = pepseq,
    hifly = hifly,
    min.peptides = min.peptides,
    norm.fun = deparse(substitute(norm.fun)),
    pep2prot = pep2prot,
    peptides = peptides,
    proteins = proteins,
    npep = npep,
    protein.method = protein.method
  )
  class(pdat) <- append(class(pdat), "proteusData")
  pdat$stats <- intensityStats(pdat)

  return(pdat)
}


#' Summary of \code{proteusData} object
#'
#' @param object \code{proteusData} object.
#' @param ... additional arguments affecting the summary produced.
#'
#' @export
#'
#' @examples
#' summary(prodat)
summary.proteusData <- function(object, ...) {
  cat("\n*** Basic statistics ***\n\n")
  cat(paste0("  content = ", object$content, "\n"))
  cat(paste0("  experiment type = ", object$type, "\n"))
  cat(paste0("  number of samples = ", nrow(object$metadata), "\n"))
  cat(paste0("  number of conditions = ", length(object$conditions), "\n"))
  cat(paste0("  number of ", object$content, "s = ", nrow(object$tab), "\n"))
  cat(paste0("  samples = ", paste0(object$metadata$sample, collapse = ", "), "\n"))
  cat(paste0("  conditions = ", paste0(object$conditions, collapse = ", "), "\n"))

  cat("\n*** Data processing ***\n\n")
  cat(paste0("  evidence columns used = ", paste0(object$values, collapse = ", "), "\n"))
  cat(paste0("  sequence = '", ifelse(object$pepseq == 'sequence', "Sequence", "Modified sequence"), "'\n"))
  cat(paste0("  normalization = ", object$norm.fun, "\n"))

  if(object$content == "protein") {
    cat("\n*** Protein data ***\n\n")
    cat(paste0("  quantification method = ", object$protein.method, "\n"))
    if(object$protein.method == "hifly") {
      cat(paste0("  number of high-flyers = ", object$hifly, "\n"))
    }
    cat(paste0("  minimum number of peptides per protein = ", object$min.peptides, "\n"))
  }
}


#' Read MaxQuant's evidence file
#'
#' \code{readEvidenceFile} reads MaxQuant's evidence file. Contaminants and
#' reverse sequences are filtered out. It can read only one intensity column
#' (specifed by \code{columns} parameter).
#'
#' @param file File name.
#' @param columns Named list with columns to read.
#' @return Data frame with selected columns from the evidence file.
#'
#' @examples
#' \dontrun{
#' evi <- readEvidenceFile("evidence.txt")
#' }
#'
#' @export
readEvidenceFile <- function(file, columns=evidenceColumns) {
  evi <- read.delim(file, header=TRUE, sep="\t", check.names=FALSE, as.is=TRUE, strip.white=TRUE)

  # check if all required columns exist
  missing <- NULL
  for(col in columns) {
    if(!(col %in% colnames(evi))) missing <- c(missing, paste0("'", col, "'"))
  }
  if(!is.null(missing))
    stop(paste0("Column(s) ", paste0(missing, collapse=", "), " not found in evidence file ", file))

  evi <- evi[, as.character(columns)]
  names(evi) <- names(columns)
  # sometimes there are only NAs and the condition doesn't work
  evi$reverse[is.na(evi$reverse)] = ''
  evi$contaminant[is.na(evi$contaminant)] = ''
  evi <- evi[which(evi$contaminant != '+' & evi$reverse != '+'),]
  evi <- evi[which(!is.na(evi$intensity)),]
}

#' Read peptides or proteinGroups file from MaxQuant
#'
#'
#' \code{readMaxQuantTable} reads a MaxQuant's output table (either peptides or
#' proteinGroups), extracts intensity data and creates a minimal
#' \code{proteusData} object.
#'
#' @param file Input file
#' @param meta Metadata
#' @param id Name of the column with identifiers (e.g. "Sequence" or "Protein
#'   ID")
#' @param columns A vector of intensity/ratio column names
#' @export
readMaxQuantTable <- function(file, meta, id, columns) {
  dat <- read.delim(file, header=TRUE, sep="\t", check.names=FALSE, as.is=TRUE, strip.white=TRUE)

  # check columns
  if(!(id %in% colnames(dat))) stop(paste0("Column '", id, "' not found in ", file))
  for(col in columns) {
    if(!(col %in% colnames(dat))) stop(paste0("Column '", col, "' not found in ", file))
  }

  tab <- as.matrix(dat[columns])
  tab[tab==0] <- NA
  tab[is.nan(tab)] <- NA
  colnames(tab) <- meta$sample
  rownames(tab) <- dat[[id]]

  # create pdat object
  pdat <- proteusData(tab, meta, "other", NULL, rownames(tab), NULL, NULL)
  class(pdat) <- append(class(pdat), "proteusData")
  return(pdat)
}


#' Create peptide table from evidence data
#'
#' \code{makePeptideTable} computes a peptide table and related data. Peptide
#' table is a matrix with columns corresponding to conditions and rows
#' corresponding to peptide sequences.
#'
#' @details
#'
#' In case of unlabelled experiments only one value is required, usually stored
#' in the "intensity" column. \code{makePeptideTable} will create a peptide
#' table with columns corresponding to samples in the metadata and rows
#' corresponding to peptide sequences. Each cell of this table is the sum of
#' peptide intensities (as defined by the default \code{fun.aggregate} function;
#' we do not recommend to change it).
#'
#' However, in SILAC and iTRAQ experiments there are multiple values, for
#' example "Ratio H/L", "Ratio M/L" and "Ratio H/M". \code{makePeptideTable}
#' will create a peptide table with columns corresponding to all combinations of
#' samples and values. The number of columns in the peptide table is the number
#' of samples in the metadata times the number of values. The cells of the
#' peptide table will be aggregate with the \code{fun.aggregate} function. We
#' recommend the median (use \code{fun.aggregate = median}) for this purpose. The
#' metadata attached to the output \code{proteusData} object will be adjusted
#' accordingly to contain new, expanded samples and conditions.
#'
#' In either case only samples from metadata are used, regardless of the content
#' of the evidence data. This makes selection of samples for downstream
#' processing easy: select only required rows in the metadata data frame.
#'
#' @param evi Evidence table created with \code{\link{readEvidenceFile}}.
#' @param meta Data frame with metadata. As a minimum, it should contain
#'   "sample" and "condition" columns.
#' @param pepseq A column name to identify peptides. Can be either "sequence" or
#'   "modseq".
#' @param values A vector of names (or one name) of column(s) with
#'   intensity/ratio data to be used.
#' @param fun.aggregate A function to aggregate pepetides with the same
#'   sequence/sample.
#' @param experiment.type Type of the experiment, "unlabelled" or "SILAC".
#' @return A \code{proteusData} object, containing peptide intensities and
#'   metadata.
#'
#' @examples
#' pepdat <- makePeptideTable(evi, meta)
#'
#' @export
makePeptideTable <- function(evi, meta, pepseq="sequence", values="intensity",
                             fun.aggregate=sum, experiment.type="unlabelled") {

  # check it ratios are in evidence data
  for(col in values) {
    if(!(col %in% names(evi))) stop(paste0("Column '", col, "' not found in evidence data."))
  }

  # zeroes in the evidence file are missing data
  evi[evi == 0] <- NA
  f <- function(x) fun.aggregate(x, na.rm=TRUE)

  # melt and recast evidence data
  eviMeasured <- evi[,c(pepseq, "experiment", values)]
  eviMelted <- reshape2::melt(eviMeasured, id.vars = c(pepseq, "experiment"))
  eviMelted$value <- as.numeric(eviMelted$value)   # integers do not work well in cast + median
  tab <- reshape2::dcast(eviMelted, paste0(pepseq, " ~ experiment + variable"), f)

  # convert to matrix, keep row names as 'peptides'
  peptides <- as.character(tab[,1])
  tab <- as.matrix(tab[,2:ncol(tab)])
  tab[tab == 0] <- NA
  rownames(tab) <- peptides


  # dcast created columns with merged names, these are our new samples
  # (e.g. samp1_ratio.HL, samp1_ratio.ML, samp2_ratio.HL, ...)
  # need to create new metadata to match these new samples

  # combine all pairs sample-value
  pairs <- expand.grid(meta$sample, values)
  colnames(pairs) <- c("sample", "value")
  # replicate naming convention from dcast (merge by "_")
  pairs$sample.new <- apply(pairs, 1, function(x) paste0(x, collapse="_"))
  # new metadata
  mt <- merge(meta, pairs, by="sample", sort=FALSE)
  selection <- mt$sample.new

  # clean-up metadata; when only one value, restore originals
  if(length(values) == 1) {
    mt <- meta
  } else {
    # need "sample" and "condition" corresponding to the new table
    # keep original columns for reference
    mt <- dplyr::rename(mt, sample.original = sample, condition.original = condition, sample = sample.new)
    # concatenate condition and value to create new condition
    mt$condition <- apply(mt[,c("condition.original", "value")], 1, function(x) paste0(x, collapse="_"))
  }

  # keep only columns from the metadata file
  # you can remove bad replicates
  tab <- tab[,as.character(selection)]
  colnames(tab) <- mt$sample

  # peptide to protein conversion
  pep2prot <- data.frame(sequence=evi$sequence, protein=evi$protein)
  pep2prot <- unique(pep2prot)
  rownames(pep2prot) <- pep2prot$sequence
  pep2prot <- pep2prot[peptides,]
  proteins <- levels(as.factor(pep2prot$protein))

  # create pdat object
  pdat <- proteusData(tab, mt, 'peptide', pep2prot, peptides, proteins, values,
                      type = experiment.type)

  return(pdat)
}



#' Make protein table
#'
#' \code{makeProteinTable} creates a protein table from the peptide table.
#' Protein intensities or ratios are aggregated from peptide data. There are
#' three ways of protein quantification: "hifly" and "sum" for intensity data
#' and "median" for SILAC data.
#'
#'
#' @details
#'
#' The "hifly" method is a follows. \enumerate{ \item For a given protein find
#' all corresponding peptides. \item In each replicate, order intensities from
#' the highest to the lowest. This is done separetly for each replicate. \item
#' Select n top rows of the ordered table. \item In each replicate, find the
#' mean of these three rows. This is the estimated protein intensity. } The
#' "sum" method simply adds all intensities for a given protein in each sample.
#' The "median", well, you can guess.
#'
#' @param pepdat A \code{proteusData} object containing peptide data, normally
#'   created by \code{\link{makePeptideTable}}
#' @param method Method to create proteins. Can be "hifly", "sum" or "median".
#' @param hifly The number of top peptides (high-flyers) to be used for protein
#'   intensity.
#' @param min.peptides Minimum number of peptides per protein.
#' @return A \code{proteusData} object containing protein intensities and
#'   metadata.
#'
#' @examples
#' prodat <- makeProteinTable(pepdat)
#'
#' @export
makeProteinTable <- function(pepdat, method="hifly", hifly=3, min.peptides=1) {
  if(!is(pepdat, "proteusData")) stop ("Input data must be of class proteusData.")
  if(!(method %in% c("hifly", "sum", "median"))) stop(paste0("Unknown method ", method))

  meta <- pepdat$metadata
  tab <- pepdat$tab

  protlist <- list()
  for(cond in pepdat$conditions) {
    w <- tab[,which(meta$condition == cond), drop=FALSE]
    samples <- colnames(w)
    protint <- NULL
    for(prot in pepdat$proteins) {
      sel <- which(pepdat$pep2prot$protein == prot)
      npep <- length(sel)
      if(npep >= min.peptides)
      {
        # ARGHHH! Took me forever to get it right
        # without drop=FALSE it will drop a dimension for one-row selection
        # and all downstream analysis goes to hell
        wp <- w[sel,, drop=FALSE]
        row <- makeProtein(wp, method, hifly)
        row <- data.frame(protein=prot, npep=npep, row)

        # this was old method using the same three peptides for all samples
        # but it doesn't work well, missing a lot of data
        # now (above) I pick three top peptides for each sample separately
        #medPep <- apply(wp, 1, function(x) {median(x, na.rm=TRUE)})  # median across replicates
        #nmed <- length(medPep)
        #top <- min(c(nmed, hifly))
        #itop <- sort.int(medPep, decreasing=TRUE, index.return=TRUE)$ix[1:top] #index of top peptides
        #row <- data.frame(protein=prot, t(colMeans(wp[itop,]))) # mean of the top peptides
        protint <- rbind(protint, row)
      }
    }
    colnames(protint) <- c("protein", "npep", samples)
    protlist[[cond]] <- protint
  }

  # dplyr join all tables
  protab <- Reduce(function(df1, df2) dplyr::full_join(df1,df2, by="protein"), protlist)

  proteins <- protab$protein
  npep <- data.frame(npep=protab$npep.x)     # join split npep into npep.x, npep.y, ... for conditions
  rownames(npep) <- proteins                 # a bit redundant, but might be useful
  protab <- as.matrix(protab[,as.character(meta$sample)])  # get rid of npep.y...

  prodat <- proteusData(protab, meta, "protein", pepdat$pep2prot, pepdat$peptides, pepdat$proteins, pepdat$values,
                        type = pepdat$type,
                        pepseq = pepdat$pepseq,
                        hifly = hifly,
                        min.peptides = min.peptides,
                        protein.method = method,
                        npep = npep)
  return(prodat)
}

#' Helper function, not exported
#'
#' Make protein from a selection of peptides
#' @param wp Intensity table with selection of peptides for this protein
#' @param method Method of creating proteins
#' @param hifly Number of high-fliers
makeProtein <- function(wp, method, hifly=3) {
  npep <- nrow(wp)
  if(npep == 0) stop("No peptides to aggregate.")
  cols <- colnames(wp)
  if(method == "hifly") {
    if(npep > 1) {
      wp[is.na(wp)] <- 0   # zeroes for sorting only
      sp <- as.matrix(apply(wp, 2, sort, decreasing=TRUE))
      sp[sp == 0] <- NA  # put NAs back
      ntop <- min(npep, hifly)
      row <- t(colMeans(sp[1:ntop,,drop=FALSE], na.rm=TRUE))
      row[is.nan(row)] <- NA    # colMeans puts NaNs where the column contains only NAs
      row <- as.data.frame(row)
    } else {
      row <- wp
    }
  } else if(method == "sum") {
    row <- as.data.frame(t(colSums(wp, na.rm=TRUE)))
    row[row==0] <- NA    # colSums puts zeroes where the column contains only NAs (!!!)
  } else if(method == "median") {
    row <- as.data.frame(t(apply(wp, 2, function(x) median(x, na.rm=TRUE))))
  } else stop(paste0("Unknown method ", method))


  colnames(row) <- cols
  return(row)
}

#' Normalize columns of a matrix to medians
#'
#' \code{normalizeMedian} normalizes the columns of a matrix to have the same
#' medians. It should be used with \code{\link{normalizeData}} function.
#'
#' @param tab Data frame with raw intensities. Normally, this is a \code{tab}
#'   field in the \code{proteusData} object (see examples below).
#' @return Normalized matrix.
#'
#' @examples
#' prodat.norm <- normalizeData(prodat, norm.fun=normalizeMedian)
#'
#' @export
normalizeMedian <- function(tab) {
  norm.fac <- apply(tab, 2, function(x) {median(x, na.rm=TRUE)})
  norm.fac <- norm.fac / mean(norm.fac, na.rm=TRUE)
  tab <- t(t(tab) / norm.fac)
  return(tab)
}


#' Normalize proteus data
#'
#' \code{normalizeData} normalizes the intensity table in a \code{proteusData}
#' object, using the provided normalizing function.
#'
#' @param pdat A \code{proteusData} object with peptide or protein intensities.
#' @param norm.fun A normalizing function.
#' @return A \code{proteusData} object with normalized intensities.
#'
#' @details The normalizng function, specified by \code{norm.fun} needs to
#'   normalize columns of a numerical matrix. The input is a matrix and the
#'   output is a normalized matrix. The default value points to
#'   \code{\link{normalizeMedian}}, Proteus's function normalizing each column
#'   to its median. Other functions can be used, for example
#'   \code{normalizeQuantiles} from the \code{limma} package.
#'
#' @examples
#' prodat.norm <- normalizeData(prodat, norm.fun=normalizeMedian)
#'
#' @export
normalizeData <- function(pdat, norm.fun=normalizeMedian) {
  if(!is(pdat, "proteusData")) stop ("Input data must be of class proteusData.")
  if(class(norm.fun) != "function") stop ("'norm.fun' has to be a function.")

  pdat$tab <- norm.fun(pdat$tab)
  pdat$norm.fun <- deparse(substitute(norm.fun))
  pdat$stats <- intensityStats(pdat)  # need to recalculate stats!
  return(pdat)
}


#' Plot correlation matrix
#'
#' \code{plotCorrelationMatrix} plots a correlation matrix for peptide or
#' protein data.
#'
#' @param pdat Peptide or protein \code{proteusData} object.
#'
#' @examples
#' plotCorrelationMatrix(pepdat)
#'
#' @export
plotCorrelationMatrix <- function(pdat) {
  if(!is(pdat, "proteusData")) stop ("Input data must be of class proteusData.")
  corr.mat <- cor(pdat$tab, use="complete.obs")
  gplots::heatmap.2(corr.mat, trace="none", density.info="none", dendrogram="none", Rowv=FALSE, Colv=FALSE, key.xlab = "Correlation coefficient")
}


#' Plot clustering dendrogram
#'
#' \code{plotClustering} plots a dendrogram of intensity data, using
#' hierarchical clustering.
#'
#' @param pdat Peptide or protein \code{proteusData} object.
#'
#' @examples
#' plotClustering(pepdat)
#'
#' @export
plotClustering <- function(pdat) {
  if(!is(pdat, "proteusData")) stop ("Input data must be of class proteusData.")
  corr.mat <- cor(pdat$tab, use="complete.obs")
  dis <- as.dist(1 - corr.mat)  # dissimilarity matrix
  dend <- as.dendrogram(hclust(dis))
  colors_to_use <- as.numeric(pdat$metadata$condition)
  colors_to_use <- colors_to_use[order.dendrogram(dend)]
  dendextend::labels_colors(dend) <- colors_to_use
  plot(dend)
}


#' Plot peptide count per sample
#'
#' \code{plotPeptideCount} makes a plot of peptide count per sample.
#'
#' @param pdat Peptide \code{proteusData} object.
#' @param x.text.size Size of text on the x-axis.
#' @return A plot of the number of peptides detected in each sample.
#'
#' @examples
#' plotPeptideCount(pepdat)
#'
#' @export
plotPeptideCount <- function(pdat, x.text.size=10){
  if(!is(pdat, "proteusData")) stop ("Input data must be of class proteusData.")
  meta <- pdat$meta
  pep.count <- apply(pdat$tab, 2, function(x) sum(!is.na(x)))
  med.count <- median(pep.count)
  df <- data.frame(x=meta$sample, y=pep.count, condition=meta$condition)
  ggplot(df, aes(x=x,y=y,fill=condition)) +
    geom_col() +
    simple_theme +
    scale_y_continuous(expand = c(0,0)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=x.text.size)) +
    labs(x='Sample', y='Peptide count') +
    labs(title = paste0("Median peptide count = ", med.count)) +
    theme(plot.title=element_text(hjust=0, size=12))
}


#' Plot distribution of intensities/ratios for each sample
#'
#' \code{plotSampleDistributions} makes a boxplot with intensity/ratio distribution for each sample.
#'
#' @param pdat A \code{proteusData} object with peptide/protein intensities.
#' @param title Title of the plot.
#' @param method "box", "violin" or "dist"
#' @param x.text.size Text size in value axis
#' @param x.text.angle Text angle in value axis
#' @param vmin Lower bound on log value
#' @param vmax Upper bound on log value
#' @param n.grid.rows Number of rows in the grid of facets
#' @param hist.bins Number of bins in histograms
#' @param logbase Base of the logarithm which will be applied to data
#' @param fill A metadata column to use for the fill of boxes
#' @param colour A metadata column to use for the outline colour of boxes
#' @param hline Logical, if true a horizontal line at zero is added
#'
#' @export
#'
#' @examples
#' plotSampleDistributions(prodat)
#' plotSampleDistributions(normalizeData(prodat))
plotSampleDistributions <-
function(pdat, title="", method="dist", x.text.size=7, n.grid.rows=3, hist.bins=100,                                  x.text.angle=90, vmin=as.numeric(NA), vmax=as.numeric(NA), logbase=10, fill=NULL,
         colour=NULL, hline=FALSE) {
  m <- reshape2::melt(pdat$tab, varnames=c("ID", "sample"))
  mt <- data.frame(pdat$metadata, row.names = pdat$metadata$sample)
  if(!is.null(fill)) m[['fill']] <- mt[m$sample, fill]
  if(!is.null(colour)) m[['colour']] <- mt[m$sample, colour]

  m$value <- log(m$value, base=logbase)

  if(method == "box" | method == "violin") {
    g <- ggplot(m, aes(x=sample, y=value)) +
      simple_theme +
      ylim(vmin, vmax) +
      theme(axis.text.x = element_text(angle = x.text.angle, hjust = 1, vjust=0.5, size=x.text.size)) +
      labs(title=title, x="sample", y=paste0("log", logbase, " value"))
      if(hline) g <- g + geom_hline(yintercept=0, colour='grey')
      if(method=="box") g <- g + geom_boxplot(outlier.shape=NA, na.rm=TRUE)
      if(method=="violin") g <- g + geom_violin(na.rm=TRUE, draw_quantiles=c(0.25, 0.5, 0.75), scale="width")
  } else if (method == "dist") {
    g <- ggplot(m, aes(x=value)) +
      #simple_theme +
      geom_histogram(bins=hist.bins) +
      facet_wrap(~sample, nrow=n.grid.rows) +
      xlim(vmin, vmax)
  } else stop("Wrong method.")

  if(!is.null(fill)) g <- g + aes(fill=fill) + scale_fill_discrete(name=fill)
  if(!is.null(colour)) g <- g + aes(colour=colour) + scale_color_discrete(name=colour)
  g
}


#' Statistics for an intensity table
#'
#' \code{intensityStats} calculates the mean, variance and number of good data
#' points for peptide or protein intensities.
#'
#' @param pdat Peptide or protein \code{proteusData} object.
#' @return A data frame with several statistics.
#'
#' @examples
#' stats <- intensityStats(prodat)
#'
#' @export
intensityStats <- function(pdat) {
  if(!is(pdat, "proteusData")) stop ("Input data must be of class proteusData.")

  meta <- pdat$metadata

  stats <- NULL
  for(cond in pdat$conditions) {
    w <- pdat$tab[,which(meta$condition == cond), drop=FALSE]
    m <- rowMeans(w, na.rm=TRUE)
    m[which(is.nan(m))] <- NA
    v <- apply(w, 1, function(v) sd(v, na.rm=TRUE)^2)
    ngood <- apply(w, 1, function(v) sum(!is.na(v)))
    stats <- rbind(stats, data.frame(condition=cond, mean=m, variance=v, ngood=ngood))
  }
  return(stats)
}


#' Plot mean-variance relationship
#'
#' \code{plotMV} makes a plot with variance of log-intensity vs mean of log-intensity.
#'
#' @param pdat Peptide or protein \code{proteusData} object.
#' @param with.loess Logical. If true, a loess line will be added.
#' @param bins Number of bins for binhex
#' @param xmin Lower limit on x-axis
#' @param xmax Upper limit on x-axis
#' @param ymin Lower limit on y-axis
#' @param ymax Upper limit on y-axis
#'
#' @examples
#' plotMV(prodat, with.loess=TRUE)
#'
#' @export
plotMV <- function(pdat, with.loess=FALSE, bins=80, xmin=5, xmax=10, ymin=7, ymax=20) {
  if(!is(pdat, "proteusData")) stop ("Input data must be of class proteusData.")
  meta <- pdat$metadata
  if(is.null(meta)) stop("No metadata found.")

  stats <- pdat$stats
  if(is.null(stats)) stats <- intensityStats(pdat)
  stats <- stats[which(!is.na(stats$mean) & !is.na(stats$variance) & stats$mean > 0 & stats$variance > 0),]
  stats$mean <- log10(stats$mean)
  stats$variance <- log10(stats$variance)
  protnum <- as.data.frame(table(stats$condition))  #number of proteins in each condition
  colnames(protnum) <- c("condition", "n")

  # has to be calculated for each condition separately
  if(with.loess) {
    ldf <- NULL
    for(cond in pdat$conditions)
    {
      st <- stats[which(stats$condition == cond),]
      ls <- loess(variance ~ mean, data=st)
      x <- seq(from=min(na.omit(st$mean)), to=max(na.omit(st$mean)), by=0.01)
      pr <- predict(ls, x)
      ldf <- rbind(ldf, data.frame(condition=cond, x=x, y=pr))
    }
  }

  g <- ggplot(stats, aes(x=mean, y=variance)) +
    simple_theme_grid +
    xlim(xmin, xmax) +
    ylim(ymin, ymax) +
    facet_wrap(~condition) +
    stat_binhex(bins=bins) +
    scale_fill_gradientn(colours=c("green","yellow", "red"), name = "count", na.value=NA) +
    geom_text(data=protnum, aes(x=xmin+0.5, y=ymax, label=paste0("n = ", n)))
  if(with.loess) g <- g + geom_line(data=ldf, aes(x,y), color='black')
  return(g)
}

#' Plot protein(s)
#'
#' \code{plotProteins} makes a plot with protein intensity as a function of the
#' condition and replicate. When multiple proteins are entered, the mean and
#' standard error is plotted.
#'
#' @param pdat Protein \code{proteusData} object.
#' @param protein Protein name (string) or a vector with protein names.
#' @param log Logical. If set TRUE a logarithm of intensity is plotted.
#' @param ymin Lower bound for y-axis
#' @param ymax Upper bound for y-axis
#' @param text.size Text size
#' @param point.size Point size
#' @param title Title of the plot (defaults to protein name)
#'
#' @examples
#' plotProteins(prodat, "sp|P16522|CDC23_YEAST", log=TRUE, title="MNT2")
#'
#' @export
plotProteins <- function(pdat, protein=protein, log=FALSE, ymin=as.numeric(NA), ymax=as.numeric(NA),
                         text.size=12, point.size=3, title=NULL) {
  # without 'as.numeric' it returns logical NA (!!!)

  if(!is(pdat, "proteusData")) stop ("Input data must be of class proteusData.")

  sel <- which(pdat$proteins %in% protein)
  if(length(sel) > 0 && sel > 0) {
    E <- if(log) log10(pdat$tab[sel,,drop=FALSE]) else pdat$tab[sel,,drop=FALSE]
    e <- colMeans(E, na.rm=TRUE)
    s <- sapply(E, function(x) sd(x, na.rm=TRUE)/sqrt(length(x)))
    n <- length(sel)

    if(is.null(title)) {
      title <- ifelse(n == 1, protein, paste0("selection of ", n, " proteins."))
    }

    meta <- pdat$metadata
    p <- data.frame(
      expr = e,
      lo = e - s,
      up = e + s,
      condition = factor(meta$condition, levels=unique(as.factor(meta$condition))),
      replicates = factor(meta$replicate)
    )
    # define shapes
    p$shape <- rep(21, length(p$expr))
    p$shape[which(p$expr==0)] <- 24
    pd <- position_dodge(width = 0.15)
    #
    if(is.na(ymin) && !log) ymin=0
    ylab <- ifelse(log, "log10 intensity", "Intensity")
    # colour-blind friendly palette
    cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    ggplot(p, aes(x=condition, y=expr, ymin=lo, ymax=up, fill=replicates, shape=shape)) +
      simple_theme_grid +
      theme(text = element_text(size=text.size), legend.position = "none") +
      ylim(ymin, ymax) +
      {if(n > 1) geom_errorbar(position=pd, width = 0.1)} +
      geom_point(position=pd, size=point.size) +
      scale_shape_identity() +  # necessary for shape mapping
      #scale_fill_manual(values=cbPalette) +
      labs(x = 'Condition', y = ylab, title=title)
  }
}

#' Differential expression with limma
#'
#' \code{limmaDE} is a simple wrapper around limma differential expression. It
#' performs differential expression on the intensity table.
#'
#' @param pdat Protein \code{proteusData} object.
#' @param formula A string with a formula for building the linear model.
#' @return limma output from eBays. See limma documentation for more details.
#'   The output from this function can be used with \code{\link{limmaTable}} to
#'   create a table with differential expression results.
#'
#' @examples
#' ebay <- limmaDE(prodat, formula="~condition")
#' res <- limmaTable(prodat, ebay)
#'
#' @export
limmaDE <- function(pdat, formula="~condition") {
  if(!is(pdat, "proteusData")) stop ("Input data must be of class proteusData.")
  design <- model.matrix(as.formula(formula), pdat$metadata)
  tab <- log10(pdat$tab)
  fit <- limma::lmFit(tab, design)
  ebay <- limma::eBayes(fit)
}

#' Create differential expression result
#'
#' \code{limmaTable} creates a table with differential expressioni results,
#' using an object created with \code{\link{limmaDE}}
#'
#' @param pdat Protein \code{proteusData} object.
#' @param ebay Output from  \code{\link{limmaDE}}.
#' @param column Which column should be used to extract data. The default value
#'   is "condition".
#' @return A data frame with DE results.
#'
#' @examples
#' ebay <- limmaDE(prodat, formula="~condition")
#' res <- limmaTable(prodat, ebay)
#'
#' @export
limmaTable <- function(pdat, ebay, column="condition") {
  if(!is(pdat, "proteusData")) stop ("Input data must be of class proteusData.")
  # levels from the column (e.g. conditions)
  levs <- levels(factor(pdat$metadata[[column]]))
  # coef is the column name + the last level, e.g. "conditionWT"
  coef <- paste0(column, levs[-1])
  res <- limma::topTable(ebay, coef=coef, adjust="BH", sort.by="none", number=1e9)
  res <- cbind(protein=rownames(res), res)
  rownames(res) <- c()
  return(res)
}

#' Fold-change intensity diagram
#'
#' \code{plotFID} makes a log fold change versus log sum intensity plot, usually
#' known as MA plot.
#'
#' @param pdat Protein \code{proteusData} object.
#' @param pair A two-element vector containing the pair of conditions to use.
#'   Can be skipped if there are only two conditions.
#' @param pvalue An optional vector with corresponding p-values to be used with
#'   an interactive plotly plot.
#' @param bins Number of bins for binhex.
#' @param marginal.histograms A logical to add marginal histograms.
#' @param xmin Lower limit on x-axis.
#' @param xmax Upper limit on x-axis.
#' @param ymax Upper limit on y-axis. If used, the lower limit is -ymax.
#' @param text.size Text size.
#' @param show.legend Logical to show legend (colour key).
#' @param plot.grid Logical to plot a grid.
#' @param binhex Logical. If TRUE, a hexagonal densit plot is made, otherwise it
#'   is a simple point plot.
#'
#' @examples
#' plotFID(prodat)
#'
#' @export
plotFID <- function(pdat, pair=NULL, pvalue=NULL, bins=80, marginal.histograms=FALSE,
                   xmin=NULL, xmax=NULL, ymax=NULL, text.size=12, show.legend=TRUE, plot.grid=TRUE,
                   binhex=TRUE) {
  if(is.null(pair)) pair <- pdat$conditions
  if(length(pair) != 2) stop("Need exactly two conditions. You might need to specify pair parameter.")
  if(!is(pdat, "proteusData")) stop ("Input data must be of class proteusData.")

  stats <- pdat$stats
  c1 <- pair[1]
  c2 <- pair[2]
  m1 <- log10(stats[which(stats$condition == c1),]$mean)
  m2 <- log10(stats[which(stats$condition == c2),]$mean)
  d <- data.frame(x = (m1 + m2) / 2, y = m2 - m1)
  d$id <- rownames(pdat$tab)
  n <- length(m1)
  if(is.null(pvalue)) {
    d$p <- rep('NA', n)
  } else {
    d$p <- signif(pvalue, 2)
  }

  title <- paste0(c1, ":", c2)
  g <- ggplot(d, aes(x=x, y=y)) +
    {if(plot.grid) simple_theme_grid else simple_theme} +
    {if(binhex) stat_binhex(bins=bins, show.legend=show.legend) else geom_point(aes(text=id))}+
    scale_fill_gradientn(colours=c("green","yellow", "red"), name = "count",na.value=NA) +
    #geom_point(na.rm=TRUE, alpha=0.6, size=0.5, aes(text=paste(id, "\nP = ", p))) +
    geom_abline(colour='red', slope=0, intercept=0) +
    labs(title=title, x=paste0(c1, '+', c2), y=paste0(c2, '-', c1)) +
    theme(text = element_text(size=text.size))
  if(!is.null(xmin) && !is.null(xmax)) g <- g + scale_x_continuous(limits = c(xmin, xmax), expand = c(0, 0))
  if(!is.null(ymax) ) g <- g + scale_y_continuous(limits = c(-ymax, ymax), expand = c(0, 0))
  if(marginal.histograms) g <- ggExtra::ggMarginal(g, size=10, type = "histogram", xparams=list(bins=100), yparams=list(bins=50))
  return(g)
}

#' Plot p-value distribution
#'
#' \code{plotPdist} makes a plot with distribution of raw p-values, obtained by
#' \code{\link{limmaDE}}.
#'
#' @param res Output table from \code{\link{limmaTable}}.
#' @param text.size Text size.
#' @param plot.grid Logical to plot grid.
#' @param bin.size Bin size for the histogram.
#'
#' @examples
#' ebay <- limmaDE(prodat)
#' res <- limmaTable(prodat, ebay)
#' plotPdist(res)
#'
#' @export
plotPdist <- function(res, bin.size=0.02, text.size=12, plot.grid=TRUE) {
  ggplot(res, aes(P.Value, ..density..)) +
    {if(plot.grid) simple_theme_grid else simple_theme} +
    geom_histogram(breaks=seq(0, 1, bin.size), colour='blue') +
    labs(x='P-value', y='Density') +
    theme(text = element_text(size=text.size))
}

#'Volcano plot
#'
#'\code{plotVolcano} makes a volcano plot from limma results. Uses
#'\code{\link{stat_binhex}} function from ggplo2 to make a hexagonal heatmap.
#'
#'@param res Result table from  \code{\link{limmaTable}}.
#'@param bins Number of bins for binhex.
#'@param xmax Upper limit on x-axis. If used, the lower limit is -xmax.
#'@param ymax Upper limit on y-axis. If used, the lower limit is -ymax.
#'@param text.size Text size.
#'@param show.legend Logical to show legend (colour key).
#'@param plot.grid Logical to plot grid.
#'@param binhex Logical. If TRUE, a hexagonal densit plot is made, otherwise it
#'  is a simple point plot.
#'
#' @examples
#' ebay <- limmaDE(prodat)
#' res <- limmaTable(prodat, ebay)
#' plotVolcano(res)
#'
#'@export
plotVolcano <- function(res, bins=80, xmax=NULL, ymax=NULL, text.size=12, show.legend=TRUE,
                        plot.grid=TRUE, binhex=TRUE) {
  g <- ggplot(res, aes(logFC, -log10(P.Value))) +
    {if(plot.grid) simple_theme_grid else simple_theme} +
    {if(binhex) stat_binhex(bins=bins, show.legend=show.legend) else geom_point(aes(text=protein))} +
    scale_fill_gradientn(colours=c("green","yellow", "red"), name = "count",na.value=NA) +
    geom_vline(colour='red', xintercept=0) +
    theme(text = element_text(size=text.size))
    # labs(title=title, x=paste0(c1, '+', c2), y=paste0(c2, '-', c1))

    if(!is.null(xmax)) g <- g + scale_x_continuous(limits = c(-xmax, xmax), expand = c(0, 0))
    if(!is.null(ymax) ) g <- g + scale_y_continuous(limits = c(0, ymax), expand = c(0, 0))
  return(g)
}


#' Protein-peptide plot
#'
#' \code{plotProtPeptides} creates a plot consisting of two panels. The top
#' panel shows peptide log intensity. Each box represents one peptide, peptide
#' numbering follows alphabetical sequence order. The bottom panel shows sample
#' intensity. Each box represents one sample. White boxes show derived protein
#' intensities (if \code{prodat} is provided).
#'
#' @param pepdat Peptide \code{proteusData} object.
#' @param protein Protein name.
#' @param prodat (optional) protein \code{proteusData} object.
#'
#' @examples
#' plotProtPeptides(pepdat, "sp|P16522|CDC23_YEAST", prodat = prodat)
#'
#' @export
plotProtPeptides <- function(pepdat, protein, prodat=NULL) {
  if(!is(pepdat, "proteusData")) stop ("Input data must be of class proteusData.")
  tab <- normalizeMedian(pepdat$tab)

  # all peptides for this protein
  selprot <- which(pepdat$pep2prot$protein == protein)
  if(length(selprot) == 0) stop(paste0("Protein '", protein, "' not found."))

  peps <- pepdat$pep2prot[selprot,'sequence']
  mat <- as.matrix(tab[peps,])
  dat <- reshape2::melt(mat, varnames=c("peptide", "sample"))
  levels(dat$sample) <- pepdat$metadata$sample # melt loses order of sample levels
  dat$pepnum <- sprintf("%02d", as.numeric(dat$peptide))  # convert sequences into numbers
  dat$intensity <- log10(dat$value)

  # add condition column (is there a simpler way?)
  s2c <- dplyr::select(pepdat$metadata, condition)
  rownames(s2c) <- pepdat$metadata$sample
  dat$condition <- s2c[dat$sample,]

  # add protein intensity column (is there a simpler way)
  if(!is.null(prodat)) {
    p2p <- data.frame(int=as.numeric(prodat$tab[protein,]))
    rownames(p2p) <- colnames(prodat$tab)
    dat$prot.intensity <- log10(p2p[dat$sample,])
  }

  g1 <- ggplot(dat, aes(x=pepnum, y=intensity, fill=condition)) +
    geom_boxplot(outlier.shape = NA)  +
    geom_jitter(width=0, size=0.5) +
    facet_wrap(~condition) +
    theme(legend.position="none") +
    labs(x="Peptide", y="log intensity", title=protein)
  g2 <- ggplot(dat, aes(x=sample, y=intensity, fill=condition)) +
    geom_boxplot(outlier.shape = NA)  +
    geom_jitter(width=0, size=0.5) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) +
    theme(legend.position="none") +
    labs(x="Sample", y="log intensity")
  if(!is.null(prodat)) g2 <- g2 + geom_point(aes(x=sample, y=prot.intensity), shape=22, size=3, fill='white')
  gridExtra::grid.arrange(g1, g2, ncol=1)
}

