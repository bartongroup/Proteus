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
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#' Evidence columns
#'
#' \code{evidenceColumns} contains default columns to be read from the evidence
#' file. The names of list elements are used internally to reference evidence
#' data.
#'
#' @examples
#' str(evidenceColumns)
#'
#' @export
evidenceColumns <- list(
  sequence = 'Sequence',
  modseq = 'Modified sequence',
  modifications = 'Modifications',
  proteins = 'Proteins',
  protein = 'Leading razor protein',
  experiment = 'Experiment',
  charge = 'Charge',
  reverse = 'Reverse',
  contaminant = 'Potential contaminant'
)

#' Measure columns
#'
#' \code{measureColumns} contains measurement columns from evidence file. In
#' case of unlabelled data, there is only one column: Intensity.
#'
#' @examples
#' str(measureColumns)
#'
#' @export
measureColumns <- list(
  intensity = 'Intensity'
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
#' @param measures A character vector with names of intensity and/or ratio
#'   columns
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
proteusData <- function(tab, metadata, content, pep2prot, peptides, proteins, measures,
                        npep=NULL, type="unlabelled", pepseq="sequence", hifly=NULL,
                        min.peptides=NULL, norm.fun=identity, protein.method=NULL) {
  stopifnot(
    ncol(tab) == nrow(metadata),
    is(tab, "matrix"),
    content %in% c("peptide", "protein", "other"),
    type %in% c("unlabelled", "SILAC", "TMT"),
    pepseq %in% c("sequence", "modseq")
  )

  # make sure to keep correct order of samples
  metadata$sample <- factor(metadata$sample, levels=metadata$sample)

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

  # number of replicates in each condition
  cnd <- metadata$condition
  cnd.fac <- as.factor(cnd)
  nrep <- tapply(cnd, cnd.fac, length) # nice trick

  pdat <- list(
    tab = tab,
    metadata = metadata,
    measures = measures,
    content = content,
    conditions = levels(cnd.fac),
    nrep = nrep,
    type = type,
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
  pdat$detect <- goodData(pdat)

  return(pdat)
}


#' Summary of \code{proteusData} object
#'
#' @param object \code{proteusData} object.
#' @param ... additional arguments affecting the summary produced.
#'
#' @return Text output with object summary
#'
#' @examples
#' library(proteusUnlabelled)
#' data(proteusUnlabelled)
#' summary(prodat)
#'
#' @export
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
  cat(paste0("  evidence columns used = ", paste0(object$measures, collapse = ", "), "\n"))
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
#' reverse sequences are filtered out.
#'
#' @details
#'
#' There are two parameters controlling which columns are read from the evidence
#' file. Parameter \code{measure.cols} selects columns with measurements: these
#' are intensities (unlabelled, TMT) or ratios (Silac). In the simplest case of
#' unlabelled data, there is only one measure column: "Intensity". Parameter
#' \code{data.columns} selects all other columns read from the evidence file.
#' There are two default lists, supplied with the package, appropriate for an
#' unlabelled experiment, \code{measureColumns} and \code{evidenceColumns}.
#'
#'
#' @param file File name.
#' @param measure.cols Named list with measure columns to read.
#' @param data.cols Named list with other columns to read (in addition to
#'   measure columns).
#' @return Data frame with selected columns from the evidence file.
#'
#' @examples
#' library(proteusUnlabelled)
#' evidenceFile <- system.file("extdata", "evidence.txt.gz", package="proteusUnlabelled")
#' evi <- readEvidenceFile(evidenceFile)
#'
#' @export
readEvidenceFile <- function(file, measure.cols=measureColumns, data.cols=evidenceColumns) {

  # check if all required columns are in the evidence file
  evi.cols <- read.delim(file, header=TRUE, sep="\t", check.names=FALSE, as.is=TRUE, strip.white=TRUE, nrows = 1)
  columns <- c(data.cols, measure.cols)
  missing <- NULL
  for(col in columns) {
    if(!(col %in% colnames(evi.cols))) missing <- c(missing, paste0("'", col, "'"))
  }
  if(!is.null(missing))
    stop(paste0("Column(s) ", paste0(missing, collapse=", "), " not found in evidence file ", file))

  # read and process evidence file
  evi <- read.delim(file, header=TRUE, sep="\t", check.names=FALSE, as.is=TRUE, strip.white=TRUE)
  evi <- evi[, as.character(columns)]
  names(evi) <- names(columns)
  # sometimes there are only NAs and the condition doesn't work
  evi$reverse[is.na(evi$reverse)] = ''
  evi$contaminant[is.na(evi$contaminant)] = ''
  evi <- evi[which(evi$contaminant != '+' & evi$reverse != '+'),]

  # replace zeroes with NAs in measure columns
  evi.meas <- evi[,names(measure.cols), drop=FALSE]
  evi.meas[evi.meas == 0] <- NA
  evi[,names(measure.cols)] <- evi.meas
  rm(evi.meas)

  # remove rows that have only NAs in measure columns
  not.empty <- which(rowSums(!is.na(evi[,names(measure.cols), drop=FALSE])) > 0)
  evi <- evi[not.empty,]
}


#' Create peptide table from evidence data
#'
#' \code{makePeptideTable} computes a peptide table and related data. Peptide
#' table is a matrix with columns corresponding to conditions and rows
#' corresponding to peptide sequences.
#'
#' @details
#'
#' The evidence file contains a column called "Experiment" and one or more
#' columns with measure values. In case of unlabelled experiment there is only
#' one measure column: "Intensity". In case of TMT experiment there are several
#' measure columns, usually called "Reporter intensity 0", "Reporter intensity
#' 1", and so on. \code{makePeptideTable} will combine "Experiment" and measure
#' columns in a way defined by the metadata (parameter \code{meta}). The name of
#' the combined column will be named using "sample" column in metadata.
#'
#' The result is a \code{proteusData} object containing a table with rows
#' corresponding to peptides and columns corresponding to samples (as defined in
#' metadata). Each cell of the table is an aggregated measure values over all
#' evidence entries corresponding to the given sequence and experiment. How
#' these measure values are aggregated is controlled by the parameter
#' \code{fun.aggregate}. We recommend using sum for unlabelled and TMT
#' experiments and median for SILAC experiments.
#'
#' Only samples from metadata are used, regardless of the content of the
#' evidence data. This makes selection of samples for downstream processing
#' easy: select only required rows in the metadata data frame.
#'
#' @param evi Evidence table created with \code{\link{readEvidenceFile}}.
#' @param meta Data frame with metadata. As a minimum, it should contain
#'   "sample" and "condition" columns.
#' @param pepseq A column name to identify peptides. Can be either "sequence" or
#'   "modseq".
#' @param measure.cols A named list of measure columns; should be the same as
#'   used in \code{\link{readEvidenceFile}}
#' @param fun.aggregate A function to aggregate pepetides with the same
#'   sequence/sample.
#' @param experiment.type Type of the experiment, "unlabelled" or "SILAC".
#' @return A \code{proteusData} object, containing peptide intensities and
#'   metadata.
#'
#' @examples
#' library(proteusUnlabelled)
#' data(proteusUnlabelled)
#' pepdat <- makePeptideTable(evi, meta)
#'
#' @export
makePeptideTable <- function(evi, meta, pepseq="sequence", measure.cols=measureColumns,
                             fun.aggregate=sum, experiment.type="unlabelled") {

  # check if measure.cols are in the evidence file
  measures <- names(measure.cols)
  for(col in measures) {
    if(!(col %in% names(evi))) stop(paste0("Column '", col, "' not found in evidence data."))
  }

  # test columns in metadata
  for(col in c("experiment", "measure", "sample", "condition")) {
    if(!(col %in% names(meta))) stop(paste0("Column '", col, "' not found in metadata."))
  }

  # make sure to keep correct order of samples
  meta$sample <- factor(meta$sample, levels=meta$sample)

  # zeroes in the evidence file are missing data
  evi[evi == 0] <- NA
  f <- function(x) fun.aggregate(x, na.rm=TRUE)

  # melt and recast evidence data
  # conflict between reshape and reshape2 requires a direct call to melt.data.frame (note three colons!)
  eviMelted <- reshape2:::melt.data.frame(evi,
                                          id.vars = c(pepseq, "experiment"),
                                          measure.vars=measures,
                                          variable.name="measure"
                                          )
  eviMelted$value <- as.numeric(eviMelted$value)   # integers do not work well in cast + median

  # translate "experiment" and "measure" into "sample" from metadata
  m2s <- meta$sample
  names(m2s) <- paste0(meta$experiment, ".", meta$measure)

  # recover raw evidence names, they are in metadata. merge with experiment
  eviMelted$expmes <- paste0(eviMelted$experiment, ".", measure.cols[eviMelted$measure])
  # select only experiment/measure combinations present in metadata
  eviMelted <- eviMelted[which(eviMelted$expmes %in% names(m2s)),]

  # translate
  eviMelted$sample <- m2s[eviMelted$expmes]

  # cast into table: sample vs sequence
  tab <- reshape2::dcast(eviMelted, paste0(pepseq, " ~ sample"), f)

  # convert to matrix, keep row names as 'peptides'
  peptides <- as.character(tab[,1])
  tab <- as.matrix(tab[,2:ncol(tab)])
  tab[tab == 0] <- NA
  rownames(tab) <- peptides

  # peptide to protein conversion
  pep2prot <- data.frame(sequence=evi[[pepseq]], protein=evi$protein)
  pep2prot <- unique(pep2prot)
  rownames(pep2prot) <- pep2prot$sequence
  pep2prot <- pep2prot[peptides,]
  proteins <- levels(as.factor(pep2prot$protein))

  # create pdat object
  pdat <- proteusData(tab, meta, 'peptide', pep2prot, peptides, proteins,
                      as.character(measure.cols), type = experiment.type)

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
#' @param ncores Number of cores for parallel processing
#' @return A \code{proteusData} object containing protein intensities and
#'   metadata.
#'
#' @examples
#' library(proteusUnlabelled)
#' data(proteusUnlabelled)
#' prodat <- makeProteinTable(pepdat.clean, ncores=2)
#'
#' @export
makeProteinTable <- function(pepdat, method="hifly", hifly=3, min.peptides=1, ncores=4) {
  if(!is(pepdat, "proteusData")) stop ("Input data must be of class proteusData.")
  if(!(method %in% c("hifly", "sum", "median"))) stop(paste0("Unknown method ", method))

  meta <- pepdat$metadata
  tab <- pepdat$tab

  # make sure to keep correct order of samples
  # this should be OK in peptide table, but just in case
  meta$sample <- factor(meta$sample, levels=meta$sample)

  protlist <- list()
  for(cond in pepdat$conditions) {
    w <- tab[,which(meta$condition == cond), drop=FALSE]
    samples <- colnames(w)
    protcond <- parallel::mclapply(pepdat$proteins, function(prot) {
      sel <- which(pepdat$pep2prot$protein == prot)
      npep <- length(sel)
      if(npep >= min.peptides)
      {
        wp <- w[sel,, drop=FALSE]
        row <- makeProtein(wp, method, hifly)
      } else {
        row <- as.data.frame(t(rep(NA, length(samples))))
        names(row) <- samples
      }
      row <- data.frame(protein=prot, npep=npep, row)
      return(row)
    }, mc.cores=ncores)
    protint <- do.call(rbind, protcond)

    colnames(protint) <- c("protein", "npep", samples)
    protlist[[cond]] <- protint
  }

  # dplyr join all tables
  protab <- Reduce(function(df1, df2) dplyr::full_join(df1,df2, by="protein"), protlist)

  # remove empty rows (happens when min.peptides > 1)
  protab <- protab[which(rowSums(!is.na(protab)) > 0), ]

  proteins <- protab$protein
  npep <- data.frame(npep=protab$npep.x)     # join split npep into npep.x, npep.y, ... for conditions
  rownames(npep) <- proteins                 # a bit redundant, but might be useful
  protab <- as.matrix(protab[,as.character(meta$sample)])  # get rid of npep.y...

  prodat <- proteusData(protab, meta, "protein", pepdat$pep2prot, pepdat$peptides,
                        proteins, pepdat$measures,
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
#' @return A data frame row with aggregated protein intensities
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


#' Annotate proteins
#'
#' \code{annotateProteins} adds protein annotations to a \code{proteusData}
#' object.
#'
#' @details
#'
#' The only information about proteins proteus package extracts from the
#' evidence file are protein identifiers. These can be in various forms,
#' depending on how MaxQuant was ran. In order to annotate them, two steps are
#' required.
#'
#' First, the user needs to create a data frame linking protein identifiers (as
#' in \code{pdat$proteins} vector) to some metadata. This data frame has to
#' contain a column called \code{protein} with the identifiers and any
#' additional columns with annotations, e.g., UniProt IDs, protein names, gene
#' names, domains, GO-terms and so on. These data can be obtained from UniProt.
#'
#' Once the annotation table is created, it can be merged into the
#' \code{proteusData} object using \code{annotateProteins} function. The order
#' of identifiers in the annotation table is not important. Also, not all
#' proteins have to be present. The merge will add only the matching proteins.
#' As a result, this function returns a \code{proteusData} object with
#' \code{annotation} field added. It contains an annotation data frame with rows
#' corresponding to the internal list of proteins. This means the rows of this
#' data frame correspond one-to-one to the rows in \code{pdat$tab} and
#' \code{pdat$detect}, if pdat contains proteins (not peptides).
#'
#' @param pdat A \code{proteusData} object containing protein data
#' @param annotation A data frame with a column \code{protein} containing
#'   protein identifiers, as in \code{pdat$proteins}
#'
#' @return A \code{proteusData} with annotation field added.
#'
#' @examples
#' library(proteusUnlabelled)
#' data(proteusUnlabelled)
#' library(UniProt.ws)
#'
#' # Extract UniProt identifiers from protein IDs
#' unis <- sapply(as.character(prodat$proteins), function(prot) {
#'  s <- unlist(strsplit(prot, "|", fixed=TRUE))
#'  s[2]
#' })
#'
#' up <- UniProt.ws(559292)
#'
#' # Fetch data from UniProt
#' unitab <- select(up, unis, c("GENES", "PROTEIN-NAMES"), "UNIPROTKB")
#' unitab$protein <- names(unis)
#' unitab <- unitab[!is.na(unitab$`PROTEIN-NAMES`),]
#' # simplify protein names
#' unitab$name <- gsub("\\s\\(.*$", "", unitab$`PROTEIN-NAMES`, perl=TRUE)
#'
#' prodat <- annotateProteins(prodat, unitab)
#'
#' @export
annotateProteins <- function(pdat, annotation) {
  if(!is(pdat, "proteusData")) stop ("Input data must be of class proteusData.")
  if(!is.data.frame(annotation)) stop ("Annotation needs to be a data frame.")
  if(is.null(annotation$protein)) stop ("Annotation table has to have a 'protein' column.")

  nover <- length(intersect(pdat$proteins, annotation$protein))
  if(nover == 0) stop("No overlap between data and annotation. Nothing to annotate.")

  mrg <- data.frame(protein=pdat$proteins)
  mrg <- merge(mrg, annotation, by="protein", all.x=TRUE)
  pdat$annotation <- mrg

  cat(paste("Annotated", nover, "proteins.\n"))

  return(pdat)
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
#' library(proteusUnlabelled)
#' data(proteusUnlabelled)
#' normtab <- normalizeMedian(prodat$tab)
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
#' library(proteusUnlabelled)
#' data(proteusUnlabelled)
#' prodat.med <- normalizeData(prodat, norm.fun=normalizeMedian)
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



#' RAS procedure for CONSTANd normalization (not exported)
#'
#' @param K Input matrix
#' @param max.iter Maximum number of iterations
#' @param eps Convergence limit
#'
#' @return Matrix normalized so row and column means equal 1/n (n - number of
#'   columns)
RAS <- function(K, max.iter=50, eps=1e-5) {
  n <- ncol(K)
  m <- nrow(K)

  # ignore rows with only NAs
  good.rows <- which(rowSums(!is.na(K)) > 0)
  K <- K[good.rows, ]

  cnt <- 1
  repeat {
    row.mult <- 1 / (n * rowMeans(K, na.rm=TRUE))
    K <- K * row.mult
    err1 <- 0.5 * sum(abs(colMeans(K, na.rm=TRUE) - 1/n))
    col.mult <- 1 / (n * colMeans(K, na.rm=TRUE))
    K <- t(t(K) * col.mult)
    err2 <- 0.5 * sum(abs(rowMeans(K, na.rm=TRUE) - 1/n))
    cnt <- cnt + 1
    if(cnt > max.iter || (err1 < eps && err2 < eps)) break
  }

  # reconstruct full table
  KF <- matrix(NA, nrow=m, ncol=n)
  KF[good.rows, ] <- K

  return(KF)
}

#' Normalize protein/peptide data using CONSTANd normalization
#'
#' @details
#'
#' \code{normalizeTMT} implements CONSTANd algorithm from Maes et al. (2016)
#' \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4974351/pdf/zjw2779.pdf}.
#' It normalizes TMT table, for each experiment separately. After normalization
#' each row shows the precentage of total row intensity and each column is
#' normalized to their mean value. Hence, both row and column means are equal
#' 1/n, where n is the number of columns (reporters).
#'
#' @param pdat A \code{proteusData} object with peptide or protein intensities.
#' @param max.iter Maximum number of iterations for the RAS procedure.
#' @param eps Convergence limit for the RAS procedure.
#'
#' @return A \code{proteusData} with normalized data.
#'
#' @examples
#' library(proteusTMT)
#' data(proteusTMT)
#' prodat.norm <- normalizeTMT(prodat)
#'
#' @export
normalizeTMT <- function(pdat, max.iter=50, eps=1e-5) {
  for(ex in levels(as.factor(pdat$meta$experiment))) {
    colsel <- which(pdat$metadata$experiment == ex)
    pdat$tab[, colsel] <- RAS(pdat$tab[, colsel])
    pdat$norm.fun <- "CONSTANd"
    pdat$stats <- intensityStats(pdat)  # need to recalculate stats!

  }
  return(pdat)
}

#' Plot distance matrix
#'
#' \code{plotDistanceMatrix} plots a distance matrix for peptide or protein
#' data.
#'
#' Computes a distance matrix and plots it as a shade-plot. The default distance
#' is Pearson's correlation coefficient. Other methods will be introduced later.
#'
#' @param pdat Peptide or protein \code{proteusData} object.
#' @param distance A method to calculate distance.
#' @param text.size Text size on axes
#' @return A \code{ggplot} object.
#'
#' @examples
#' library(proteusUnlabelled)
#' data(proteusUnlabelled)
#' plotDistanceMatrix(pepdat)
#'
#' @export
plotDistanceMatrix <- function(pdat, distance=c("correlation"), text.size=10) {
  if(!is(pdat, "proteusData")) stop ("Input data must be of class proteusData.")
  distance <- match.arg(distance)

  corr.mat <- cor(pdat$tab, use="complete.obs")
  m <- reshape2::melt(corr.mat, varnames=c("Sample1", "Sample2"))
  m$Sample1 <- factor(m$Sample1, levels=pdat$metadata$sample)
  m$Sample2 <- factor(m$Sample2, levels=pdat$metadata$sample)
  #gplots::heatmap.2(corr.mat, trace="none", density.info="none", dendrogram="none", Rowv=FALSE, Colv=FALSE, key.xlab = "Correlation coefficient")
  ggplot(m, aes(x=Sample1, y=Sample2)) +
    geom_tile(aes(fill=value)) +
    #scale_fill_manual(values=palette) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=text.size),
      axis.text.y = element_text(size=text.size)
    ) +
    labs(x='Sample', y='Sample', fill="Correlation")
}


#' Plot clustering dendrogram
#'
#' \code{plotClustering} plots a dendrogram of intensity data, using
#' hierarchical clustering.
#'
#' @param pdat Peptide or protein \code{proteusData} object.
#' @return Creates a plot
#'
#' @examples
#' library(proteusUnlabelled)
#' data(proteusUnlabelled)
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


#' Plot peptide or protein count per sample
#'
#' \code{plotCount} makes a plot of peptide/protein count per sample.
#'
#' @param pdat A \code{proteusData} object.
#' @param x.text.size Size of text on the x-axis.
#' @param palette Palette of colours
#' @return A plot (\code{ggplot} object) of the number of peptides/proteins
#'   detected in each sample.
#'
#' @examples
#' library(proteusUnlabelled)
#' data(proteusUnlabelled)
#' plotCount(pepdat)
#'
#' @export
plotCount <- function(pdat, x.text.size=10, palette=cbPalette){
  if(!is(pdat, "proteusData")) stop ("Input data must be of class proteusData.")
  meta <- pdat$metadata
  entry.count <- apply(pdat$tab, 2, function(x) sum(!is.na(x)))
  med.count <- median(entry.count)
  df <- data.frame(x=meta$sample, y=entry.count, condition=meta$condition)
  g <- ggplot(df, aes(x=x,y=y,fill=condition)) +
    geom_col(colour='grey60') +
    simple_theme +
    scale_y_continuous(expand = c(0,0)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=x.text.size)) +
    labs(x="Sample", y="Count") +
    labs(title = paste0("Median count = ", med.count)) +
    theme(plot.title=element_text(hjust=0, size=12))
  if(nlevels(as.factor(df$condition)) <= length(palette)) g <- g + scale_fill_manual(values=palette)
  g
}

#' Jaccard similarity
#'
#' Computes Jaccard similarity between "detections" in two vectors of the same
#' length. A detection is a value, as opposed to NA.
#'
#' @param x A vector
#' @param y A vector
#'
#' @return Jaccard similarity between x and y
#'
#' @examples
#' library(proteusUnlabelled)
#' data(proteusUnlabelled)
#' sim <- jaccardSimilarity(pepdat$tab[,1], pepdat$tab[,2])
#'
#' @export
#'
jaccardSimilarity <- function(x, y) {
  stopifnot(length(x) == length(y))

  intersection <- which(!is.na(x) & !is.na(y))
  union <- which(!is.na(x) | !is.na(y))
  jaccard <- ifelse(length(union) > 0, length(intersection) / length(union), 0)
  return(jaccard)
}

#' Detection Jaccard similarity
#'
#' \code{plotDetectionSimilarit} plots a distribution of (peptide) detection
#' similarities between samples. Jaccard similarity between two sets is defined
#' as the size of the intersection divided by the size of the union. This plot
#' can be used to assess quality of data.
#'
#' @param pdat A \code{proteusData} object with peptides (or proteins).
#' @param text.size Text size.
#' @param plot.grid Logical to plot grid.
#' @param bin.size Bin size for the histogram.
#' @param hist.colour Colour of the histogram.
#' @return A \code{ggplot} object.
#'
#' @examples
#' library(proteusUnlabelled)
#' data(proteusUnlabelled)
#' plotDetectionSimilarity(pepdat)
#'
#'
#' @export
#'
plotDetectionSimilarity <- function(pdat, bin.size=0.01, text.size=12, plot.grid=TRUE, hist.colour='grey30') {
  if(!is(pdat, "proteusData")) stop ("Input data must be of class proteusData.")

  n <- ncol(pdat$tab)
  # indices of upper triangular: all pair-wise combinations of columns in tab
  pairs <- which(upper.tri(matrix(1, n, n)) == TRUE, arr.ind=TRUE)
  # Jaccard simliarity for each pair of columns in tab
  sim <- apply(pairs, 1, function(i) jaccardSimilarity(pdat$tab[,i[1]], pdat$tab[,i[2]]))

  ggplot(data.frame(sim), aes(sim, ..density..)) +
  {if(plot.grid) simple_theme_grid else simple_theme} +
    geom_histogram(breaks=seq(0, 1, bin.size), colour=hist.colour, fill=hist.colour) +
    labs(x='Jaccard similarity', y='Density') +
    theme(text = element_text(size=text.size))
}


#' Plot distribution of intensities/ratios for each sample
#'
#' \code{plotSampleDistributions} makes a boxplot with intensity/ratio
#' distribution for each sample.
#'
#' @param pdat A \code{proteusData} object with peptide/protein intensities.
#' @param title Title of the plot.
#' @param method "box", "violin" or "dist"
#' @param palette Palette of colours
#' @param x.text.size Text size in value axis
#' @param x.text.angle Text angle in value axis
#' @param vmin Lower bound on log value
#' @param vmax Upper bound on log value
#' @param n.grid.rows Number of rows in the grid of facets
#' @param hist.bins Number of bins in histograms
#' @param log.scale Logical, to plot in logarithmic scale
#' @param log.base Base of the logarithm which will be applied to data
#' @param fill A metadata column to use for the fill of boxes
#' @param colour A metadata column to use for the outline colour of boxes
#' @param hline Logical, if true a horizontal line at zero is added
#'
#' @return A \code{ggplot} object.
#'
#' @export
#'
#' @examples
#' library(proteusUnlabelled)
#' data(proteusUnlabelled)
#' plotSampleDistributions(prodat)
#' plotSampleDistributions(normalizeData(prodat))
#'
plotSampleDistributions <-
function(pdat, title="", method=c("violin", "dist", "box"), x.text.size=7, n.grid.rows=3,
         hist.bins=100, x.text.angle=90, vmin=as.numeric(NA), vmax=as.numeric(NA), log.scale=TRUE,
         log.base=10, palette=cbPalette, fill=NULL, colour=NULL, hline=FALSE) {
  method <- match.arg(method)

  m <- reshape2::melt(pdat$tab, varnames=c("ID", "sample"))
  mt <- data.frame(pdat$metadata, row.names = pdat$metadata$sample)
  if(!is.null(fill)) m[['fill']] <- as.factor(mt[m$sample, fill])
  if(!is.null(colour)) m[['colour']] <- as.factor(mt[m$sample, colour])

  if(log.scale > 0) m$value <- log(m$value, base=log.base)
  lg <- ifelse(log.scale, paste("log", log.base), "")

  if(method == "box" | method == "violin") {
    g <- ggplot(m, aes(x=sample, y=value)) +
      simple_theme +
      ylim(vmin, vmax) +
      theme(axis.text.x = element_text(angle = x.text.angle, hjust = 1, vjust=0.5, size=x.text.size)) +
      labs(title=title, x="sample", y=paste0(lg, " value"))
    if(hline) g <- g + geom_hline(yintercept=0, colour='grey')
    if(method=="box") g <- g + geom_boxplot(outlier.shape=NA, na.rm=TRUE)
    if(method=="violin") g <- g + geom_violin(na.rm=TRUE, draw_quantiles=c(0.25, 0.5, 0.75), scale="width")
  } else if (method == "dist") {
    g <- ggplot(m, aes(x=value)) +
      geom_histogram(bins=hist.bins) +
      facet_wrap(~sample, nrow=n.grid.rows) +
      xlim(vmin, vmax)
  } else stop("Wrong method.")

  if(!is.null(fill)) {
    g <- g + aes(fill=fill)
    if(nlevels(m$fill) <= length(palette)) g <- g + scale_fill_manual(name=fill, values=palette)
  }
  if(!is.null(colour)) {
    g <- g + aes(colour=colour)
    if(nlevels(m$colour) <= length(palette)) g <- g + scale_color_manual(name=colour, values=palette)
  }
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
#' library(proteusUnlabelled)
#' data(proteusUnlabelled)
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
    stats <- rbind(stats, data.frame(id=rownames(w), condition=cond, mean=m, variance=v, ngood=ngood))
    rownames(stats) <- NULL
  }
  return(stats)
}


#' Create a table indicating good data.
#'
#' @param pdat Peptide or protein \code{proteusData} object.
#'
#' @return A data frame with rows corresponding to peptides/proteins and columns
#'   corresponding to conditions. Logical values indicate if at least one
#'   replicate is available (TRUE) or if all replicates are missing (FALSE).
#'
#' @examples
#' library(proteusUnlabelled)
#' data(proteusUnlabelled)
#' detect <- goodData(prodat)
#'
#' @export
goodData <- function(pdat) {
  if(!is(pdat, "proteusData")) stop ("Input data must be of class proteusData.")

  meta <- pdat$metadata
  G <- lapply(pdat$conditions, function(cond){
    w <- pdat$tab[,which(meta$condition == cond), drop=FALSE]
    rowSums(!is.na(w)) > 0
  })
  good <- as.data.frame(do.call(cbind, G))
  names(good) <- pdat$conditions
  return(good)
}


#' Plot mean-variance relationship
#'
#' \code{plotMV} makes a plot with variance of log-intensity vs mean of
#' log-intensity.
#'
#' @param pdat Peptide or protein \code{proteusData} object.
#' @param with.loess Logical. If true, a loess line will be added.
#' @param bins Number of bins for binhex
#' @param xmin Lower limit on x-axis
#' @param xmax Upper limit on x-axis
#' @param ymin Lower limit on y-axis
#' @param ymax Upper limit on y-axis
#' @param with.n Logical to add a text with the number of proteins to the plot
#' @param mid.gradient A parameter to control the midpoint of colour gradient
#'   (between 0 and 1)
#' @param text.size Text size on axes
#' @return A \code{ggplot} object.
#'
#' @examples
#' library(proteusUnlabelled)
#' data(proteusUnlabelled)
#' plotMV(prodat, with.loess=TRUE)
#'
#' @export
plotMV <- function(pdat, with.loess=FALSE, bins=80, xmin=5, xmax=10, ymin=7, ymax=20, with.n=FALSE, mid.gradient=0.3, text.size=10) {
  if(!is(pdat, "proteusData")) stop ("Input data must be of class proteusData.")
  meta <- pdat$metadata
  if(is.null(meta)) stop("No metadata found.")

  stats <- pdat$stats
  if(is.null(stats)) stats <- intensityStats(pdat)
  stats <- stats[which(!is.na(stats$mean) & !is.na(stats$variance) & stats$mean > 0 & stats$variance > 0),]
  stats$mean <- log10(stats$mean)
  stats$variance <- log10(stats$variance)
  protnum <- as.data.frame(table(stats$condition))  # number of proteins in each condition
  colnames(protnum) <- c("condition", "n")
  protnum$labn <- paste0("n = ", protnum$n)

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
    simple_theme +
    theme(panel.border = element_rect(fill=NA, color='black')) +
    xlim(xmin, xmax) +
    ylim(ymin, ymax) +
    facet_wrap(~condition) +
    stat_binhex(bins=bins) +
    theme(text = element_text(size=text.size)) +
    scale_fill_gradientn(colours=c("seagreen","yellow", "red"), values=c(0, mid.gradient, 1), name="count", na.value=NA)
  if(with.n) g <- g + geom_text(data=protnum, aes(x=xmin+0.5, y=ymax, label=labn))
  if(with.loess) g <- g + geom_line(data=ldf, aes(x,y), color='black')
  return(g)
}

#' Plot protein/peptide intensities
#'
#' \code{plotIntensities} makes a plot with peptide/protein intensity as a
#' function of the condition and replicate. When multiple IDs are entered, the
#' mean and standard error is plotted.
#'
#' @param pdat A \code{proteusData} object.
#' @param id Protein name, peptide sequence or a vector with these.
#' @param log Logical. If set TRUE a logarithm of intensity is plotted.
#' @param ymin Lower bound for y-axis
#' @param ymax Upper bound for y-axis
#' @param text.size Text size
#' @param point.size Point size
#' @param title Title of the plot (defaults to protein name)
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' library(proteusUnlabelled)
#' data(proteusUnlabelled)
#' prodat.med <- normalizeData(prodat)
#' plotIntensities(prodat.med, id='sp|P26263|PDC6_YEAST', log=TRUE)
#'
#' @export
plotIntensities <- function(pdat, id=NULL, log=FALSE, ymin=as.numeric(NA), ymax=as.numeric(NA),
                         text.size=12, point.size=3, title=NULL) {
  # without 'as.numeric' it returns logical NA (!!!)

  if(!is(pdat, "proteusData")) stop ("Input data must be of class proteusData.")
  if(is.null(id)) stop ("Need peptide/protein id.")

  if(pdat$content == "peptide") {
    selcol <- pdat$peptides
  } else if (pdat$content == "protein") {
    selcol <- pdat$proteins
  } else {
    stop("Unrecognized content in pdat object.")
  }
  sel <- which(selcol %in% id)

  if(length(sel) > 0 && sel > 0) {
    E <- if(log) log10(pdat$tab[sel,,drop=FALSE]) else pdat$tab[sel,,drop=FALSE]
    e <- colMeans(E, na.rm=TRUE)
    s <- sapply(E, function(x) sd(x, na.rm=TRUE)/sqrt(length(x)))
    n <- length(sel)

    if(is.null(title)) {
      title <- ifelse(n == 1, id, paste0("selection of ", n, " ", pdat$content, "s."))
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
    ggplot(p, aes(x=condition, y=expr, ymin=lo, ymax=up, fill=replicates, shape=shape)) +
      simple_theme_grid +
      theme(text = element_text(size=text.size), legend.position = "none") +
      ylim(ymin, ymax) +
      {if(n > 1) geom_errorbar(position=pd, width = 0.1)} +
      geom_point(position=pd, size=point.size) +
      scale_shape_identity() +  # necessary for shape mapping
      #scale_fill_manual(values=palette) +
      labs(x = 'Condition', y = ylab, title=title)
  }
}

#' Simple differential expression with limma
#'
#' \code{limmaDE} is a wrapper around \code{\link{limma}} to perform a
#' differential expression between a pair of conditions.
#'
#' @details
#'
#' Before \code{limma} is called, intensity data are transformed using the
#' \code{transform.fun} function. The default for this transformation is
#' \code{log10}. Therefore, by default, the column "logFC" in the output data
#' frame contains log10 fold change. If you need log2-based fold change, you can
#' use \code{transform.fun=log2}.
#'
#' \code{limmaDE} is only a simple wrapper around \code{\link{limma}}, to
#' perform differential expression between two conditions. For more complicated
#' designs we recommend using \code{\link{limma}} functions directly.
#'
#'
#' @param pdat Protein \code{proteusData} object.
#' @param formula A string with a formula for building the linear model.
#' @param conditions A character vector with two conditions for differential
#'   expression. Can be omitted if there are only two condition in \code{pdat}.
#' @param transform.fun A function to transform data before differential
#'   expression.
#' @param sig.level Significance level for rejecting the null hypothesis.
#' @return A data frame with DE results. "logFC" colum is a log-fold-change
#'   (using the \code{transform.fun}). Two columns wiht mean log-intensity
#'   (again, using \code{transform.fun}) are added. Attributes contain
#'   additional information about the transformation function, significance
#'   level, formula and conditions.
#'
#' @examples
#' library(proteusUnlabelled)
#' data(proteusUnlabelled)
#' prodat.med <- normalizeData(prodat)
#' res <- limmaDE(prodat.med)
#'
#' @export
limmaDE <- function(pdat, formula="~condition", conditions=NULL, transform.fun=log10, sig.level=0.05) {
  if(!is(pdat, "proteusData")) stop ("Input data must be of class proteusData.")

  meta <- pdat$metadata
  tab <- transform.fun(pdat$tab)

  # default conditions
  if(!is.null(conditions)) {
    for(cond in conditions ) {
      if(!(cond %in% meta$condition)) stop(paste("Condition", cond, "not found in metadata."))
    }
    sel <- which(meta$condition %in% conditions)
    meta <- droplevels(meta[sel,])
    tab <- tab[,sel]
  }

  if(nlevels(meta$condition) != 2) stop("This function requires exactly two conditions. Use the parameter conditions.")

  # limma analysis
  design <- model.matrix(as.formula(formula), meta)
  fit <- limma::lmFit(tab, design)
  ebay <- limma::eBayes(fit)
  coef <- colnames(ebay$design)[2]
  res <- limma::topTable(ebay, coef=coef, adjust="BH", sort.by="none", number=1e9)
  res <- cbind(rownames(res), res)
  colnames(res)[1] <- pdat$content
  res$significant <- res$adj.P.Val <= sig.level
  rownames(res) <- c()

  # add columns with mean intensity
  for(cond in levels(meta$condition)) {
    cname <- paste0("mean_", cond)
    w <- transform.fun(pdat$tab[,which(meta$condition == cond), drop=FALSE])
    m <- rowMeans(w, na.rm=TRUE)
    m[which(is.nan(m))] <- NA
    res[, cname] <- m
  }

  attr(res, "transform.fun") <- deparse(substitute(transform.fun))
  attr(res, "sig.level") <- sig.level
  attr(res, "formula") <- formula
  attr(res, "conditions") <- paste(levels(meta$condition), collapse=",")

  return(res)
}


#' Fold-change intensity diagram
#'
#' \code{plotFID} makes a log10 fold change versus log10 sum intensity plot,
#' usually known as MA plot.
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
#' @return A \code{ggplot} object.
#'
#' @examples
#' library(proteusUnlabelled)
#' data(proteusUnlabelled)
#' prodat.med <- normalizeData(prodat)
#' plotFID(prodat.med)
#'
#' @export
plotFID <- function(pdat, pair=NULL, pvalue=NULL, bins=80, marginal.histograms=FALSE,
                   xmin=NULL, xmax=NULL, ymax=NULL, text.size=12, show.legend=TRUE, plot.grid=TRUE,
                   binhex=TRUE) {
  if(binhex & marginal.histograms) {
    warning("Cannot plot with both binhex=TRUE and marginal.histograms=TRUE. Ignoring binhex.")
    binhex=FALSE
  }

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
    {if(binhex) stat_binhex(bins=bins, show.legend=show.legend) else geom_point()}+
    scale_fill_gradientn(colours=c("seagreen","yellow", "red"), name = "count",na.value=NA) +
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
#' @param res Output table from \code{\link{limmaDE}}.
#' @param text.size Text size.
#' @param plot.grid Logical to plot grid.
#' @param bin.size Bin size for the histogram.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' library(proteusUnlabelled)
#' data(proteusUnlabelled)
#' prodat.med <- normalizeData(prodat)
#' res <- limmaDE(prodat.med)
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
#'@param res Result table from  \code{\link{limmaDE}}.
#'@param bins Number of bins for binhex.
#'@param xmax Upper limit on x-axis. If used, the lower limit is -xmax.
#'@param ymax Upper limit on y-axis. If used, the lower limit is -ymax.
#'@param marginal.histograms A logical to add marginal histograms.
#'@param text.size Text size.
#'@param show.legend Logical to show legend (colour key).
#'@param plot.grid Logical to plot grid.
#'@param binhex Logical. If TRUE, a hexagonal densit plot is made, otherwise it
#'  is a simple point plot.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' library(proteusUnlabelled)
#' data(proteusUnlabelled)
#' prodat.med <- normalizeData(prodat)
#' res <- limmaDE(prodat.med)
#' plotVolcano(res)
#'
#'@export
plotVolcano <- function(res, bins=80, xmax=NULL, ymax=NULL, marginal.histograms=FALSE, text.size=12, show.legend=TRUE,
                        plot.grid=TRUE, binhex=TRUE) {
  if(binhex & marginal.histograms) {
    warning("Cannot plot with both binhex=TRUE and marginal.histograms=TRUE. Ignoring binhex.")
    binhex=FALSE
  }

  tr <- attr(res, "transform.fun")
  xlab <- ifelse(is.null(tr), "FC", paste(tr, "FC"))
  id <- names(res)[1]

  g <- ggplot(res, aes(logFC, -log10(P.Value))) +
    {if(plot.grid) simple_theme_grid else simple_theme} +
    {if(binhex) stat_binhex(bins=bins, show.legend=show.legend) else geom_point(aes_string(text=id))} +
    scale_fill_gradientn(colours=c("seagreen","yellow", "red"), name = "count", na.value=NA) +
    geom_vline(colour='red', xintercept=0) +
    theme(text = element_text(size=text.size)) +
    labs(x=xlab, y="-log10 P")

    if(!is.null(xmax)) g <- g + scale_x_continuous(limits = c(-xmax, xmax), expand = c(0, 0))
    if(!is.null(ymax) ) g <- g + scale_y_continuous(limits = c(0, ymax), expand = c(0, 0))

    if(marginal.histograms) g <- ggExtra::ggMarginal(g, size=10, type = "histogram", xparams=list(bins=100), yparams=list(bins=50))
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
#' @param palette Palette of colours
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' library(proteusUnlabelled)
#' data(proteusUnlabelled)
#' prodat.med <- normalizeData(prodat)
#' plotProtPeptides(pepdat.clean, 'sp|P26263|PDC6_YEAST', prodat.med)
#'
#' @export
plotProtPeptides <- function(pepdat, protein, prodat=NULL, palette=cbPalette) {
  if(!is(pepdat, "proteusData")) stop ("Input data must be of class proteusData.")
  tab <- normalizeMedian(pepdat$tab)

  # all peptides for this protein
  selprot <- which(pepdat$pep2prot$protein == protein)
  if(length(selprot) == 0) stop(paste0("Protein '", protein, "' not found."))

  # need as.character, or indexing by factor is wrong!
  peps <- as.character(pepdat$pep2prot[selprot,'sequence'])
  mat <- as.matrix(tab[peps,])
  dat <- reshape2::melt(mat, varnames=c("peptide", "sample"))
  dat$sample <- factor(dat$sample, levels=pepdat$metadata$sample)
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
    scale_fill_manual(values=palette) +
    geom_boxplot(outlier.shape = NA)  +
    geom_jitter(width=0, size=0.5) +
    facet_wrap(~condition) +
    theme(legend.position="none") +
    labs(x="Peptide", y="log intensity", title=protein)
  g2 <- ggplot(dat, aes(x=sample, y=intensity, fill=condition)) +
    scale_fill_manual(values=palette) +
    geom_boxplot(outlier.shape = NA)  +
    geom_jitter(width=0, size=0.5) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
    theme(legend.position="none") +
    labs(x="Sample", y="log intensity")
  if(!is.null(prodat)) g2 <- g2 + geom_point(aes(x=sample, y=prot.intensity), shape=22, size=3, fill='white')
  gridExtra::grid.arrange(g1, g2, ncol=1)
}

