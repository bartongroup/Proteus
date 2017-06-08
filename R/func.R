#' Evidence columns
#'
#' Columns to be read from the evidence file. This list contains names of columns to be read.
#' The names of list elements are used internally to reference evidence data.
evidenceColumns <- list(
  sequence = 'Sequence',
  modseq = 'Modified sequence',
  modifications = 'Modifications',
  proteins = 'Proteins',
  protein = 'Leading razor protein',
  type = 'Type',
  experiment = 'Experiment',
  charge = 'Charge',
  numpoints = 'Number of data points',
  numscans = 'Number of scans',
  pep = 'PEP',
  score = 'Score',
  intensity = 'Intensity',
  reverse = 'Reverse',
  contaminant = 'Potential contaminant'
)

#' Read MaxQuant's evidence file
#'
#' Works only with unlabelled data, that is, the evidence file with one Intensity column.
#' Contaminants and reverse sequences are filtered out.
#'
#' @param file File name.
#' @param columns Named list with columns to read.
#' @return Data frame with selected columns from the evidence file.
readEvidenceFile <- function(file, columns=evidenceColumns) {
  evi <- read.delim(file, header=TRUE, sep="\t", check.names=FALSE, as.is=TRUE, strip.white=TRUE)
  evi <- evi[, as.character(columns)]
  names(evi) <- names(columns)
  # sometimes there are only NAs and the condition doesn't work
  evi$reverse[is.na(evi$reverse)] = ''
  evi$contaminant[is.na(evi$contaminant)] = ''
  evi <- evi[which(evi$contaminant != '+' & evi$reverse != '+'),]
  evi <- evi[which(!is.na(evi$intensity)),]
}

#' Read MaxQuant's peptide file
#'
#' @param file File name.
#' @return Data frame with selected columns from the evidence file.
readPeptideFile <- function(file) {
  pep <- read.delim(file, header=TRUE, sep="\t", check.names=FALSE, as.is=TRUE, strip.white=TRUE)
  pep$Reverse[is.na(pep$Reverse)] = ''
  pep$`Potential contaminant`[is.na(pep$`Potential contaminant`)] = ''
  pep <- pep[which(pep$`Potential contaminant` != '+' & pep$Reverse != '+'),]
  attr(pep, "norm") <- "none"
}


#' Calculate normalizing factors
#'
#' @param peptab peptide intensity table for all conditions.
#' @param method method of normalization, only 'median' at the moment.
#' @return A vector with normalizing factors.
normalizingFactors <- function(peptab, method="median") {
  if(method == 'median') {
    norm <- apply(peptab, 2, function(x) {median(x, na.rm=TRUE)})
    norm <- norm / mean(norm)
  } else {
    stop("Unknown normalization method.")
  }
  attr(norm, "method") <- method
  return(norm)
}

######

#' Normalize peptide table
#'
#' Normalize peptide table with normalizing factors. Attributes are retained.
#' @param peptab Peptab table.
#' @param normfac A vector of normalizing factors.
normalizeTable <- function(peptab, normfac) {
  atr <- attributes(peptab)
  peptab <- data.frame(t(t(peptab) / normfac))  # this clears all attributes!
  attributes(peptab) <- atr
  attr(peptab, "norm") = attr(normfac, "method")
  return(peptab)
}

######

#' Create peptide table from evidence data
#'
#' Creates a table with columns corresponding to samples (experiments) and rows corresponding to peptides.
#' Each cell is a sum of all intensities for this sample/peptide in the input evidence data.
#' @param evi Evidence table created with readEvidenceFile.
#' @param meta Data frame with metadata. As minimum, it should contain "sample" and "condition".
#' @param id A column name to identify peptides. The default is "sequence". Can be "modseq".
#' @param value A column name to use for results. The default is "intensity".
#' @return Peptide intensity table.
makePeptideTable <- function(evi, meta, id="sequence", value="intensity") {
  if(!id %in% c("sequence", "modseq")) stop("Incorrect id. Has to be 'sequence' or 'modseq'.")
  form <- as.formula(paste0(id, " ~ experiment"))
  peptab.raw <- cast(evi, form, sum, value = value)
  peptab.raw[peptab.raw == 0] <- NA
  rownames(peptab.raw) <- peptab.raw[,1]
  peptab.raw <- peptab.raw[,2:ncol(peptab.raw)]
  # cast creates a 'cast_data_frame' or something and this screws things up
  peptab.raw <- as.data.frame(peptab.raw)
  attr(peptab.raw, "metadata") <- meta
  attr(peptab.raw, "norm") <- "none"
  attr(peptab.raw, "logflag") <- FALSE
  attr(peptab.raw, "id") <- id
  attr(peptab.raw, "value") <- value
  return(peptab.raw)
}

#######

#' Split peptide table into conditions
#'
#' Split peptide table (containing all samples) into separate tables, one per condition.
#' Condition information is provided in the meta table.
#' @param intab Intensity table
#' @return A list of intensity tables, one per condition.
splitConditions <- function(intab) {
  meta <- attr(intab, "metadata")
  conditions <- levels(meta$condition)
  inlist <- list()
  for(condition in conditions) {
    colnames <- as.character(meta$sample[which(meta$condition == condition)])
    tab <- intab[,colnames]
    inlist[[condition]] <- tab
  }
  return(inlist)
}

#######

#' Normality test wrapper
#'
#' Performs a normality test for the intensity table list (all conditions).
#' @param intab Intensity table, unnormalized raw data
#' @param norm How to normalize intensities (default: "median")
#' @return A list of data frames with mean, p-value and adjusted p-value
testNormalityConditions <- function(intab, norm="median") {
  normfac <- normalizingFactors(intab, norm)
  intab <- normalizeTable(intab, normfac)
  intab <- log10(intab)
  inlist <- splitConditions(intab)
  norm.test <- list()
  for(condition in names(inlist)) {
    w <- inlist[[condition]]
    P <- c()
    # struggled with apply here, so just an old-fashined loop
    for(i in 1:nrow(w)) {
      p <- 1
      x <- na.omit(as.numeric(w[i,]))
      if(length(x) > 7) {
        #sw <- shapiro.test(x)
        ad <- ad.test(x)
        p <- ad$p.value
      }
      P <- c(P, p)
    }
    norm.test[[condition]] <- data.frame(
      mean = rowMeans(w, na.rm=TRUE),
      p = P,
      p.adj = p.adjust(P, method="BH")
    )
  }
  return(norm.test)
}

#####

#' Anderson-Darling test for normality
#'
#' Performs a normality test for the intensity table, row by row.
#' @param w Intensity table (columns: replicates, rows: peptides/proteins/genes etc.)
#' @return Data frame with mean intensity, p-value and BH-adjusted p-value
testNormality <- function(w) {
  P <- c()
  # struggled with apply here, so just an old-fashined loop
  for(i in 1:nrow(w)) {
    p <- 1
    x <- na.omit(as.numeric(w[i,]))
    if(length(x) > 7) {
      ad <- ad.test(x)
      p <- ad$p.value
    }
    P <- c(P, p)
  }
  data.frame(
    mean = rowMeans(w, na.rm=TRUE),
    p = P,
    p.adj = p.adjust(P, method="BH")
  )
}


#####

#' Standardize a vector
#'
#' Calculate Z-scores for a vector, NAs are ignored.
#' @param v Input numeric vector
#' @param trim A number of points to be trimmed on each side for mean and sd
#' @return Z-score vector
standardize <- function(v, trim=0){
  v <- as.numeric(v)
  x <- v[!is.na(v)]
  n <- length(x)
  #if(n < 2) stop("Need at lest 2 points in standardize")

  if(trim > 0) {
   if(n - 2*trim < 2) stop("Too much trim")
    x <- sort(x)
    lo <- trim + 1
    up <- n - trim
    x <- x[lo:up]
  }
  s <- (v - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
  return(s)
}

#' Trim a vector
#'
#' Trim a vector by replacing "trim" top and bottom elements with NA.
#'
#' @param v Input vector
#' @param trim number of elements to remove on each side
#' @return Trimmed vector
trimVector <- function(v, trim) {
  if(trim == 0) return(v)
  n <- length(which(!is.na(v)))
  if(n - 2*trim < 2) stop("Too much trim")
  # index of sorted elements. na.last indicates to put NAs at the end. This works only with radix.
  idx <- sort.int(x, index.return = TRUE, na.last=TRUE, method="radix")$ix
  lo <- trim
  up <- n-trim+1
  last <- length(idx)
  v[idx[1:lo]] <- NA
  v[idx[up:last]] <- NA
  return(v)
}

#####
#' Plot peptide count per sample
#'
#' @param peptab Peptide table (with metadata attribute)
#' @return A plot of the number of peptides detected in each sample
plotPeptideCount <- function(peptab){
  meta <- attr(peptab, "metadata")
  if(is.null(meta)) stop("Attribute 'metadata' is missing from table")
  pep.count <- sapply(peptab, function(x) sum(!is.na(x)))
  df <- data.frame(x=meta$sample, y=pep.count, condition=meta$condition)
  ggplot(df, aes(x=x,y=y,color=condition)) + geom_col() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=7)) + labs(x='Sample', y='Peptide count')
}

#' Detect downliers
#'
#' Detect a downlier based on Z-score.
#' @param v An input numeric vector
#' @param sigma Z-score limit, default value is 5
#' @param trim Number of trimmed points on each side for finding Z-score
#' @return Integer vector with index of all downliers. If there are no downliers, the vector is empty.
downlierSigma <- function(v, sigma=5, trim=0) {
  n <- length(which(!is.na(v)))
  if(n < 6) return(integer(0))

  Z <- standardize(v, trim)
  down <- which(Z <= -sigma)
  return(down)
}

#' Statistics for an intensity table
#'
#' Calculate mean, standard deviation, ... for an intensity table (that is peptide or protein,
#' for all conditions)
#' @param intab Intensity table
#' @return List of data frames with statistics, one per conditions
intensityStats <- function(intab) {
  meta <- attr(intab, "metadata")
  logflag <- attr(intab, "logflag")
  if(is.null(logflag)) logflag <- FALSE

  intlist <- splitConditions(intab)
  stats <- list()
  for(condition in conditions) {
    w <- intlist[[condition]]
    if(!logflag) w <- log10(w)
    m <- rowMeans(w, na.rm=TRUE)
    v <- apply(w, 1, function(v) sd(v, na.rm=TRUE)^2)
    stats[[condition]] <- data.frame(mean=m, variance=v)
  }
  return(stats)
}


#' Plot mean-variance relationship
#'
#' Plot variance of log-intensity vs mean of log-intensity.
#' @param intab Intensity table, preferably with added stats.
#' @param with.loess Logical. If true, a loess line will be added.
plotMV <- function(intab, with.loess=FALSE, xmin=5, xmax=10, ymin=7, ymax=20) {
  meta <- attr(intab, "metadata")
  if(is.null(meta)) stop("No metadata found.")
  conditions <- meta$condition

  stats <- attr(intab, "stats")
  if(is.null(stats)) stats <- intensityStats(intab)

  P <- list()
  for(condition in conditions) {
    st <- stats[[condition]]
    df <- data.frame(mean=log10(st$mean), variance=log10(st$variance))
    P[[condition]] <- ggplot(df, aes(mean, variance)) + geom_point(alpha=0.3) +
      labs(x="log mean", y="log variance", title=condition) +
      xlim(xmin, xmax) +
      ylim(ymin, ymax)
    if(with.loess) {
      ls <- loess(variance ~ mean, data=df)
      x <- seq(from=min(na.omit(df$mean)), to=max(na.omit(df$mean)), by=0.05)
      pr <- predict(ls, x)
      pf <- data.frame(x=x, y=pr)
      P[[condition]] <- P[[condition]] + geom_line(data=pf, aes(x,y), color='yellow')
    }
  }
  grid.arrange(grobs=P, ncol=length(P))
}
