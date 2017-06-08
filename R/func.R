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

library(ggplot2)
simple_theme <- theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_line(colour = "grey90"),
    panel.grid.minor = element_line(colour = "grey95"),
    axis.line = element_line(colour = "black")
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
#' @param tab Intensity table for all conditions.
#' @param method Method of normalization, only 'median' at the moment.
#' @return A vector with normalizing factors.
normalizingFactors <- function(tab, method="median") {
  if(method == 'median') {
    norm <- apply(tab, 2, function(x) {median(x, na.rm=TRUE)})
    norm <- norm / mean(norm)
  } else {
    stop("Unknown normalization method.")
  }
  attr(norm, "method") <- method
  return(norm)
}

######

#' Normalize intensity table
#'
#' Normalize intensity table with normalizing factors. Attributes are retained.
#' @param tab Peptab table.
#' @param normfac A vector of normalizing factors.
normalizeTable <- function(tab, normfac) {
  atr <- attributes(tab)
  tab <- data.frame(t(t(tab) / normfac))  # this clears all attributes!
  attributes(tab) <- atr
  return(tab)
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
#' @return An intensity table structure
makePeptideTable <- function(evi, meta, id="sequence", value="intensity") {
  if(!id %in% c("sequence", "modseq")) stop("Incorrect id. Has to be 'sequence' or 'modseq'.")
  form <- as.formula(paste0(id, " ~ experiment"))
  tab <- cast(evi, form, sum, value = value)
  tab[tab == 0] <- NA
  rownames(tab) <- tab[,1]
  tab <- tab[,2:ncol(tab)]
  # cast creates a 'cast_data_frame' or something and this screws things up; need df
  peptab <- list(
    tab = as.data.frame(tab),
    metadata = meta,
    norm = "none",
    logflag = FALSE,
    id = id,
    value = value
  )
  return(peptab)
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
#' @param dat Intensity data, unnormalized raw data
#' @param norm How to normalize intensities (default: "median")
#' @return A data frame with test results
testNormalityConditions <- function(dat, norm="median") {
  if(dat$norm == "none") {
    normfac <- normalizingFactors(dat$tab, norm)
    tab <- normalizeTable(dat$tab, normfac)
  } else {
    tab <- dat$tab
  }
  if(!dat$logflag) tab <- log10(tab)

  conditions <- dat$meta$condition
  norm.test <- NULL
  for(cond in levels(conditions)) {
    w <- tab[,which(conditions == cond)]
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
    norm.test <- rbind(norm.test, data.frame(
      condition = cond,
      mean = rowMeans(w, na.rm=TRUE),
      p = P,
      p.adj = p.adjust(P, method="BH")
    ))
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
  med.count <- median(pep.count)
  df <- data.frame(x=meta$sample, y=pep.count, condition=meta$condition)
  ggplot(df, aes(x=x,y=y,color=condition)) +
    geom_col() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=7)) +
    labs(x='Sample', y='Peptide count') +
    labs(title = paste0("Median peptide count = ", med.count)) +
    theme(plot.title=element_text(hjust=0, size=10))
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
#' @param pdat Intensity dat tructure
#' @return A data frame with several statistics
intensityStats <- function(pdat) {
  meta <- pdat$metadata
  logflag <- pdat$logflag
  if(is.null(logflag)) logflag <- FALSE

  conditions <- levels(meta$condition)
  stats <- NULL
  for(cond in conditions) {
    w <- pdat$tab[,which(meta$condition == cond)]
    if(logflag) w <- 10^w
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
#' Plot variance of log-intensity vs mean of log-intensity.
#' @param pdat Intensity data structure.
#' @param with.loess Logical. If true, a loess line will be added.
plotMV <- function(pdat, with.loess=FALSE, xmin=5, xmax=10, ymin=7, ymax=20) {
  meta <- pdat$metadata
  if(is.null(meta)) stop("No metadata found.")
  conditions <- meta$condition

  stats <- pdat$stats
  if(is.null(stats)) stats <- intensityStats(pdat)
  stats$mean <- log10(stats$mean)
  stats$variance <- log10(stats$variance)

  # has to be calculated for each condition separately
  if(with.loess) {
    ldf <- NULL
    for(cond in levels(conditions))
    {
      st <- stats[which(stats$condition == cond),]
      ls <- loess(variance ~ mean, data=st)
      x <- seq(from=min(na.omit(st$mean)), to=max(na.omit(st$mean)), by=0.01)
      pr <- predict(ls, x)
      ldf <- rbind(ldf, data.frame(condition=cond, x=x, y=pr))
    }
  }

  ggplot(s, aes(x=log10(mean), y=log10(variance))) +
    simple_theme +
    xlim(xmin, xmax) +
    ylim(ymin, ymax) +
    facet_wrap(~condition) +
    stat_binhex(bins=80) +
    scale_fill_gradientn(colours=c("green","yellow", "red"), name = "count",na.value=NA) +
    if(with.loess) {geom_line(data=ldf, aes(x,y), color='black')}

}


#' Make protein table
#'
#' Create a protein table from the peptide table using high-flyers.
#' @param evi Evidence table, used only to cross-reference peptides with proteins.
#' @param peptab Peptide intensity table.
#' @param hifly The number of top peptides (high-flyers) to be used for protein intensity.
#' @param norm Normalization to use on peptides before converting into proteins.
#' @param min.peptides Minimum number of peptides per protein.
#' @return Protein intensity table.
makeProteinTable <- function(evi, peptab, hifly=3, norm="median", min.peptides=1, verbose=FALSE) {
  meta <- peptab$metadata
  normfac <- normalizingFactors(peptab$tab, method=norm)
  tab <- normalizeTable(peptab$tab, normfac)

  pep2prot <- select(evi, sequence, protein)
  pep2prot <- unique(pep2prot)

  # make sure order of peptides is as in intensity tables
  peptides <- data.frame(sequence=rownames(tab))
  peptides <- merge(peptides, pep2prot, by="sequence")
  proteins <- levels(as.factor(peptides$protein))

  protlist <- list()
  for(cond in levels(meta$condition)) {
    w <- tab[,which(meta$condition == cond)]
    samples <- colnames(w)
    protint <- NULL
    for(prot in proteins) {
      sel <- which(peptides$protein == prot)
      if(length(sel) >= min.peptides)
      {
        # ARGHHH! Took me forever to get it right
        # without drop=FALSE it will drop a dimension for one-row selection
        # and all downstream analysis goes to hell
        wp <- w[sel,, drop=FALSE]
        medPep <- apply(wp, 1, function(x) {median(x, na.rm=TRUE)})  # median across replicates
        nmed <- length(medPep)
        top <- min(c(nmed, hifly))
        itop <- sort.int(medPep, decreasing=TRUE, index.return=TRUE)$ix[1:top] #index of top peptides
        row <- data.frame(protein=prot, t(colMeans(wp[itop,]))) # mean of the top peptides
        protint <- rbind(protint, row)
      }
    }
    colnames(protint) <- c("protein", samples)
    protlist[[cond]] <- protint
  }

  # dplyr join all tables
  protlist %>% Reduce(function(df1, df2) full_join(df1,df2, by="protein"), .) -> protab
  rownames(protab) <- protab$protein
  protab <- protab[,2:ncol(protab)]
  prodat <- list(
    tab = protab,
    metadata = meta,
    pep2prot = pep2prot,
    norm = norm,
    logflag = FALSE
  )
  prodat$stats <- intensityStats(prodat)

  return(prodat)
}

