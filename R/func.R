#' Evidence columns
#'
#' Columns to be read from the evidence file. This matrix contains two rows: first with the name used
#' in the package and second with the original name.
evidenceColumns <- matrix(c(
  'sequence', 'Sequence',
  'modseq', 'Modified sequence',
  'modifications', 'Modifications',
  'proteins', 'Proteins',
  'protein', 'Leading razor protein',
  'type', 'Type',
  'experiment', 'Experiment',
  'charge', 'Charge',
  'numpoints', 'Number of data points',
  'numscans', 'Number of scans',
  'pep', 'PEP',
  'score', 'Score',
  'intensity', 'Intensity',
  'reverse', 'Reverse',
  'contaminant', 'Potential contaminant'
),2
)

#' Read MaxQuant's evidence file
#'
#' Works only with unlabelled data, that is, the evidence file with one Intensity column.
#' Contaminants and reverse sequences are filtered out.
#'
#' @param file File name.
#' @return Data frame with selected columns from the evidence file.
readEvidenceFile <- function(file) {
  evi <- read.delim(file, header=TRUE, sep="\t", check.names=FALSE, as.is=TRUE, strip.white=TRUE)
  evi <- evi[,evidenceColumns[2,]]
  names(evi) <- evidenceColumns[1,]
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
  return(norm)
}

######

#' Normalize peptide table
#'
#' Normalize peptide table with normalizing factors. Column and row names are retained.
#' @param peptab Peptab table.
#' @param normfac A vector of normalizing factors.
normalizeTable <- function(peptab, normfac) {
  rows <- rownames(peptab)
  cols <- colnames(peptab)
  peptab <- t(t(peptab) / normfac)
  rownames(peptab) <- rows
  colnames(peptab) <- cols
  return(peptab)
}

######

#' Create peptide table from evidence data
#'
#' Creates a table with columns corresponding to samples (experiments) and rows corresponding to peptides.
#' Each cell is a sum of all intensities for this sample/peptide in the input evidence data.
#' @param evi Evidence table created with readEvidenceFile.
#' @return Peptide intensity table.
makePeptideTable <- function(evi) {
  peptab.raw <- cast(evi, sequence ~ experiment, sum, value = 'intensity')
  peptab.raw[peptab.raw == 0] <- NA
  rownames(peptab.raw) <- peptab.raw$sequence
  peptab.raw <- peptab.raw[,2:ncol(peptab.raw)]
  # cast creates a 'cast_data_frame' or something and this screws things up
  peptab.raw <- as.data.frame(peptab.raw)
  return(peptab.raw)
}

#######

#' Split peptide table into conditions
#'
#' Split peptide table (containing all samples) into separate tables, one per condition.
#' Condition information is provided in the meta table.
#' @param meta Metadata data frame, needs columns 'sample' and 'condition'.
#' @params peptab Peptide table
#' @return A list of intensity tables, one per condition.
splitConditions <- function(peptab, meta) {
  conditions <- levels(meta$condition)
  int <- list()
  for(condition in conditions) {
    colnames <- as.character(meta$sample[which(meta$condition == condition)])
    tab <- peptab[,colnames]
    int[[condition]] <- tab
  }
  return(int)
}

#######

#' Normality test wrapper
#'
#' Performs a normality test for the intensity table list (all conditions).
#' @param int Intensity table list (created with splitConditions)
#' @return A list of data frames with mean, p-value and adjusted p-value
testNormalityConditions <- function(int) {
  norm.test <- list()
  for(condition in names(int)) {
    w <- int[[condition]]
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
#' Calculate Z-scores for a vectore, NAs are ignored.
#' @param x Input numeric vector
#' @return Z-score vector
standardize <- function(x){
  x <- as.numeric(x)
  s <- (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
  return(s)
}

#####
#' Plot peptide count per sample
#'
#' @param peptab Peptide table
#' @meta meta Metadata table (needed to separate conditions)
#' @return A plot of the number of peptides detected in each sample
plotPeptideCount <- function(peptab, meta){
  pep.count <- sapply(peptab, function(x) sum(!is.na(x)))
  df <- data.frame(x=samples, y=pep.count, condition=meta$condition)
  ggplot(df, aes(x=x,y=y,color=condition)) + geom_col() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=7)) + labs(x='Sample', y='Peptide count')
}

#' Detect downliers
#'
#' Detect a downlier based on Z-score.
#' @param v An input numeric vector
#' @param sigma Z-score limit, default value is 5
#' @return Index of the downlier if detected (that is the lowest intensity Z-score is less than sigma), otherwise zero.
downlierSigma <- function(v, sigma=5) {
  v <- na.omit(v)
  n <- length(v)
  if(n < 7) return(0)

  Z <- standardize(v)
  down <- which.min(Z)

  if(Z[down] <= -sigma) {
    return(down)
  } else {
    return(0)
  }
}
