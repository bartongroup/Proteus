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
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )
simple_theme_grid <- theme_bw() +
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
  return(norm)
}


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


#' Create peptide table from evidence data
#'
#' Creates a table with columns corresponding to samples (experiments) and rows corresponding to peptides.
#' Each cell is a sum of all intensities for this sample/peptide in the input evidence data.
#' @param evi Evidence table created with readEvidenceFile.
#' @param meta Data frame with metadata. As minimum, it should contain "sample" and "condition".
#' @param pepseq A column name to identify peptides. Can be "sequence" or "modseq".
#' @param intensity A column name to use for results. The default is "intensity".
#' @return An intensity table structure
makePeptideTable <- function(evi, meta, pepseq="sequence", intensity="intensity") {

  if(!(pepseq %in% c("sequence", "modseq"))) stop("Incorrect pepseq. Has to be 'sequence' or 'modseq'.")

  # cast evidence data (long format) into peptide table (wide format)
  form <- as.formula(paste0(pepseq, " ~ experiment"))
  tab <- cast(evi, form, sum, value = intensity)
  tab[tab == 0] <- NA
  rownames(tab) <- tab[,1]
  tab <- tab[,2:ncol(tab)]
  # cast creates a 'cast_data_frame' or something and this screws things up; need data frame
  tab <- as.data.frame(tab)
  # keep only columns from the metadata file
  # you can remove bad replicates
  tab <- tab[,meta$sample]

  # peptide to protein conversion
  peptides <- rownames(tab)
  pep2prot <- select(evi, sequence, protein)
  pep2prot <- unique(pep2prot)
  rownames(pep2prot) <- pep2prot$sequence
  pep2prot <- pep2prot[peptides,]
  proteins <- levels(as.factor(pep2prot$protein))

  # create pdat object
  pdat <- list(
    tab = tab,
    content = 'peptide',
    type = 'unlabelled',
    metadata = meta,
    norm = "none",
    logflag = FALSE,
    pepseq = pepseq,
    intensity = intensity,
    pep2prot = pep2prot,
    peptides = peptides,
    proteins = proteins
  )
  return(pdat)
}


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


#' Plot correlation matrix
#'
#' Plot correlation matrix for peptide or protein data.
#' @param pdat Peptide/protein data object
plotCorrelationMatrix <- function(pdat) {
  corr.mat <- cor(pdat$tab, use="complete.obs")
  heatmap.2(corr.mat, trace="none", density.info="none", dendrogram="none", Rowv=FALSE, Colv=FALSE, key.xlab = "Correlation coefficient")
}

#' Plot clustering dendrogram
#'
#' @param pdat Peptide/protein data object.
plotClustering <- function(pdat) {
  corr.mat <- cor(pdat$tab, use="complete.obs")
  dis <- as.dist(1 - corr.mat)  # dissimilarity matrix
  dend <- as.dendrogram(hclust(dis))
  colors_to_use <- as.numeric(pdat$metadata$condition)
  colors_to_use <- colors_to_use[order.dendrogram(dend)]
  labels_colors(dend) <- colors_to_use
  plot(dend)
}


#' Plot peptide count per sample
#'
#' @param pdat Peptide data object
#' @param x.text.size Size of text on the x-axis
#' @return A plot of the number of peptides detected in each sample
plotPeptideCount <- function(pdat, x.text.size=7){
  meta <- pdat$meta
  pep.count <- sapply(pdat$tab, function(x) sum(!is.na(x)))
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
#' @param bins Number of bins for binhex
#' @param xmin Lower limit on x-axis
#' @param xmax Upper limit on x-axis
#' @param ymin Lower limit on y-axis
#' @param ymax Upper limit on y-axis
plotMV <- function(pdat, with.loess=FALSE, bins=80, xmin=5, xmax=10, ymin=7, ymax=20) {
  meta <- pdat$metadata
  if(is.null(meta)) stop("No metadata found.")
  conditions <- meta$condition

  stats <- pdat$stats
  if(is.null(stats)) stats <- intensityStats(pdat)
  stats <- stats[which(!is.na(stats$mean) & !is.na(stats$variance)),]
  stats$mean <- log10(stats$mean)
  stats$variance <- log10(stats$variance)
  protnum <- as.data.frame(table(stats$condition))  #number of proteins in each condition
  colnames(protnum) <- c("condition", "n")

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

  ggplot(stats, aes(x=mean, y=variance)) +
    simple_theme_grid +
    xlim(xmin, xmax) +
    ylim(ymin, ymax) +
    facet_wrap(~condition) +
    stat_binhex(bins=bins) +
    scale_fill_gradientn(colours=c("green","yellow", "red"), name = "count",na.value=NA) +
    geom_text(data=protnum, aes(x=xmin+0.5, y=ymax, label=paste0("n = ", n)))
   # if(with.loess) {geom_line(data=ldf, aes(x,y), color='black')}

}


#' Make protein table
#'
#' Create a protein table from the peptide table using high-flyers.
#' @param pepdat Peptide intensity object.
#' @param hifly The number of top peptides (high-flyers) to be used for protein intensity.
#' @param norm Normalization to use on peptides before converting into proteins.
#' @param min.peptides Minimum number of peptides per protein.
#' @return Protein intensity object.
makeProteinTable <- function(pepdat, hifly=3, norm="median", min.peptides=1, verbose=FALSE) {
  meta <- pepdat$metadata
  normfac <- normalizingFactors(pepdat$tab, method=norm)
  tab <- normalizeTable(pepdat$tab, normfac)

  protlist <- list()
  for(cond in levels(meta$condition)) {
    w <- tab[,which(meta$condition == cond)]
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
        if(npep > 1) {
          wp[is.na(wp)] <- 0   # zeroes for sorting only
          sp <- as.matrix(apply(wp, 2, sort, decreasing=TRUE))
          sp[sp == 0] <- NA  # put NAs back
          ntop <- min(npep, hifly)
          row <- t(colMeans(sp[1:ntop,,drop=FALSE], na.rm=TRUE))
          row[is.nan(row)] <- NA    # colMeans puts NaNs where the column contains only NAs
        } else {
          row <- wp
        }
        row <- data.frame(protein=prot, row)

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
    colnames(protint) <- c("protein", samples)
    protlist[[cond]] <- protint
  }

  # dplyr join all tables
  protlist %>% Reduce(function(df1, df2) full_join(df1,df2, by="protein"), .) -> protab
  rownames(protab) <- protab$protein
  protab <- protab[,2:ncol(protab)]
  prodat <- list(
    tab = protab,
    content = "protein",
    type = pepdat$type,
    pepseq = pepdat$pepseq,
    intensity = pepdat$intensity,
    metadata = meta,
    hifly = hifly,
    mmin.peptides = min.peptides,
    norm = norm,
    logflag = FALSE,
    pep2prot = pepdat$pep2prot,
    peptides = pepdat$peptides,
    proteins = pepdat$proteins
  )
  prodat$stats <- intensityStats(prodat)

  return(prodat)
}


#' Plot protein(s)
#'
#' Plot protein intensity as a function of the condition and replicate. When
#' multiple proteins are entered, the mean and standard error is plotted.
#' @param pdat Protein intensity structure.
#' @param protein Protein name (string) or a vector with protein names.
#' @param log Logical. If set TRUE a logarithm of intensity is plotted.
#' @param ymin Lower bound for y-axis
#' @param ymax Upper bound for y-axis
#' @param text.size Text size
#' @param point.size Point size
#' @param title Title of the plot (defaults to protein name)
# without 'as.numeric' it returns logical NA (!!!)
plotProteins <- function(pdat, protein=protein, log=FALSE, ymin=as.numeric(NA), ymax=as.numeric(NA),
                         text.size=12, point.size=3, title=NULL) {
  sel <- which(pdat$proteins %in% protein)
  if(length(sel) > 0 && sel > 0) {
    E <- if(log) log10(pdat$tab[sel,]) else pdat$tab[sel,]
    e <- colMeans(E, na.rm=TRUE)
    s <- sapply(E, function(x) sd(x, na.rm=TRUE)/sqrt(length(x)))
    n <- length(sel)

    if(is.null(title)) {
      title <- ifelse(n == 1, protein, paste0("selection of ", n, " proteins."))
    }

    p <- data.frame(
      expr = e,
      lo = e - s,
      up = e + s,
      condition = factor(meta$condition, levels=unique(meta$condition)),
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
#' Calls limma to perform differential expression.
#' @param pdat Protein intensity structure.
#' @param formula A string with a formula for building the linear model. The default value is "~condition".
#' @return limma output from eBays. See limma documentation for more details.
limmaDE <- function(pdat, formula="~condition") {
  design <- model.matrix(as.formula(formula), pdat$metadata)
  tab <- log10(pdat$tab)
  fit <- lmFit(tab, design)
  ebay <- eBayes(fit)
}

#' Create differential expression result.
#'
#' Use limmaDE output to create a table with DE results.
#' @param pdat Protein intensity structure.
#' @param ebay Output from "limmaDE"
#' @param column Which column should be used to extract data. The default value is "condition".
#' @return A data frame with DE results.
limmaTable <- function(pdat, ebay, column="condition") {
  # levels from the column (e.g. conditions)
  levs <- levels(factor(pdat$metadata[[column]]))
  # coef is the column name + the last level, e.g. "conditionWT"
  coef <- paste0(column, levs[-1])
  res <- topTable(ebay, coef=coef, adjust="BH", sort.by="none", number=1e9)
  res <- cbind(protein=rownames(res), res)
  rownames(res) <- c()
  return(res)
}

#' MA plot
#'
#' log fold change versus log sum plot.
#' @param pdat Protein data structure.
#' @param pair A two-element vector containing the pair of conditions to use. Can be skipped if there are only two conditions.
#' @param pvalue A vector with corresponding p-values for an interactive plotly plot
#' @param bins Number of bins for binhex
#' @param marginal.histograms A logical to add marginal histograms
#' @param xmin Lower limit on x-axis
#' @param xmax Upper limit on x-axis
#' @param ymax Upper limit on y-axis. If used, the lower limit is -ymax
#' @param text.size Text size
#' @param show.legend Logical to show legend (colour key)
#' @param plot.grid Logical to plot grid
plotMA <- function(pdat, pair=NULL, pvalue=NULL, bins=80, marginal.histograms=FALSE, classic=FALSE,
                   xmin=NULL, xmax=NULL, ymax=NULL, text.size=12, show.legend=TRUE, plot.grid=TRUE,
                   binhex=TRUE) {
  if(is.null(pair)) pair <- levels(pdat$metadata$condition)
  if(length(pair) != 2) stop("Need exactly two conditions. You might need to specify pair.")

  stats <- pdat$stats
  c1 <- pair[1]
  c2 <- pair[2]
  m1 <- log10(stats[which(stats$condition == c1),]$mean)
  m2 <- log10(stats[which(stats$condition == c2),]$mean)
  d <- data.frame(x = (m1 + m2) / 2, y = m2 - m1)
  d$id <- pdat$proteins
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
#' Plot distribution of raw p-value, obtained by limmaDE.
#' @param res Result table from limmaTable
#' @param text.size Text size
#' @param plot.grid Logical to plot grid
plotPdist <- function(res, bin.size=0.02, text.size=12, plot.grid=TRUE) {
  ggplot(res, aes(P.Value, ..density..)) +
    {if(plot.grid) simple_theme_grid else simple_theme} +
    geom_histogram(breaks=seq(0, 1, bin.size), colour='blue') +
    labs(x='P-value', y='Density') +
    theme(text = element_text(size=text.size))
}

#' Volcano plot
#'
#' Volcano plot from limma results.
#' @param res Result table from limmaTable
#' @param bins Number of bins for binhex
#' @param xmax Upper limit on x-axis. If used, the lower limit is -xmax
#' @param ymax Upper limit on y-axis. If used, the lower limit is -ymax
#' @param text.size Text size
#' @param show.legend Logical to show legend (colour key)
#' @param plot.grid Logical to plot grid
plotVolcano <- function(res, bins=80, xmax=NULL, ymax=NULL, text.size=12, show.legend=TRUE, plot.grid=TRUE) {
  g <- ggplot(res, aes(logFC, -log10(P.Value))) +
    {if(plot.grid) simple_theme_grid else simple_theme} +
    stat_binhex(bins=bins, show.legend=show.legend) +
    scale_fill_gradientn(colours=c("green","yellow", "red"), name = "count",na.value=NA) +
    geom_vline(colour='red', xintercept=0) +
    theme(text = element_text(size=text.size))
    # labs(title=title, x=paste0(c1, '+', c2), y=paste0(c2, '-', c1))

    if(!is.null(xmax)) g <- g + scale_x_continuous(limits = c(-xmax, xmax), expand = c(0, 0))
    if(!is.null(ymax) ) g <- g + scale_y_continuous(limits = c(0, ymax), expand = c(0, 0))
  return(g)
}



plotProtPeptides <- function(pepdat, protein, prodat=NULL) {
  # all peptides for this protein
  peps <- pepdat$pep2prot[which(pepdat$pep2prot$protein == protein),'sequence']
  mat <- as.matrix(pepdat$tab[peps,])
  dat <- melt(mat, varnames=c("peptide", "sample"))
  dat$pepnum <- sprintf("%02d", as.numeric(dat$peptide))  # convert sequences into numbers
  dat$intensity <- log10(dat$value)

  # add condition column (is there a simpler way?)
  s2c <- select(pepdat$metadata, condition)
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
    theme(legend.position="none")
  g2 <- ggplot(dat, aes(x=sample, y=intensity, fill=condition)) +
    geom_boxplot(outlier.shape = NA)  +
    geom_jitter(width=0, size=0.5) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) +
    theme(legend.position="none")
  if(!is.null(prodat)) g2 <- g2 + geom_point(aes(x=sample, y=prot.intensity), shape=22, size=3, fill='white')
  grid.arrange(g1, g2, ncol=1)
}
