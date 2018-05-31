library(testthat)

# expected result
tab.1 <- structure(
  c(7, 7, 4, 2, 5, 5, 5, 5, 8, 7,
    8, NA, 2, 4, 4, 9, 2, 2, NA, NA,
    9, 5, 4, NA, 7, NA, 2, 4, 5, 3),
  .Dim = c(10L, 3L),
  .Dimnames = list(
    c("AA", "AB", "AC", "AD", "AE", "BA", "BB", "BC", "BD", "BE"),
    c("WT1", "WT2", "KO1")
  )
)

tab.2 <- structure(
  c(3, 8, 2, 7, 2, 8, 5, 2, 2, NA,
    4, 1, 2, NA, 1, 7, NA, 2, NA, NA,
    1, NA, 6, NA, 5, 6, 1, 7, 3, 2),
  .Dim = c(10L, 3L),
  .Dimnames = list(
    c("AA", "AB", "AC", "AD", "AE", "BA", "BB", "BC", "BD", "BE"),
    c("WT1", "WT2", "KO1")
  )
)

# input data:
evi <- read.table("../testdata/data_makePeptide_evi.txt", header=TRUE, sep="\t")
meta <- read.table("../testdata/data_makePeptide_meta.txt", header=TRUE, sep="\t")

# make sure even if levels order is not correct, makePeptideTable sorts this out:
meta.ordered <- meta
meta.ordered$sample <- factor(meta.ordered$sample, levels=meta.ordered$sample)

context("Casting evidence into peptide table")

test_that("Test makePeptide unlabelled", {
  pep <- makePeptideTable(evi, meta, ncores=1)
  expect_equal(pep$tab, tab.1)
  expect_equal(pep$metadata[,c("experiment", "measure", "sample", "condition", "batch")], meta.ordered)
  expect_equal(pep$content, "peptide")
  expect_equal(pep$measures, "Intensity")
  expect_equal(pep$peptides, row.names(tab.1))
  expect_equal(pep$proteins, c("A", "B"))
  expect_true(is(pep, "proteusData"))
})

meta$measure <- "Ratio"
test_that("Test makePeptide SILAC", {
  pep <- makePeptideTable(evi, meta, aggregate.fun = aggregateMedian, measure.cols=c(ratio="Ratio"), experiment.type="SILAC", ncores=1)
  expect_equal(pep$tab, tab.2)
  expect_equal(pep$metadata$batch, c(1, 2, 2))
  expect_equal(pep$measures, "Ratio")
  expect_equal(pep$peptides, row.names(tab.2))
  expect_equal(pep$proteins, c("A", "B"))
  expect_true(is(pep, "proteusData"))
})
