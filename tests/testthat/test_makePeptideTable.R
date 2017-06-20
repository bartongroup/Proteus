library(testthat)

# expected result
tab <- structure(
  list(
    WT1 = c(3, 7, 7, 2, 5, 5, 5, 5, 8, 7),
    WT2 = c(8, NA, 2, 4, 4, 9, 2, 2, NA, NA),
    KO1 = c(9, 5, 4, NA, 7, NA, 2, 4, 5, 3)
  ),
  .Names = c("WT1", "WT2", "KO1"),
  row.names = c("AA", "AB", "AC", "AD", "AE", "BA", "BB", "BC", "BD", "BE"),
  class = "data.frame"
)

# input data:
evi <- read.table("../testdata/data_makePeptide_evi.txt", header=TRUE, sep="\t")
meta <- read.table("../testdata/data_makePeptide_meta.txt", header=TRUE, sep="\t")


context("Casting evidence into peptide table")

test_that("Test makePeptide", {
  pep <- makePeptideTable(evi, meta)
  expect_equal(pep$tab, tab)
  expect_equal(pep$metadata, meta)
  expect_equal(pep$content, "peptide")
  expect_equal(pep$intensity, "intensity")
  expect_equal(pep$peptides, row.names(tab))
  expect_equal(pep$proteins, c("A", "B"))
  expect_true(is(pep, "proteusData"))
})
