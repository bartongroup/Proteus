library(testthat)

# expected result
tab <- structure(
  c(3, 7, 7, 2, 5, 5, 5, 5, 8, 7, 8, NA, 2, 4, 4, 9, 2, 2, NA, NA, 9, 5, 4, NA, 7, NA, 2, 4, 5, 3),
  .Dim = c(10L, 3L),
  .Dimnames = list(
    c("AA", "AB", "AC", "AD", "AE", "BA", "BB", "BC", "BD", "BE"),
    c("WT1", "WT2", "KO1")
  )
)

# median of 3 and 4 is 3.5
tab.med <- tab
tab.med['AC', 'WT1'] <- 3.5

# input data:
evi <- read.table("../testdata/data_makePeptide_evi.txt", header=TRUE, sep="\t")
meta <- read.table("../testdata/data_makePeptide_meta.txt", header=TRUE, sep="\t")


context("Casting evidence into peptide table")

test_that("Test makePeptide", {
  pep <- makePeptideTable(evi, meta)
  expect_equal(pep$tab, tab)
  expect_equal(pep$metadata, meta)
  expect_equal(pep$content, "peptide")
  expect_equal(pep$value, "intensity")
  expect_equal(pep$peptides, row.names(tab))
  expect_equal(pep$proteins, c("A", "B"))
  expect_true(is(pep, "proteusData"))
})

test_that("Test makePeptide SILAC", {
  pep <- makePeptideTable(evi, meta, aggregate.fun = median, value="ratio")
  expect_equal(pep$tab, tab.med)
  expect_equal(pep$value, "ratio")
})
