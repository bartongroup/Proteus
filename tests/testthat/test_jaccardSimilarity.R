library(testthat)

# test data
x1 <- c(10, 20, 30, NA, 40, NA)
y1 <- c(NA, 10, 20, 30, 99, NA)

# expected result
sim1 <- 3/5

# no-intersection data
x2 <- c(10, 20, 30, NA, NA, NA)
y2 <- c(NA, NA, NA, 30, 99, 12)

# expected result
sim2 <- 0

# no-union data
x3 <- c(10, 20, 30, NA, NA, NA)
y3 <- c(NA, NA, NA, NA, NA, NA)

# expected result
sim3 <- 0

# wrong sizes
x4 <- c(10, 20, 30, NA, 40)
y4 <- c(NA, 10, 20, 30, 99, NA)



context("Jaccard similarity")

test_that("Correct results", {
  expect_equal(jaccardSimilarity(x1, y1), sim1)
  expect_equal(jaccardSimilarity(x2, y2), sim2)
  expect_equal(jaccardSimilarity(x3, y3), sim3)
})

test_that("Wrong size", {
  expect_error(jaccardSimilarity(x4, y4))
})
