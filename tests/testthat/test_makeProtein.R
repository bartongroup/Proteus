library(testthat)

cols <- c("R1", "R2", "R3", "R4")
test_wp <- as.matrix(rbind(
  c(NA, NA, NA, 1),
  c(NA, 10, 10, 5),
  c(NA, NA, 9, 10),
  c(NA, NA, 11, 8),
  c(NA, NA, NA, 8)
))
colnames(test_wp) <- cols

context("Aggregating peptides into proteins")

test_that("Test aggregateSum", {
  p <- aggregateSum(test_wp)
  exp_p <- c(NA, 10, 30, 32)
  expect_equal(p, exp_p)
})

test_that("Test aggregateMedian", {
  p <- aggregateMedian(test_wp)
  exp_p <- c(NA, 10, 10, 8)
  expect_equal(p, exp_p)
})

test_that("Test aggregateHifly hifly longer", {
  p <- aggregateHifly(test_wp, hifly=3)
  exp_p <- c(NA, 10, 10, mean(c(10,8,8)))
  expect_equal(p, exp_p)
})

test_that("Test aggregateHifly hifly same", {
  p <- aggregateHifly(test_wp[1:3,], hifly=3)
  exp_p <- c(NA, 10, 9.5, mean(c(1,5,10)))
  expect_equal(p, exp_p)
})

test_that("Test aggregateHifly hifly one", {
  wp <- as.matrix(t(c(NA, 1, 5, 10)))
  p <- aggregateHifly(wp, hifly=3)
  expect_equal(p, as.vector(wp))
})

test_that("Test aggregateHifly no peptides",{
  wp <- test_wp[integer(0),]
  expect_error(aggregateHifly(wp, hifly=3))
})
