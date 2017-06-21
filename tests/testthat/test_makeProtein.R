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

test_that("Test makeProtein sum", {
  p <- makeProtein(test_wp, method="sum")
  exp_p <- data.frame(t(c(NA, 10, 30, 32)))
  colnames(exp_p) <- cols
  expect_equal(p, exp_p)
})

test_that("Test makeProtein hifly longer", {
  p <- makeProtein(test_wp, method="hifly", hifly=3)
  exp_p <- data.frame(t(c(NA, 10, 10, mean(c(10,8,8)))))
  colnames(exp_p) <- cols
  expect_equal(p, exp_p)
})

test_that("Test makeProtein hifly same", {
  p <- makeProtein(test_wp[1:3,], method="hifly", hifly=3)
  exp_p <- data.frame(t(c(NA, 10, 9.5, mean(c(1,5,10)))))
  colnames(exp_p) <- cols
  expect_equal(p, exp_p)
})

test_that("Test makeProtein hifly one", {
  wp <- as.data.frame(t(c(NA, 1, 5, 10)))
  colnames(wp) <- cols
  p <- makeProtein(wp, method="hifly", hifly=3)
  expect_equal(p, wp)
})

test_that("Test makeProtein no peptides",{
  wp <- test_wp[integer(0),]
  expect_error(makeProtein())
})
