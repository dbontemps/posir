test_that("a correct grid is validated", {
  logger::log_debug(paste("test-check_grid:",
                         "a correct grid is validated."),
                   namespace = "posir")
  grille <- seq(20, 1, -1) / 20
  N <- 20
  expect_equal(check_grid(N, grille, .01), grille * N)
})

test_that("grid upper bound", {
  logger::log_debug(paste("test-check_grid:",
                         "grid upper bound."),
                   namespace = "posir")
  expect_error(check_grid(10, seq(20, 1, -1) / 10, .01))
})

test_that("grid lower bound", {
  logger::log_debug(paste("test-check_grid:",
                         "grid lower bound."),
                   namespace = "posir")
  expect_error(check_grid(20, seq(20, 0, -1) / 20, .01))
})

test_that("grid subset of (1/N)Z", {
  logger::log_debug(paste("test-check_grid:",
                         "grid subset of (1/N)Z."),
                   namespace = "posir")
  expect_error(check_grid(10, seq(20, 1, -1) / 20, .01))
})

test_that("grid sorted in decreasing order", {
  logger::log_debug(paste("test-check_grid:",
                         "grid sorted in decreasing order."),
                   namespace = "posir")
  expect_error(check_grid(10, c(1, .5, .6, .1), .01))
})
