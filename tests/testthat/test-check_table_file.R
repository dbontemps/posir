write_a_table <- function(dirpath) {
  newone <- paste(dirpath, "mytemptable.txt", sep = "/")
  X <- matrix(1, 50, 6)
  colnames(X) <- seq(7, 2, -1) / 8
  write.table(X, newone)
  return(newone)
}

test_that("Validation of a correct table", {
  logger::log_debug(paste("test-check_table_file:",
                         "Validation of a correct table."),
                   namespace = "posir")
  sim_dir = local_init_testing(dirNameBase = "posir_test")
  newone <- write_a_table(sim_dir)
  expect_true(check_table_file(newone, seq(7, 2, -1) / 8, 50))
})

test_that("Non-validation if nrow differs", {
  logger::log_debug(paste("test-check_table_file:",
                         "Non-validation if nrow differs."),
                   namespace = "posir")
  sim_dir = local_init_testing(dirNameBase = "posir_test")
  newone <- write_a_table(sim_dir)
  expect_false(check_table_file(newone, seq(7, 2, -1) / 8, 100))
})

test_that("Non-validation if grid length differs", {
  logger::log_debug(paste("test-check_table_file:",
                         "Non-validation if grid length differs."),
                   namespace = "posir")
  sim_dir = local_init_testing(dirNameBase = "posir_test")
  newone <- write_a_table(sim_dir)
  expect_false(check_table_file(newone, seq(8, 2, -1) / 8, 50))
})

test_that("Non-validation if grid content differs", {
  logger::log_debug(paste("test-check_table_file:",
                         "Non-validation if grid content differs."),
                   namespace = "posir")
  sim_dir = local_init_testing(dirNameBase = "posir_test")
  newone <- write_a_table(sim_dir)
  check_table_file(newone, c(1, seq(6, 2, -1) / 8), 50)
  expect_false(check_table_file(newone, c(1, seq(6, 2, -1) / 8), 50))
})

test_that("Non-validation if no file", {
  logger::log_debug(paste("test-check_table_file:",
                         "Non-validation if no file."),
                   namespace = "posir")
  sim_dir = local_init_testing(dirNameBase = "posir_test")
  newone <- paste(sim_dir, "mytemptable.txt", sep = "/")
  if (file.exists(newone)) file.remove(newone)
  expect_false(check_table_file(newone, seq(7, 2, -1) / 8, 50))
})
