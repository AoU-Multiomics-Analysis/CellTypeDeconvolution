args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L) {
  testthat::test_dir("tests/testthat", reporter = "summary")
} else {
  purrr::walk(args, ~ testthat::test_file(.x, reporter = "summary"))
}
