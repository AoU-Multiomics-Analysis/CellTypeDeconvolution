testthat::test_that("dtangle CLI reads direct CPM BED and GTF", {
  script_path <- testthat::test_path("../..", "scripts", "run_dtangle.R")
  text <- paste(readLines(script_path, warn = FALSE), collapse = "\n")
  testthat::expect_match(text, '"--expression"', fixed = TRUE)
  testthat::expect_match(text, '"--gtf"', fixed = TRUE)
})

testthat::test_that("dtangle CLI rejects incomplete direct CPM input mode", {
  script_path <- testthat::test_path("../..", "scripts", "run_dtangle.R")

  result <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c(script_path, "--expression", "expression.bed"),
    stdout = TRUE,
    stderr = TRUE
  ))

  testthat::expect_true(!is.null(attr(result, "status")))
  testthat::expect_match(
    paste(result, collapse = "\n"),
    "requires both --expression and --gtf",
    fixed = TRUE
  )
})

testthat::test_that("dtangle CLI rejects conflicting bulk and direct input modes", {
  script_path <- testthat::test_path("../..", "scripts", "run_dtangle.R")

  result <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c(
      script_path,
      "--bulk-log", "bulk.tsv",
      "--expression", "expression.bed",
      "--gtf", "annotation.gtf"
    ),
    stdout = TRUE,
    stderr = TRUE
  ))

  testthat::expect_true(!is.null(attr(result, "status")))
  testthat::expect_match(
    paste(result, collapse = "\n"),
    "either --bulk-log alone or --expression plus --gtf",
    fixed = TRUE
  )
})
