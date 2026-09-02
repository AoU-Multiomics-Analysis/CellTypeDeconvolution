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
    c(
      script_path,
      "--expression", "expression.bed",
      "--lm22", "lm22.tsv",
      "--output-dir", "outputs"
    ),
    stdout = TRUE,
    stderr = TRUE
  ))

  testthat::expect_true(!is.null(attr(result, "status")))
  testthat::expect_match(
    paste(result, collapse = "\n"),
    "Missing required options: gtf",
    fixed = TRUE
  )
})

testthat::test_that("dtangle CLI rejects the retired prepared-log option", {
  script_path <- testthat::test_path("../..", "scripts", "run_dtangle.R")

  result <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c(
      script_path,
      "--bulk-log", "bulk.tsv"
    ),
    stdout = TRUE,
    stderr = TRUE
  ))

  testthat::expect_true(!is.null(attr(result, "status")))
  testthat::expect_match(
    paste(result, collapse = "\n"),
    'long flag "bulk-log" is invalid',
    fixed = TRUE
  )
})
