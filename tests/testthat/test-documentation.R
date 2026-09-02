testthat::test_that("both examples use CPM and GTF workflow inputs", {
  files <- list.files(
    testthat::test_path("../..", "examples"),
    pattern = "[.]inputs[.]json$",
    full.names = TRUE
  )
  testthat::expect_length(files, 2L)

  parsed <- purrr::map(files, jsonlite::read_json)
  testthat::expect_true(all(purrr::map_lgl(parsed, ~ all(c(
    "cell_type_deconvolution.expression",
    "cell_type_deconvolution.gtf"
  ) %in% names(.x)))))
})

testthat::test_that("Terra guidance gives the input scale and LM22 contract", {
  text <- paste(
    readLines(testthat::test_path("../..", "docs", "terra.md")),
    collapse = "\n"
  )

  testthat::expect_match(text, "user-supplied")
  testthat::expect_match(text, "linear LM22")
  testthat::expect_match(text, "log2\\(CPM\\)")
  testthat::expect_match(text, "strictly positive")
  testthat::expect_false(grepl("expression_type|log2_tpm", text))
})
