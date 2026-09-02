source(testthat::test_path("helper-load.R"), local = .GlobalEnv)

validation_path <- testthat::test_path(
  "../..", "R", "workflow_validation.R"
)
if (file.exists(validation_path)) {
  source(validation_path, local = .GlobalEnv)
}

testthat::test_that("proportion mode requires exactly one optional input", {
  testthat::expect_true(exists("validate_proportion_mode", mode = "function"))
  if (!exists("validate_proportion_mode", mode = "function")) {
    return(invisible(NULL))
  }

  dtangle <- validate_proportion_mode(TRUE, FALSE)
  precomputed <- validate_proportion_mode(FALSE, TRUE)
  testthat::expect_identical(
    dtangle,
    list(selected_mode = "dtangle", estimate_proportions = TRUE)
  )
  testthat::expect_identical(
    precomputed,
    list(selected_mode = "precomputed", estimate_proportions = FALSE)
  )
  testthat::expect_error(
    validate_proportion_mode(TRUE, TRUE),
    "exactly one.*lm22.*precomputed_proportions"
  )
  testthat::expect_error(
    validate_proportion_mode(FALSE, FALSE),
    "exactly one.*lm22.*precomputed_proportions"
  )
})

run_validation_cli <- function(lm22_defined, precomputed_defined) {
  script <- testthat::test_path(
    "../..", "scripts", "validate_proportion_mode.R"
  )
  if (!file.exists(script)) {
    return(list(status = 127L, output = "validation script is absent"))
  }
  output_directory <- tempfile("proportion-mode-")
  dir.create(output_directory)
  output <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c(
      script,
      "--lm22-defined", tolower(as.character(lm22_defined)),
      "--precomputed-defined", tolower(as.character(precomputed_defined)),
      "--output-dir", output_directory
    ),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(output, "status")
  if (is.null(status)) {
    status <- 0L
  }
  list(status = as.integer(status), output = paste(output, collapse = "\n"))
}

testthat::test_that("mode validation CLI rejects both optional inputs without Docker", {
  result <- run_validation_cli(TRUE, TRUE)

  testthat::expect_false(identical(result$status, 0L))
  testthat::expect_match(
    result$output,
    "exactly one.*lm22.*precomputed_proportions"
  )
})

testthat::test_that("mode validation CLI rejects neither optional input without Docker", {
  result <- run_validation_cli(FALSE, FALSE)

  testthat::expect_false(identical(result$status, 0L))
  testthat::expect_match(
    result$output,
    "exactly one.*lm22.*precomputed_proportions"
  )
})
