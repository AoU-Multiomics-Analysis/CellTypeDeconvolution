wdl_text <- function(path) {
  paste(
    readLines(testthat::test_path("../..", path), warn = FALSE),
    collapse = "\n"
  )
}

expect_task_contract <- function(path, task, inputs, outputs, defaults) {
  text <- wdl_text(path)
  testthat::expect_match(text, paste0("task ", task))
  purrr::walk(inputs, ~ testthat::expect_match(text, .x, fixed = TRUE))
  purrr::walk(outputs, ~ testthat::expect_match(text, .x, fixed = TRUE))
  purrr::walk(defaults, ~ testthat::expect_match(text, .x, fixed = TRUE))
  testthat::expect_match(text, "set -euo pipefail", fixed = TRUE)
  testthat::expect_match(text, "stage=", fixed = TRUE)
  testthat::expect_match(text, "start_time=", fixed = TRUE)
  testthat::expect_match(text, "completion_time=", fixed = TRUE)
  testthat::expect_match(text, "dimensions=", fixed = TRUE)
  testthat::expect_match(text, "outputs=", fixed = TRUE)
  testthat::expect_match(text, "error_status=", fixed = TRUE)
  testthat::expect_match(text, "preemptible_attempts", fixed = TRUE)
  testthat::expect_match(text, "max_retries", fixed = TRUE)
}

testthat::test_that("the expression WDL requires CPM and GTF", {
  text <- wdl_text("workflows/tasks/expression.wdl")
  testthat::expect_match(text, "File expression", fixed = TRUE)
  testthat::expect_match(text, "File gtf", fixed = TRUE)
  testthat::expect_match(text, "prepared_cpm", fixed = TRUE)
  testthat::expect_match(text, "prepared_log2_cpm", fixed = TRUE)
  testthat::expect_false(grepl("expression_type|gene_length", text))
})

testthat::test_that("all modular WDL tasks declare CLI-aligned interfaces", {
  expect_task_contract(
    "workflows/tasks/expression.wdl", "PrepareExpression",
    c("File expression", "File gtf", "String docker_image"),
    c("File prepared_cpm", "File prepared_log2_cpm", "File mapping_report",
      "File excluded_genes", "File log"),
    c("Int cpu = 4", "String memory = \"64 GB\"", "Int disk_gb = 400",
      "Int preemptible_attempts = 2", "Int max_retries = 2")
  )
  expect_task_contract(
    "workflows/tasks/dtangle.wdl", "RunDtangle",
    c("File prepared_log2_cpm", "File lm22", "Float min_overlap",
      "Float marker_fraction", "String marker_method", "Boolean quantile_normalize"),
    c("File proportions", "File markers", "File metadata", "File overlap_report",
      "File transformed_lm22", "File shared_bulk", "File log"),
    c("Int cpu = 4", "String memory = \"32 GB\"", "Int disk_gb = 100",
      "Int preemptible_attempts = 2", "Int max_retries = 2")
  )
  expect_task_contract(
    "workflows/tasks/proportions.wdl", "ProcessProportions",
    c("File proportions", "Float mean_threshold", "Float zero_floor"),
    c("File original", "File combined", "File tca_weights", "File filter_report",
      "File log"),
    c("Int cpu = 2", "String memory = \"16 GB\"", "Int disk_gb = 50",
      "Int preemptible_attempts = 2", "Int max_retries = 2")
  )
  tca_text <- wdl_text("workflows/tasks/tca.wdl")
  testthat::expect_match(tca_text, "task FitTca", fixed = TRUE)
  testthat::expect_match(tca_text, "task ExtractTcaShard", fixed = TRUE)
  purrr::walk(c(
    "File prepared_log2_cpm", "File tca_weights", "File? covariates",
    "Int shard_size", "Int max_iters", "Int random_seed", "File model",
    "File model_log", "File tca_expression", "File excluded_genes",
    "File shard_manifest", "Array[File] shards", "File shard", "File shard_hdf5",
    "Int cpu = 16", "String memory = \"192 GB\"", "Int disk_gb = 750",
    "Int preemptible_attempts = 0", "Int max_retries = 1",
    "Int cpu = 8", "String memory = \"64 GB\"", "Int disk_gb = 200",
    "Int preemptible_attempts = 2", "Int max_retries = 2"
  ), ~ testthat::expect_match(tca_text, .x, fixed = TRUE))
  qc_text <- wdl_text("workflows/tasks/qc.wdl")
  purrr::walk(c(
    "task AssembleTca", "task BuildManifest", "Array[File] shard_hdf5",
    "File shard_manifest", "File tca_expression", "File model", "File tca_weights",
    "File? covariates", "Boolean write_tsv", "String pipeline_version",
    "Array[File] group_hdf5", "Array[File] group_tsv",
    "File reconstruction_by_sample", "File assembly_qc", "File output_manifest",
    "File qc_summary", "File qc_plots", "File provenance",
    "Int cpu = 8", "String memory = \"128 GB\"", "Int disk_gb = 500",
    "Int preemptible_attempts = 0", "Int max_retries = 1",
    "Int cpu = 4", "String memory = \"32 GB\"", "Int disk_gb = 100",
    "Int preemptible_attempts = 1", "Int max_retries = 2"
  ), ~ testthat::expect_match(qc_text, .x, fixed = TRUE))
})

testthat::test_that("the logging checker rejects incomplete command logging", {
  fixture <- tempfile(fileext = ".wdl")
  checker_path <- testthat::test_path(
    "../..", "tools", "check_wdl_logging.py"
  )
  writeLines(c(
    "version 1.1", "task MissingLogging {", "  command <<<", "    echo hello", "  >>>", "}"
  ), fixture)
  status <- system2(
    "python3", c(checker_path, fixture), stdout = FALSE,
    stderr = FALSE
  )
  testthat::expect_false(identical(status, 0L))
})

testthat::test_that("the logging checker checks brace-form commands", {
  incomplete_fixture <- tempfile(fileext = ".wdl")
  complete_fixture <- tempfile(fileext = ".wdl")
  checker_path <- testthat::test_path(
    "../..", "tools", "check_wdl_logging.py"
  )
  writeLines(c(
    "version 1.1", "task MissingBraceLogging {", "  command {",
    "    echo 'stage=brace start_time=now'", "  }", "}"
  ), incomplete_fixture)
  writeLines(c(
    "version 1.1", "task CompleteBraceLogging {", "  command {",
    "    echo 'stage=brace start_time=now completion_time=now dimensions=1 outputs=result'",
    "  }", "}"
  ), complete_fixture)

  incomplete_status <- system2(
    "python3", c(checker_path, incomplete_fixture),
    stdout = FALSE, stderr = FALSE
  )
  complete_status <- system2(
    "python3", c(checker_path, complete_fixture),
    stdout = FALSE, stderr = FALSE
  )
  miniwdl_status <- system2(
    "miniwdl", c("check", complete_fixture), stdout = FALSE, stderr = FALSE
  )

  testthat::expect_false(identical(incomplete_status, 0L))
  testthat::expect_identical(complete_status, 0L)
  testthat::expect_identical(miniwdl_status, 0L)
})
