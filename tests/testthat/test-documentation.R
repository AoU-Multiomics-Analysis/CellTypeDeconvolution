testthat::test_that("both examples use BED and GTF workflow inputs", {
  files <- list.files(
    testthat::test_path("../..", "examples"),
    pattern = "[.]inputs[.]json$",
    full.names = TRUE
  )
  testthat::expect_length(files, 2L)
  testthat::expect_setequal(
    basename(files),
    c("bed.inputs.json", "precomputed-proportions.inputs.json")
  )

  parsed <- purrr::map(files, jsonlite::read_json)
  testthat::expect_true(all(purrr::map_lgl(parsed, ~ all(c(
    "cell_type_deconvolution.expression",
    "cell_type_deconvolution.gtf"
  ) %in% names(.x)))))
  removed_parameters <- c(
    "tca_shard_size", "write_tsv", "hdf5_gene_chunk_max",
    "hdf5_sample_chunk_max", "hdf5_gzip_level", "extract_cpu",
    "extract_memory", "extract_disk_gb", "extract_preemptible_attempts",
    "extract_max_retries", "assemble_cpu", "assemble_memory",
    "assemble_disk_gb", "assemble_preemptible_attempts",
    "assemble_max_retries"
  )
  purrr::walk(parsed, function(inputs) {
    testthat::expect_false(any(paste0(
      "cell_type_deconvolution.", removed_parameters
    ) %in% names(inputs)))
  })
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

testthat::test_that("active dependencies and guides use direct BED outputs", {
  environment <- paste(readLines(
    testthat::test_path("../..", "envs", "environment.yml"),
    warn = FALSE
  ), collapse = "\n")
  guides <- paste(vapply(
    c("README.md", "docs/terra.md", "docs/data-dictionary.md"),
    function(path) paste(readLines(
      testthat::test_path("../..", path), warn = FALSE
    ), collapse = "\n"),
    character(1)
  ), collapse = "\n")
  testthat::expect_false(grepl("r-hdf5r", environment, fixed = TRUE))
  testthat::expect_match(guides, "#chr.*start.*end.*gene_id")
  testthat::expect_match(guides, "(?i)one.*BED[.]GZ.*retained")
  testthat::expect_match(guides, "log2_cpm", fixed = TRUE)
  testthat::expect_match(guides, "does not produce HDF5", fixed = TRUE)
  testthat::expect_false(grepl(
    "group_hdf5|tensor_shards|gene_shard_manifest|write_tsv",
    guides
  ))
})
