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
    image <- inputs[["cell_type_deconvolution.docker_image"]]
    testthat::expect_match(
      image,
      "@sha256:REPLACE_WITH_64_LOWERCASE_HEX_DIGEST$"
    )
    testthat::expect_false(grepl(":latest", image, fixed = TRUE))
  })
})

testthat::test_that("smoke inputs use the simplified global runtime controls", {
  fixtures <- c("dtangle.inputs.json", "restart.inputs.json")
  removed_inputs <- c(
    "pipeline_version",
    "prepare_cpu", "prepare_memory", "prepare_disk_gb",
    paste0(
      c("prepare", "dtangle", "proportions", "fit", "export", "manifest"),
      rep(c("_preemptible_attempts", "_max_retries"), each = 6L)
    )
  )

  purrr::walk(fixtures, function(fixture) {
    inputs <- jsonlite::read_json(testthat::test_path("..", "fixtures", fixture))
    testthat::expect_identical(
      inputs$cell_type_deconvolution.preemptible_attempts,
      2L
    )
    testthat::expect_identical(inputs$cell_type_deconvolution.max_retries, 2L)
    testthat::expect_false("cell_type_deconvolution.pipeline_version" %in% names(inputs))
    testthat::expect_false(any(paste0(
      "cell_type_deconvolution.", removed_inputs
    ) %in% names(inputs)))
    testthat::expect_identical(
      inputs$cell_type_deconvolution.docker_image,
      "celltype-deconvolution:test"
    )
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
  testthat::expect_match(text, "64 lowercase hexadecimal")
  testthat::expect_match(text, "celltype-deconvolution:test.*only")
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

testthat::test_that("data dictionary uses the implemented overlap interval", {
  text <- paste(
    readLines(testthat::test_path("../..", "docs", "data-dictionary.md")),
    collapse = "\n"
  )

  testthat::expect_match(
    text,
    "min_lm22_overlap.*greater than 0 and no greater than 1"
  )
  testthat::expect_false(grepl("Finite value from 0 through 1", text, fixed = TRUE))
})
