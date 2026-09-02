source(testthat::test_path("helper-load.R"), local = .GlobalEnv)

expression_path <- testthat::test_path("../..", "R", "expression.R")
if (file.exists(expression_path)) {
  source(expression_path, local = .GlobalEnv)
}
expression_bed_path <- testthat::test_path("../..", "R", "expression_bed.R")
if (file.exists(expression_bed_path)) {
  source(expression_bed_path, local = .GlobalEnv)
}
bed_outputs_path <- testthat::test_path("../..", "R", "bed_outputs.R")
if (file.exists(bed_outputs_path)) {
  source(bed_outputs_path, local = .GlobalEnv)
}
qc_path <- testthat::test_path("../..", "R", "qc.R")
if (file.exists(qc_path)) {
  source(qc_path, local = .GlobalEnv)
}

testthat::test_that("cell-type BED files preserve coordinates and sample order", {
  coordinates <- tibble::tibble(
    `#chr` = c("chr2", "chr1"), start = c(20L, 10L),
    end = c(30L, 15L), gene_id = c("g2", "g1")
  )
  tensor <- list(
    "CD4 T cells" = matrix(
      c(1, 2, 3, 4), nrow = 2, byrow = TRUE,
      dimnames = list(c("g2", "g1"), c("S2", "S1"))
    ),
    "CD8 T cells" = matrix(
      c(5, 6, 7, 8), nrow = 2, byrow = TRUE,
      dimnames = list(c("g2", "g1"), c("S2", "S1"))
    )
  )
  result <- write_cell_type_beds(tensor, coordinates, tempfile())
  observed <- readr::read_tsv(
    result$paths[["CD4 T cells"]],
    show_col_types = FALSE,
    progress = FALSE
  )
  testthat::expect_identical(
    names(observed),
    c("#chr", "start", "end", "gene_id", "S2", "S1")
  )
  testthat::expect_identical(observed$gene_id, c("g2", "g1"))
  testthat::expect_equal(
    unname(as.matrix(observed[c("S2", "S1")])),
    unname(tensor[[1L]])
  )
  testthat::expect_identical(
    result$inventory$cell_group,
    c("CD4 T cells", "CD8 T cells")
  )
  testthat::expect_identical(
    result$inventory$path,
    c("cd4_t_cells.bed.gz", "cd8_t_cells.bed.gz")
  )
  testthat::expect_true(all(result$inventory$scale == "log2_cpm"))
})

testthat::test_that("direct CPM BED coordinates align to filtered TCA genes", {
  source_coordinates <- tibble::tibble(
    `#chr` = c("chr2", "chr1", "chr3"),
    start = c(20L, 10L, 30L),
    end = c(30L, 15L, 40L),
    gene_id = c("g2", "g1", "g3")
  )
  source_cpm <- matrix(
    c(4, 8, 16, 32, 64, 128),
    nrow = 3L,
    byrow = TRUE,
    dimnames = list(source_coordinates$gene_id, c("S1", "S2"))
  )
  expression_path <- tempfile(fileext = ".bed.gz")
  write_expression_bed(expression_path, source_coordinates, source_cpm)

  expression <- read_expression_bed(expression_path)
  modeled_gene_ids <- c("g2", "g3")
  coordinate_index <- match(modeled_gene_ids, expression$coordinates$gene_id)
  aligned_coordinates <- expression$coordinates[coordinate_index, , drop = FALSE]

  testthat::expect_identical(aligned_coordinates$gene_id, modeled_gene_ids)
  testthat::expect_identical(aligned_coordinates$start, c(20L, 30L))
})

testthat::test_that("TCA export reads a direct CPM BED", {
  text <- paste(readLines(
    testthat::test_path("../..", "scripts", "export_tca_beds.R"),
    warn = FALSE
  ), collapse = "\n")
  testthat::expect_match(text, '"--expression"', fixed = TRUE)
  testthat::expect_match(text, "read_expression_bed(options$expression)", fixed = TRUE)
})

testthat::test_that("tensor validation rejects group and identifier mismatches", {
  tensor <- list(A = matrix(
    1:4, nrow = 2,
    dimnames = list(c("g1", "g2"), c("S1", "S2"))
  ))
  testthat::expect_error(
    validate_tensor_contract(tensor, c("g1", "g2"), c("S1", "S2"), c("A", "B")),
    "cell groups"
  )
})

testthat::test_that("tensor validation rejects missing and empty source names", {
  source_matrix <- matrix(
    1:4,
    nrow = 2L,
    dimnames = list(c("g1", "g2"), c("S1", "S2"))
  )
  unnamed <- list(source_matrix)
  empty_named <- stats::setNames(list(source_matrix), "")

  testthat::expect_error(
    validate_tensor_contract(unnamed, c("g1", "g2"), c("S1", "S2"), NULL),
    "non-empty"
  )
  testthat::expect_error(
    validate_tensor_contract(empty_named, c("g1", "g2"), c("S1", "S2"), ""),
    "non-empty"
  )
})

testthat::test_that("full tensor extraction preserves model source and matrix order", {
  testthat::skip_if_not_installed("TCA")
  testthat::expect_identical(
    as.character(utils::packageVersion("TCA")), "1.2.1"
  )
  set.seed(20260901)
  data <- TCA::test_data(12, 16, 3, 0, 0, 0.01)
  model <- TCA::tca(
    data$X, data$W, refit_W = FALSE, max_iters = 2, verbose = FALSE
  )

  tensor <- extract_full_tensor(data$X, model, 1L, tempfile())

  testthat::expect_identical(names(tensor), colnames(model$W))
  testthat::expect_true(all(purrr::map_lgl(
    tensor,
    ~ identical(rownames(.x), rownames(data$X)) &&
      identical(colnames(.x), colnames(data$X)) &&
      identical(dim(.x), dim(data$X))
  )))
})

testthat::test_that("reconstruction includes the fitted C2 term", {
  tensor <- list(
    A = matrix(
      c(1, 3, 2, 4), nrow = 2,
      dimnames = list(c("g1", "g2"), c("S1", "S2"))
    ),
    B = matrix(
      c(5, 7, 6, 8), nrow = 2,
      dimnames = list(c("g1", "g2"), c("S1", "S2"))
    )
  )
  weights <- matrix(
    c(0.25, 0.75, 0.60, 0.40), nrow = 2, byrow = TRUE,
    dimnames = list(c("S1", "S2"), c("A", "B"))
  )
  C2 <- matrix(c(2, 3), ncol = 1, dimnames = list(c("S1", "S2"), "batch"))
  deltas_hat <- matrix(c(0.5, -1), ncol = 1, dimnames = list(c("g1", "g2"), "batch"))
  expected <- tensor$A |>
    sweep(2L, weights[, "A"], "*")
  expected <- expected + sweep(tensor$B, 2L, weights[, "B"], "*")
  expected <- expected + t(C2 %*% t(deltas_hat))

  observed <- reconstruct_tensor(tensor, weights, C2, deltas_hat)

  testthat::expect_equal(observed, expected)
})

testthat::test_that("reconstruction statistics combine genes by sample", {
  observed <- matrix(
    c(1, 2, 2, 4, 3, 5, 6, 10), nrow = 4, byrow = TRUE,
    dimnames = list(paste0("g", 1:4), c("S1", "S2"))
  )
  reconstructed <- observed + matrix(
    c(0, 1, 0, -1, 1, 0, -1, 0), nrow = 4, byrow = TRUE,
    dimnames = dimnames(observed)
  )
  stats <- initialize_reconstruction_stats(c("S1", "S2")) |>
    update_reconstruction_stats(observed, reconstructed)
  metrics <- finalize_reconstruction_stats(stats)

  testthat::expect_identical(metrics$sample_id, c("S1", "S2"))
  testthat::expect_equal(metrics$rmse, c(sqrt(2 / 4), sqrt(2 / 4)))
  testthat::expect_true(all(is.finite(metrics$correlation)))
})

testthat::test_that("manifest records hashes and dimensions", {
  output_path <- tempfile(fileext = ".bed.gz")
  writeBin(charToRaw("bed-output"), output_path)
  outputs <- tibble::tibble(
    logical_name = "cd4_t_cells_expression",
    path = output_path,
    n_genes = 2L,
    n_samples = 3L,
    scale = "log2_cpm",
    cell_group = "CD4 T cells"
  )

  manifest <- build_output_manifest(
    outputs = outputs,
    tca_version = "1.2.1",
    parameters = list(scale = "log2_cpm"),
    container_image = paste0(
      "example.org/pipeline@sha256:",
      paste(rep("a", 64L), collapse = "")
    )
  )

  testthat::expect_false("pipeline_version" %in% names(manifest))
  testthat::expect_identical(manifest$tca_version, "1.2.1")
  testthat::expect_identical(manifest$outputs[[1L]]$dimensions, c(2L, 3L))
  testthat::expect_identical(manifest$outputs[[1L]]$path, basename(output_path))
  testthat::expect_identical(manifest$outputs[[1L]]$file_name, basename(output_path))
  testthat::expect_identical(
    manifest$outputs[[1L]]$sha256,
    digest::digest(file = output_path, algo = "sha256", serialize = FALSE)
  )
  testthat::expect_match(manifest$outputs[[1L]]$sha256, "^[0-9a-f]{64}$")
})

testthat::test_that("manifest metadata rejects empty fields and fractional dimensions", {
  output_path <- tempfile(fileext = ".bed.gz")
  writeBin(charToRaw("bed-output"), output_path)
  valid <- tibble::tibble(
    logical_name = "a_expression",
    path = output_path,
    n_genes = 2,
    n_samples = 3,
    scale = "log2_cpm",
    cell_group = "A"
  )

  for (field in c("logical_name", "scale", "cell_group")) {
    invalid <- valid
    invalid[[field]] <- ""
    testthat::expect_error(validate_manifest_outputs(invalid), "non-empty")
  }
  invalid_dimensions <- valid
  invalid_dimensions$n_genes <- 2.5
  testthat::expect_error(
    validate_manifest_outputs(invalid_dimensions),
    "integer"
  )
})

testthat::test_that("manifest provenance accepts digests and approved image tags", {
  testthat::expect_true(exists("validate_container_image", mode = "function"))
  if (!exists("validate_container_image", mode = "function")) {
    return(invisible(NULL))
  }
  valid_digest <- paste0(
    "ghcr.io/example/pipeline@sha256:",
    paste(rep("0123456789abcdef", 4L), collapse = "")
  )
  testthat::expect_identical(validate_container_image(valid_digest), valid_digest)
  testthat::expect_identical(
    validate_container_image("celltype-deconvolution:test"),
    "celltype-deconvolution:test"
  )
  testthat::expect_identical(
    validate_container_image(
      "ghcr.io/aou-multiomics-analysis/celltypedeconvolution:latest"
    ),
    "ghcr.io/aou-multiomics-analysis/celltypedeconvolution:latest"
  )
  purrr::walk(
    c(
      "ghcr.io/example/pipeline:latest",
      "ghcr.io/example/pipeline@sha256:abc",
      paste0("ghcr.io/example/pipeline@sha256:", paste(rep("A", 64L), collapse = "")),
      "example:test"
    ),
    ~ testthat::expect_error(validate_container_image(.x), "immutable.*SHA-256")
  )
})

testthat::test_that("constant exclusions are derived from input and modeled genes", {
  testthat::expect_true(exists(
    "count_excluded_constant_genes",
    mode = "function"
  ))
  if (!exists("count_excluded_constant_genes", mode = "function")) {
    return(invisible(NULL))
  }
  coordinates <- tibble::tibble(
    `#chr` = rep("chr1", 3L),
    start = c(0L, 10L, 20L),
    end = c(5L, 15L, 25L),
    gene_id = c("g1", "g2", "g3")
  )

  testthat::expect_identical(
    count_excluded_constant_genes(coordinates, c("g1", "g3")),
    1L
  )
  testthat::expect_error(
    count_excluded_constant_genes(coordinates, c("g1", "missing")),
    "modeled gene"
  )
})

testthat::test_that("pipeline QC records validation, row sums, and convergence", {
  testthat::expect_true(exists("build_pipeline_qc_summary", mode = "function"))
  if (!exists("build_pipeline_qc_summary", mode = "function")) {
    return(invisible(NULL))
  }
  export_summary <- tibble::tibble(
    metric = c("gene_count", "sample_count"),
    value = c(3, 2)
  )
  original <- matrix(
    rep(1 / 22, 44L),
    nrow = 2L,
    dimnames = list(c("S1", "S2"), lm22_cell_types())
  )
  original[1L, 1L] <- original[1L, 1L] + 5e-7
  combined <- matrix(
    c(0.6, 0.4, 0.5, 0.5),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(c("S1", "S2"), c("A", "B"))
  )
  weights <- matrix(
    c(0.600001, 0.399999, 0.5, 0.5),
    nrow = 2L,
    byrow = TRUE,
    dimnames = dimnames(combined)
  )
  filter_report <- tibble::tibble(
    cell_group = c("A", "B"),
    retained = c(TRUE, TRUE),
    zero_count_before = c(1, 0),
    zero_floor = c(1e-6, 1e-6)
  )
  model <- list(tau_hat = 0.25)
  tca_log <- c(
    "INFO Iteration 1 out of 10 internal iterations...",
    "INFO Iteration 2 out of 10 internal iterations...",
    "INFO Internal loop converged."
  )
  dtangle_metadata <- list(
    lm22_qc = list(
      gene_count = 66L,
      cell_type_count = 22L,
      value_min = 0.25,
      value_max = 61.25,
      validation_status = "passed"
    )
  )

  summary <- build_pipeline_qc_summary(
    export_summary = export_summary,
    original_proportions = original,
    combined_proportions = combined,
    tca_weights = weights,
    filter_report = filter_report,
    tca_model = model,
    tca_log_lines = tca_log,
    dtangle_metadata = dtangle_metadata
  )
  expected_metric_prefix <- c(
    "gene_count",
    "sample_count",
    "lm22_gene_count",
    "lm22_cell_type_count",
    "lm22_value_min",
    "lm22_value_max",
    "lm22_value_validation",
    "input_proportion_max_row_sum_error",
    "combined_proportion_max_row_sum_error",
    "adjusted_weight_max_row_sum_error",
    "normalization_adjustment_max_abs",
    "zero_values_adjusted"
  )
  expected_metric_suffix <- c(
    "tca_internal_iterations",
    "tca_max_internal_iterations",
    "tca_convergence",
    "tca_tau_hat"
  )
  testthat::expect_identical(
    summary$metric,
    c(expected_metric_prefix, expected_metric_suffix)
  )
  metric_value <- stats::setNames(summary$value, summary$metric)
  metric_status <- stats::setNames(summary$status, summary$metric)

  testthat::expect_equal(metric_value[["lm22_value_min"]], 0.25)
  testthat::expect_equal(metric_value[["lm22_value_max"]], 61.25)
  testthat::expect_identical(metric_status[["lm22_value_validation"]], "passed")
  testthat::expect_equal(
    metric_value[["input_proportion_max_row_sum_error"]],
    5e-7,
    tolerance = 1e-12
  )
  testthat::expect_equal(metric_value[["combined_proportion_max_row_sum_error"]], 0)
  testthat::expect_equal(metric_value[["adjusted_weight_max_row_sum_error"]], 0)
  testthat::expect_equal(metric_value[["normalization_adjustment_max_abs"]], 1e-6)
  testthat::expect_equal(metric_value[["zero_values_adjusted"]], 1)
  testthat::expect_equal(metric_value[["tca_internal_iterations"]], 2)
  testthat::expect_equal(metric_value[["tca_max_internal_iterations"]], 10)
  testthat::expect_identical(metric_status[["tca_convergence"]], "converged")
  testthat::expect_equal(metric_value[["tca_tau_hat"]], 0.25)

})

testthat::test_that("manifest CLI hashes localized files and publishes basenames", {
  working_directory <- tempfile("manifest cli inputs ")
  dir.create(working_directory)
  bed_path <- file.path(working_directory, "a group.bed.gz")
  writeBin(charToRaw("localized-bed-content"), bed_path)
  inventory_path <- file.path(working_directory, "checksum inventory.tsv")
  readr::write_tsv(tibble::tibble(
    logical_name = "a_group_expression",
    path = bed_path,
    n_genes = 2L,
    n_samples = 2L,
    scale = "log2_cpm",
    cell_group = "A group"
  ), inventory_path)
  export_qc_path <- file.path(working_directory, "export qc.tsv")
  readr::write_tsv(tibble::tibble(
    metric = c("gene_count", "sample_count"),
    value = c(2, 2)
  ), export_qc_path)
  original <- matrix(
    rep(1 / 22, 44L),
    nrow = 2L,
    dimnames = list(c("S1", "S2"), lm22_cell_types())
  )
  combined <- matrix(
    c(0.6, 0.4, 0.5, 0.5),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(c("S1", "S2"), c("A", "B"))
  )
  original_path <- file.path(working_directory, "original proportions.tsv")
  combined_path <- file.path(working_directory, "combined proportions.tsv")
  weights_path <- file.path(working_directory, "tca weights.tsv")
  write_numeric_matrix(original, original_path, "sample_id")
  write_numeric_matrix(combined, combined_path, "sample_id")
  write_numeric_matrix(combined, weights_path, "sample_id")
  filter_path <- file.path(working_directory, "filter report.tsv")
  readr::write_tsv(tibble::tibble(
    cell_group = c("A", "B"),
    cohort_mean = c(0.55, 0.45),
    threshold = c(0.0001, 0.0001),
    retained = c(TRUE, TRUE),
    filter_reason = c("retained", "retained"),
    zero_count_before = c(0L, 0L),
    zero_floor = c(1e-6, 1e-6)
  ), filter_path)
  model_path <- file.path(working_directory, "model.rds")
  saveRDS(list(tau_hat = 0.25), model_path)
  model_log_path <- file.path(working_directory, "model log.txt")
  writeLines(c(
    "Iteration 1 out of 10 internal iterations...",
    "Internal loop converged."
  ), model_log_path)
  manifest_path <- file.path(working_directory, "output manifest.json")
  qc_path <- file.path(working_directory, "final qc.tsv")
  log_path <- file.path(working_directory, "manifest log.txt")
  script <- testthat::test_path("../..", "scripts", "build_manifest.R")
  arguments <- c(
    script,
    "--outputs", inventory_path,
    "--export-qc-summary", export_qc_path,
    "--original-proportions", original_path,
    "--combined-proportions", combined_path,
    "--tca-weights", weights_path,
    "--filter-report", filter_path,
    "--model", model_path,
    "--model-log", model_log_path,
    "--tca-version", "1.2.1",
    "--container-image", "celltype-deconvolution:test",
    "--output", manifest_path,
    "--qc-output", qc_path,
    "--log-file", log_path
  )

  status <- system2("Rscript", shQuote(arguments))

  testthat::expect_identical(status, 0L)
  manifest <- jsonlite::read_json(manifest_path, simplifyVector = FALSE)
  testthat::expect_identical(
    manifest$outputs[[1L]]$path,
    basename(bed_path)
  )
  testthat::expect_identical(
    manifest$outputs[[1L]]$sha256,
    digest::digest(file = bed_path, algo = "sha256", serialize = FALSE)
  )
  final_qc <- readr::read_tsv(qc_path, show_col_types = FALSE)
  testthat::expect_true("tca_convergence" %in% final_qc$metric)
  testthat::expect_false("pipeline_version" %in% names(manifest))
  testthat::expect_false(
    "duplicate_gene_symbol_input_row_count" %in% final_qc$metric
  )
  testthat::expect_false("duplicate_gene_symbol_count" %in% final_qc$metric)
})

testthat::test_that("QC plots use a minimal theme without titles", {
  weights <- matrix(
    c(0.2, 0.8, 0.7, 0.3), nrow = 2, byrow = TRUE,
    dimnames = list(c("S1", "S2"), c("A", "B"))
  )
  plots <- make_qc_plots(weights, c(1, 2, 3, 4), c(1.1, 1.9, 3.2, 3.8))

  testthat::expect_length(plots, 3L)
  testthat::expect_true(all(purrr::map_lgl(
    plots,
    ~ is.null(.x$labels$title) && is.null(.x$labels$subtitle)
  )))
})
