source(testthat::test_path("helper-load.R"), local = .GlobalEnv)

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
  testthat::expect_true(all(result$inventory$scale == "log2_cpm"))
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

testthat::test_that("manifest records hashes, dimensions, and provenance", {
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
    pipeline_version = "test",
    tca_version = "1.2.1",
    parameters = list(scale = "log2_cpm"),
    container_image = "example.org/pipeline@sha256:abc"
  )

  testthat::expect_identical(manifest$pipeline_version, "test")
  testthat::expect_identical(manifest$tca_version, "1.2.1")
  testthat::expect_identical(manifest$outputs[[1L]]$dimensions, c(2L, 3L))
  testthat::expect_match(manifest$outputs[[1L]]$sha256, "^[0-9a-f]{64}$")
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
