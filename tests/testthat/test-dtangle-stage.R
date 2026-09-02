source(testthat::test_path("helper-load.R"), local = .GlobalEnv)

dtangle_stage_path <- testthat::test_path("../..", "R", "dtangle_stage.R")
if (file.exists(dtangle_stage_path)) {
  source(dtangle_stage_path, local = .GlobalEnv)
}

make_synthetic_lm22 <- function(markers_per_type = 3L, baseline = 4, marker = 1024) {
  genes <- paste0("G", seq_len(22L * markers_per_type))
  reference <- matrix(
    baseline,
    nrow = length(genes),
    ncol = 22L,
    dimnames = list(genes, lm22_cell_types())
  )
  for (cell_index in seq_len(22L)) {
    first <- (cell_index - 1L) * markers_per_type + 1L
    last <- cell_index * markers_per_type
    reference[first:last, cell_index] <- marker
  }
  reference
}

testthat::test_that("standard LM22 is logged without a pseudocount", {
  lm22 <- matrix(
    rep(c(0.25, 4, 16), each = 22), nrow = 3, byrow = TRUE,
    dimnames = list(c("G1", "G2", "G3"), lm22_cell_types())
  )

  observed <- transform_lm22(lm22)

  testthat::expect_equal(observed, log2(lm22))
})

testthat::test_that("LM22 rejects zero and missing cell types", {
  lm22 <- matrix(
    1,
    nrow = 3,
    ncol = 21,
    dimnames = list(c("G1", "G2", "G3"), lm22_cell_types()[-1])
  )

  testthat::expect_error(validate_lm22(lm22), "22 standard LM22 columns")
})

testthat::test_that("LM22 rejects duplicate trimmed gene symbols", {
  lm22 <- matrix(
    1,
    nrow = 3,
    ncol = 22,
    dimnames = list(c(" G1", "G1 ", "G2"), lm22_cell_types())
  )

  testthat::expect_error(validate_lm22(lm22), "unique")
})

testthat::test_that("LM22 rejects negative values", {
  lm22 <- matrix(
    1,
    nrow = 3,
    ncol = 22,
    dimnames = list(c("G1", "G2", "G3"), lm22_cell_types())
  )
  lm22[1, 1] <- -1

  testthat::expect_error(validate_lm22(lm22), "positive")
})

testthat::test_that("LM22 rejects zero values", {
  lm22 <- matrix(
    1,
    nrow = 3,
    ncol = 22,
    dimnames = list(c("G1", "G2", "G3"), lm22_cell_types())
  )
  lm22[1, 1] <- 0

  testthat::expect_error(validate_lm22(lm22), "positive")
})

testthat::test_that("LM22 rejects nonfinite values", {
  lm22 <- matrix(
    1,
    nrow = 3,
    ncol = 22,
    dimnames = list(c("G1", "G2", "G3"), lm22_cell_types())
  )
  lm22[1, 1] <- Inf

  testthat::expect_error(validate_lm22(lm22), "finite")
})

testthat::test_that("dtangle input preparation stops below the overlap threshold", {
  lm22 <- make_synthetic_lm22()
  bulk_log <- matrix(
    log2(5),
    nrow = 52,
    ncol = 1,
    dimnames = list(rownames(lm22)[seq_len(52)], "S1")
  )

  testthat::expect_error(
    prepare_dtangle_inputs(bulk_log, lm22, min_overlap = 0.80),
    "overlap.*below"
  )
})

testthat::test_that("dtangle input preparation preserves LM22 and official orders", {
  lm22 <- make_synthetic_lm22()
  lm22 <- lm22[, rev(lm22_cell_types()), drop = FALSE]
  bulk_log <- matrix(
    log2(5),
    nrow = nrow(lm22),
    ncol = 2,
    dimnames = list(rev(rownames(lm22)), c("S2", "S1"))
  )

  inputs <- prepare_dtangle_inputs(bulk_log, lm22, min_overlap = 0.80)

  testthat::expect_identical(colnames(inputs$Y), rownames(lm22))
  testthat::expect_identical(rownames(inputs$Y), c("S2", "S1"))
  testthat::expect_identical(rownames(inputs$references), lm22_cell_types())
  testthat::expect_identical(colnames(inputs$references), rownames(lm22))
})

testthat::test_that("joint quantile normalization uses joined log-scale profiles", {
  lm22 <- make_synthetic_lm22()
  bulk_log <- matrix(
    log2(5),
    nrow = nrow(lm22),
    ncol = 2,
    dimnames = list(rownames(lm22), c("S1", "S2"))
  )

  inputs <- prepare_dtangle_inputs(
    bulk_log,
    lm22,
    min_overlap = 0.80,
    quantile_normalize = TRUE
  )
  expected <- limma::normalizeBetweenArrays(cbind(log2(lm22), bulk_log))

  testthat::expect_equal(inputs$references, t(expected[, lm22_cell_types()]))
  testthat::expect_equal(inputs$Y, t(expected[, c("S1", "S2")]))
})

testthat::test_that("dtangle returns normalized proportions and nonempty markers", {
  set.seed(20260901)
  reference <- make_synthetic_lm22(markers_per_type = 3L, baseline = 4, marker = 1024)
  weights <- matrix(
    rexp(8 * 22, rate = 1),
    nrow = 8,
    dimnames = list(paste0("S", 1:8), lm22_cell_types())
  )
  weights <- weights / rowSums(weights)
  bulk_linear <- reference %*% t(weights)
  bulk_log <- log2(bulk_linear + 1)

  inputs <- prepare_dtangle_inputs(bulk_log, reference, 0.80, FALSE)
  fit <- estimate_dtangle(inputs, marker_fraction = 0.10)
  marker_counts <- dplyr::count(fit$markers, .data$cell_type, name = "marker_count")

  testthat::expect_equal(dim(fit$proportions), c(8L, 22L))
  testthat::expect_identical(rownames(fit$proportions), paste0("S", 1:8))
  testthat::expect_identical(colnames(fit$proportions), lm22_cell_types())
  testthat::expect_true(all(fit$proportions >= 0))
  testthat::expect_equal(
    unname(rowSums(fit$proportions)),
    rep(1, 8),
    tolerance = 1e-8
  )
  marker_counts <- marker_counts |>
    dplyr::arrange(match(.data$cell_type, lm22_cell_types()))
  testthat::expect_identical(marker_counts$cell_type, lm22_cell_types())
  testthat::expect_true(all(marker_counts$marker_count > 0L))
})

testthat::test_that("the dtangle CLI writes its six declared outputs", {
  lm22 <- make_synthetic_lm22()
  bulk_log <- matrix(
    log2(5),
    nrow = nrow(lm22),
    ncol = 1,
    dimnames = list(rownames(lm22), "S1")
  )
  lm22_path <- tempfile(fileext = ".tsv")
  bulk_path <- tempfile(fileext = ".tsv")
  output_dir <- tempfile()
  dir.create(output_dir)
  write_numeric_matrix(lm22, lm22_path, "gene_symbol")
  write_numeric_matrix(bulk_log, bulk_path, "gene_symbol")
  original_working_directory <- setwd(pipeline_root)
  on.exit(setwd(original_working_directory), add = TRUE)

  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(
      "scripts/run_dtangle.R", "--bulk-log", bulk_path, "--lm22", lm22_path,
      "--output-dir", output_dir
    )
  )

  testthat::expect_equal(status, 0L)
  testthat::expect_true(all(file.exists(file.path(output_dir, c(
    "dtangle_proportions.tsv", "dtangle_markers.tsv", "dtangle_metadata.json",
    "dtangle_overlap.tsv", "dtangle_lm22_log.tsv.gz", "dtangle_shared_bulk.tsv.gz"
  )))))
})
