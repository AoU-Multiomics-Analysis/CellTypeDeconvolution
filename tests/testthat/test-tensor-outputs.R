source(testthat::test_path("helper-load.R"), local = .GlobalEnv)

tensor_outputs_path <- testthat::test_path("../..", "R", "tensor_outputs.R")
if (file.exists(tensor_outputs_path)) {
  source(tensor_outputs_path, local = .GlobalEnv)
}
qc_path <- testthat::test_path("../..", "R", "qc.R")
if (file.exists(qc_path)) {
  source(qc_path, local = .GlobalEnv)
}

make_shard_manifest <- function(genes, shard_ids) {
  tibble::tibble(
    gene_index = seq_along(genes),
    gene_name = genes,
    shard_id = as.integer(shard_ids),
    shard_name = sprintf("shard_%05d", shard_ids),
    index_within_shard = ave(
      seq_along(genes), shard_ids, FUN = seq_along
    ) |> as.integer()
  )
}

testthat::test_that("assembled HDF5 uses the log2 CPM contract", {
  testthat::skip_if_not_installed("hdf5r")
  samples <- c("S1", "S2")
  shard <- list(
    "B cells" = matrix(
      1:4,
      nrow = 2,
      dimnames = list(c("G1", "G2"), samples)
    )
  )
  shard_path <- tempfile(fileext = ".h5")
  write_tensor_shard(shard_path, shard, 1L)
  shard_h5 <- hdf5r::H5File$new(shard_path, mode = "r")
  testthat::expect_true(
    "expression" %in% names(shard_h5[["sources/b_cells"]])
  )
  shard_h5$close_all()
  manifest <- make_shard_manifest(c("G1", "G2"), c(1L, 1L))

  output <- assemble_hdf5_shards(
    shard_path, manifest, tempfile(), "test", "1.2.1"
  )
  h5 <- hdf5r::H5File$new(output[["B cells"]], mode = "r")
  on.exit(h5$close_all())

  testthat::expect_true(all(c(
    "expression", "gene_name", "sample_id"
  ) %in% names(h5)))
  testthat::expect_identical(h5$attr_open("scale")$read(), "log2_cpm")
  testthat::expect_identical(h5$attr_open("cell_group")$read(), "B cells")
  testthat::expect_identical(h5$attr_open("pipeline_version")$read(), "test")
  testthat::expect_identical(h5$attr_open("tca_version")$read(), "1.2.1")
  testthat::expect_identical(h5[["gene_name"]][], c("G1", "G2"))
  testthat::expect_identical(h5[["sample_id"]][], samples)
  testthat::expect_equal(
    unname(h5[["expression"]][, ]),
    unname(shard[[1L]])
  )
  testthat::expect_true(all(h5[["expression"]]$chunk_dims <= c(500L, 256L)))
  creation_properties <- h5[["expression"]]$get_create_plist()
  on.exit(creation_properties$close(), add = TRUE)
  compression <- creation_properties$get_filter(0L)
  testthat::expect_identical(compression$name, "deflate")
  testthat::expect_identical(compression$cd_values, 6L)
})

testthat::test_that("optional TSV assembly streams shard rows in gene order", {
  testthat::skip_if_not_installed("hdf5r")
  samples <- c("S1", "S2")
  shards <- list(
    list("B cells" = matrix(
      1:4,
      nrow = 2,
      dimnames = list(c("G1", "G2"), samples)
    )),
    list("B cells" = matrix(
      5:8,
      nrow = 2,
      dimnames = list(c("G3", "G4"), samples)
    ))
  )
  paths <- c(tempfile(fileext = ".h5"), tempfile(fileext = ".h5"))
  write_tensor_shard(paths[[1L]], shards[[1L]], 1L)
  write_tensor_shard(paths[[2L]], shards[[2L]], 2L)
  manifest <- make_shard_manifest(
    paste0("G", 1:4), c(1L, 1L, 2L, 2L)
  )

  output <- assemble_hdf5_shards(
    paths, manifest, tempfile(), "test", "1.2.1", write_tsv = TRUE
  )
  tsv_path <- attr(output, "tsv_paths")[["B cells"]]
  observed <- readr::read_tsv(
    tsv_path,
    show_col_types = FALSE,
    progress = FALSE
  )

  testthat::expect_identical(observed$gene_name, paste0("G", 1:4))
  testthat::expect_identical(names(observed), c("gene_name", samples))
  testthat::expect_equal(
    unname(as.matrix(dplyr::select(observed, -"gene_name"))),
    rbind(shards[[1L]][[1L]], shards[[2L]][[1L]]) |>
      unname()
  )
})

testthat::test_that("assembly rejects missing shards and duplicated genes", {
  manifest <- make_shard_manifest(
    c("G1", "G1", "G3"), c(1L, 1L, 2L)
  )
  testthat::expect_error(
    assemble_hdf5_shards(
      character(), manifest, tempfile(), "test", "1.2.1"
    ),
    "duplicate.*gene"
  )

  manifest$gene_name[[2L]] <- "G2"
  testthat::expect_error(
    assemble_hdf5_shards(
      character(), manifest, tempfile(), "test", "1.2.1"
    ),
    "missing.*shard"
  )
})

testthat::test_that("shard writing rejects nonfinite values", {
  bad <- list(
    "B cells" = matrix(
      c(1, Inf),
      nrow = 1,
      dimnames = list("G1", c("S1", "S2"))
    )
  )
  testthat::expect_error(
    write_tensor_shard(tempfile(fileext = ".h5"), bad, 1L),
    "finite"
  )
})

testthat::test_that("tensor sources, genes, and samples follow the model", {
  testthat::skip_if_not_installed("TCA")
  testthat::expect_identical(
    as.character(utils::packageVersion("TCA")), "1.2.1"
  )
  set.seed(20260901)
  data <- TCA::test_data(12, 16, 3, 0, 0, 0.01)
  model <- TCA::tca(
    data$X,
    data$W,
    refit_W = FALSE,
    max_iters = 2,
    verbose = FALSE
  )
  genes <- rownames(data$X)[1:4]

  tensor <- extract_tensor_shard(
    data$X, model, genes, 1L, tempfile()
  )

  testthat::expect_identical(names(tensor), colnames(data$W))
  testthat::expect_true(all(purrr::map_lgl(
    tensor,
    ~ identical(rownames(.x), genes)
  )))
  testthat::expect_true(all(purrr::map_lgl(
    tensor,
    ~ identical(colnames(.x), colnames(data$X))
  )))
})

testthat::test_that("assembly rejects shard sample-order mismatches", {
  testthat::skip_if_not_installed("hdf5r")
  first <- list(
    "B cells" = matrix(
      1:4,
      nrow = 2,
      dimnames = list(c("G1", "G2"), c("S1", "S2"))
    )
  )
  second <- list(
    "B cells" = matrix(
      5:8,
      nrow = 2,
      dimnames = list(c("G3", "G4"), c("S2", "S1"))
    )
  )
  paths <- c(tempfile(fileext = ".h5"), tempfile(fileext = ".h5"))
  write_tensor_shard(paths[[1L]], first, 1L)
  write_tensor_shard(paths[[2L]], second, 2L)
  manifest <- make_shard_manifest(
    paste0("G", 1:4), c(1L, 1L, 2L, 2L)
  )

  testthat::expect_error(
    assemble_hdf5_shards(
      paths, manifest, tempfile(), "test", "1.2.1"
    ),
    "sample order"
  )
})

testthat::test_that("assembly rejects misordered shard files", {
  testthat::skip_if_not_installed("hdf5r")
  samples <- c("S1", "S2")
  paths <- c(tempfile(fileext = ".h5"), tempfile(fileext = ".h5"))
  write_tensor_shard(paths[[1L]], list(
    "B cells" = matrix(
      1:4,
      nrow = 2,
      dimnames = list(c("G1", "G2"), samples)
    )
  ), 1L)
  write_tensor_shard(paths[[2L]], list(
    "B cells" = matrix(
      5:8,
      nrow = 2,
      dimnames = list(c("G3", "G4"), samples)
    )
  ), 2L)
  manifest <- make_shard_manifest(
    paste0("G", 1:4), c(1L, 1L, 2L, 2L)
  )

  testthat::expect_error(
    assemble_hdf5_shards(
      rev(paths), manifest, tempfile(), "test", "1.2.1"
    ),
    "shard order"
  )
})

testthat::test_that("assembly rejects shard dimension mismatches", {
  testthat::skip_if_not_installed("hdf5r")
  first <- list(
    "B cells" = matrix(
      1:4,
      nrow = 2,
      dimnames = list(c("G1", "G2"), c("S1", "S2"))
    )
  )
  second <- list(
    "B cells" = matrix(
      5:6,
      nrow = 1,
      dimnames = list("G3", c("S1", "S2"))
    )
  )
  paths <- c(tempfile(fileext = ".h5"), tempfile(fileext = ".h5"))
  write_tensor_shard(paths[[1L]], first, 1L)
  write_tensor_shard(paths[[2L]], second, 2L)
  manifest <- make_shard_manifest(
    paste0("G", 1:4), c(1L, 1L, 2L, 2L)
  )

  testthat::expect_error(
    assemble_hdf5_shards(
      paths, manifest, tempfile(), "test", "1.2.1"
    ),
    "dimension"
  )
})

testthat::test_that("assembly rejects genes that do not follow manifest order", {
  testthat::skip_if_not_installed("hdf5r")
  shard <- list(
    "B cells" = matrix(
      1:4,
      nrow = 2,
      dimnames = list(c("G2", "G1"), c("S1", "S2"))
    )
  )
  path <- tempfile(fileext = ".h5")
  write_tensor_shard(path, shard, 1L)
  manifest <- make_shard_manifest(c("G1", "G2"), c(1L, 1L))

  testthat::expect_error(
    assemble_hdf5_shards(
      path, manifest, tempfile(), "test", "1.2.1"
    ),
    "gene order"
  )
})

testthat::test_that("reconstruction includes the fitted C2 term", {
  tensor <- list(
    A = matrix(
      c(1, 3, 2, 4),
      nrow = 2,
      dimnames = list(c("G1", "G2"), c("S1", "S2"))
    ),
    B = matrix(
      c(5, 7, 6, 8),
      nrow = 2,
      dimnames = list(c("G1", "G2"), c("S1", "S2"))
    )
  )
  weights <- matrix(
    c(0.25, 0.75, 0.60, 0.40),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("S1", "S2"), c("A", "B"))
  )
  C2 <- matrix(
    c(2, 3),
    ncol = 1,
    dimnames = list(c("S1", "S2"), "batch")
  )
  deltas_hat <- matrix(
    c(0.5, -1),
    ncol = 1,
    dimnames = list(c("G1", "G2"), "batch")
  )
  expected <- tensor$A |>
    sweep(2L, weights[, "A"], "*")
  expected <- expected + sweep(tensor$B, 2L, weights[, "B"], "*")
  expected <- expected + t(C2 %*% t(deltas_hat))

  observed <- reconstruct_tensor_shard(
    tensor, weights, C2, deltas_hat
  )

  testthat::expect_equal(observed, expected)
})

testthat::test_that("reconstruction statistics combine shards by sample", {
  observed_1 <- matrix(
    c(1, 2, 2, 4),
    nrow = 2,
    dimnames = list(c("G1", "G2"), c("S1", "S2"))
  )
  fitted_1 <- observed_1 + matrix(c(0, 1, 0, -1), nrow = 2)
  observed_2 <- matrix(
    c(3, 5, 6, 10),
    nrow = 2,
    dimnames = list(c("G3", "G4"), c("S1", "S2"))
  )
  fitted_2 <- observed_2 + matrix(c(1, 0, -1, 0), nrow = 2)
  stats <- initialize_reconstruction_stats(c("S1", "S2")) |>
    update_reconstruction_stats(observed_1, fitted_1) |>
    update_reconstruction_stats(observed_2, fitted_2)

  metrics <- finalize_reconstruction_stats(stats)
  expected_correlations <- c(
    stats::cor(c(1, 2, 3, 5), c(1, 3, 4, 5)),
    stats::cor(c(2, 4, 6, 10), c(2, 3, 5, 10))
  )

  testthat::expect_identical(metrics$sample_id, c("S1", "S2"))
  testthat::expect_equal(metrics$correlation, expected_correlations)
  testthat::expect_equal(metrics$rmse, c(sqrt(2 / 4), sqrt(2 / 4)))
  testthat::expect_true(all(is.finite(metrics$correlation)))
  testthat::expect_true(all(is.finite(metrics$rmse)))
})

testthat::test_that("manifest records hashes, dimensions, and provenance", {
  output_path <- tempfile(fileext = ".h5")
  writeBin(charToRaw("tensor-output"), output_path)
  outputs <- tibble::tibble(
    logical_name = "b_cells_expression",
    path = output_path,
    n_genes = 2L,
    n_samples = 3L,
    scale = "log2_cpm",
    cell_group = "B cells"
  )

  manifest <- build_output_manifest(
    outputs = outputs,
    pipeline_version = "test",
    tca_version = "1.2.1",
    parameters = list(shard_size = 500L, compression = 6L),
    container_image = "example.org/pipeline@sha256:abc"
  )

  testthat::expect_identical(manifest$pipeline_version, "test")
  testthat::expect_identical(manifest$tca_version, "1.2.1")
  testthat::expect_identical(
    manifest$container_image, "example.org/pipeline@sha256:abc"
  )
  testthat::expect_identical(manifest$outputs[[1L]]$dimensions, c(2L, 3L))
  testthat::expect_identical(manifest$outputs[[1L]]$scale, "log2_cpm")
  testthat::expect_match(manifest$outputs[[1L]]$sha256, "^[0-9a-f]{64}$")
})

testthat::test_that("QC plots use a minimal theme without titles", {
  weights <- matrix(
    c(0.2, 0.8, 0.7, 0.3),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("S1", "S2"), c("A", "B"))
  )
  observed <- c(1, 2, 3, 4)
  reconstructed <- c(1.1, 1.9, 3.2, 3.8)

  plots <- make_qc_plots(weights, observed, reconstructed)

  testthat::expect_length(plots, 3L)
  testthat::expect_true(all(purrr::map_lgl(
    plots,
    ~ is.null(.x$labels$title) && is.null(.x$labels$subtitle)
  )))
  testthat::expect_true(all(purrr::map_lgl(
    plots,
    ~ inherits(.x$theme, "theme")
  )))
})

testthat::test_that("tensor output CLIs log contextual failures", {
  original_working_directory <- setwd(pipeline_root)
  on.exit(setwd(original_working_directory), add = TRUE)
  scripts <- c(
    extract = "scripts/extract_tca_shard.R",
    assemble = "scripts/assemble_tca_outputs.R",
    manifest = "scripts/build_manifest.R"
  )

  purrr::iwalk(scripts, function(script, stage) {
    output <- suppressWarnings(system2(
      file.path(R.home("bin"), "Rscript"),
      script,
      stdout = TRUE,
      stderr = TRUE
    ))
    status <- attr(output, "status")
    testthat::expect_true(!is.null(status) && status != 0L)
    testthat::expect_match(
      paste(output, collapse = "\n"),
      sprintf("stage=%s status=failed", stage)
    )
  })
})

testthat::test_that("assembly and manifest CLIs write declared outputs", {
  testthat::skip_if_not_installed("hdf5r")
  samples <- c("S1", "S2", "S3")
  genes <- paste0("G", 1:4)
  X <- matrix(
    seq_len(12L),
    nrow = 4,
    dimnames = list(genes, samples)
  )
  weights <- matrix(
    c(0.25, 0.75, 0.50, 0.50, 0.75, 0.25),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(samples, c("A", "B"))
  )
  manifest <- make_shard_manifest(genes, c(1L, 1L, 2L, 2L))
  input_dir <- tempfile("assembly-inputs-")
  output_dir <- tempfile("assembly-outputs-")
  dir.create(input_dir)
  dir.create(output_dir)
  expression_path <- file.path(input_dir, "expression.tsv.gz")
  weights_path <- file.path(input_dir, "weights.tsv")
  manifest_path <- file.path(input_dir, "manifest.tsv")
  model_path <- file.path(input_dir, "model.rds")
  shard_paths <- file.path(input_dir, c("shard_1.h5", "shard_2.h5"))
  shard_list_path <- file.path(input_dir, "shards.txt")
  write_numeric_matrix(X, expression_path, "gene_name")
  write_numeric_matrix(weights, weights_path, "sample_id")
  readr::write_tsv(manifest, manifest_path)
  model <- list(
    W = weights,
    C2 = matrix(
      numeric(),
      nrow = length(samples),
      ncol = 0L,
      dimnames = list(samples, character())
    ),
    deltas_hat = matrix(
      numeric(),
      nrow = length(genes),
      ncol = 0L,
      dimnames = list(genes, character())
    )
  )
  saveRDS(model, model_path)
  write_tensor_shard(shard_paths[[1L]], list(
    A = X[1:2, , drop = FALSE],
    B = X[1:2, , drop = FALSE]
  ), 1L)
  write_tensor_shard(shard_paths[[2L]], list(
    A = X[3:4, , drop = FALSE],
    B = X[3:4, , drop = FALSE]
  ), 2L)
  writeLines(shard_paths, shard_list_path)
  original_working_directory <- setwd(pipeline_root)
  on.exit(setwd(original_working_directory), add = TRUE)

  assembly_status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(
      "scripts/assemble_tca_outputs.R",
      "--shard-list", shard_list_path,
      "--manifest", manifest_path,
      "--expression-log", expression_path,
      "--model", model_path,
      "--weights", weights_path,
      "--pipeline-version", "test",
      "--write-tsv",
      "--output-dir", output_dir
    )
  )
  expected_assembly_files <- c(
    "a.h5", "b.h5", "a.tsv.gz", "b.tsv.gz",
    "reconstruction_by_sample.tsv", "qc_summary.tsv", "qc_plots.pdf",
    "assembled_outputs.tsv", "assemble_tca_outputs.log"
  )

  testthat::expect_equal(assembly_status, 0L)
  testthat::expect_true(all(file.exists(file.path(
    output_dir,
    expected_assembly_files
  ))))
  reconstruction <- readr::read_tsv(
    file.path(output_dir, "reconstruction_by_sample.tsv"),
    show_col_types = FALSE
  )
  testthat::expect_equal(reconstruction$correlation, rep(1, 3))
  testthat::expect_equal(reconstruction$rmse, rep(0, 3))

  output_manifest_path <- file.path(output_dir, "output_manifest.json")
  manifest_status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(
      "scripts/build_manifest.R",
      "--outputs", file.path(output_dir, "assembled_outputs.tsv"),
      "--pipeline-version", "test",
      "--container-image", "example.org/pipeline@sha256:abc",
      "--output", output_manifest_path
    )
  )
  output_manifest <- jsonlite::read_json(
    output_manifest_path,
    simplifyVector = FALSE
  )

  testthat::expect_equal(manifest_status, 0L)
  testthat::expect_identical(output_manifest$pipeline_version, "test")
  testthat::expect_length(output_manifest$outputs, 4L)
})
