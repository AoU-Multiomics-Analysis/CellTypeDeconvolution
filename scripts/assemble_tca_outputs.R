#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "bootstrap.R"))

assemble_log_path <- NULL

read_shard_paths <- function(shard_list, shards) {
  paths <- if (!is.null(shard_list) && nzchar(shard_list)) {
    readLines(shard_list, warn = FALSE)
  } else if (!is.null(shards) && nzchar(shards)) {
    strsplit(shards, ",", fixed = TRUE)[[1L]]
  } else {
    character()
  }
  paths <- trimws(paths)
  paths[nzchar(paths)]
}

sample_qc_values <- function(observed, reconstructed, maximum = 5000L) {
  observed_values <- as.vector(observed)
  reconstructed_values <- as.vector(reconstructed)
  keep_count <- min(length(observed_values), maximum)
  keep <- unique(as.integer(round(seq(
    1,
    length(observed_values),
    length.out = keep_count
  ))))
  tibble::tibble(
    observed = observed_values[keep],
    reconstructed = reconstructed_values[keep]
  )
}

run_assemble_tca_outputs <- function() {
  message(sprintf("stage=assemble utc_start=%s", tensor_utc_time()))
  option_list <- list(
    optparse::make_option(
      "--shard-list",
      dest = "shard_list",
      type = "character",
      default = NULL,
      help = "File with one tensor shard HDF5 path per line."
    ),
    optparse::make_option(
      "--shards",
      type = "character",
      default = NULL,
      help = "Comma-separated tensor shard HDF5 paths."
    ),
    optparse::make_option(
      "--manifest",
      type = "character",
      help = "Gene shard manifest TSV."
    ),
    optparse::make_option(
      "--expression-log",
      dest = "expression_log",
      type = "character",
      help = "TCA gene-by-sample log2(CPM) TSV."
    ),
    optparse::make_option(
      "--model",
      type = "character",
      help = "Cohort-wide fitted TCA model RDS."
    ),
    optparse::make_option(
      "--weights",
      type = "character",
      help = "Sample-by-group TCA weight TSV."
    ),
    optparse::make_option(
      "--covariates",
      type = "character",
      default = NULL,
      help = "Optional sample-by-covariate C2 TSV."
    ),
    optparse::make_option(
      "--pipeline-version",
      dest = "pipeline_version",
      type = "character",
      help = "Pipeline version for HDF5 attributes."
    ),
    optparse::make_option(
      "--tca-version",
      dest = "tca_version",
      type = "character",
      default = "1.2.1",
      help = "TCA version for HDF5 attributes."
    ),
    optparse::make_option(
      "--write-tsv",
      dest = "write_tsv",
      action = "store_true",
      default = FALSE,
      help = "Stream compressed TSV copies while assembling HDF5 files."
    ),
    optparse::make_option(
      "--output-dir",
      dest = "output_dir",
      type = "character",
      help = "Directory for assembled matrices and QC files."
    ),
    optparse::make_option(
      "--log-file",
      dest = "log_file",
      type = "character",
      default = NULL,
      help = "Assembly log path."
    )
  )
  options <- optparse::parse_args(optparse::OptionParser(
    option_list = option_list
  ))
  required_options <- c(
    "manifest", "expression_log", "model", "weights",
    "pipeline_version", "output_dir"
  )
  missing_options <- required_options[vapply(
    options[required_options],
    function(value) is.null(value) || !nzchar(value),
    logical(1)
  )]
  if (length(missing_options) > 0L) {
    stop(
      sprintf(
        "Missing required options: %s",
        paste(missing_options, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  shard_paths <- read_shard_paths(options$shard_list, options$shards)
  if (length(shard_paths) == 0L) {
    stop("At least one tensor shard path is required", call. = FALSE)
  }

  dir.create(options$output_dir, recursive = TRUE, showWarnings = FALSE)
  assemble_log_path <<- if (is.null(options$log_file)) {
    file.path(options$output_dir, "assemble_tca_outputs.log")
  } else {
    options$log_file
  }
  manifest <- readr::read_tsv(
    options$manifest,
    col_types = readr::cols(
      gene_index = readr::col_integer(),
      gene_name = readr::col_character(),
      shard_id = readr::col_integer(),
      shard_name = readr::col_character(),
      index_within_shard = readr::col_integer()
    ),
    show_col_types = FALSE,
    progress = FALSE
  ) |>
    validate_shard_manifest()
  X <- read_numeric_matrix(options$expression_log, "gene_name")
  weights <- read_numeric_matrix(options$weights, "sample_id")
  model <- readRDS(options$model)
  C2 <- if (is.null(options$covariates) || !nzchar(options$covariates)) {
    model$C2
  } else {
    read_numeric_matrix(options$covariates, "sample_id")
  }
  if (!identical(rownames(X), manifest$gene_name)) {
    stop("Expression gene order must match the shard manifest exactly",
      call. = FALSE
    )
  }
  if (!identical(colnames(X), rownames(weights))) {
    stop("Weight sample order must match expression sample order exactly",
      call. = FALSE
    )
  }
  if (!identical(colnames(weights), colnames(model$W))) {
    stop("Weight source order must match the TCA model exactly", call. = FALSE)
  }
  if (!is.null(C2) && ncol(C2) > 0L &&
      !identical(rownames(C2), colnames(X))) {
    stop("C2 sample order must match expression sample order exactly",
      call. = FALSE
    )
  }
  dimensions_message <- sprintf(
    paste0(
      "stage=assemble input_dimensions=genes:%d samples:%d sources:%d ",
      "shards:%d covariates:%d scale=log2_cpm"
    ),
    nrow(X),
    ncol(X),
    ncol(weights),
    length(shard_paths),
    if (is.null(C2)) 0L else ncol(C2)
  )
  paths_message <- sprintf(
    "stage=assemble input_paths=%s output_paths=%s",
    paste(
      c(
        shard_paths, options$manifest, options$expression_log,
        options$model, options$weights, options$covariates
      ),
      collapse = ","
    ),
    options$output_dir
  )
  message(dimensions_message)
  message(paths_message)
  append_tensor_log(assemble_log_path, dimensions_message)
  append_tensor_log(assemble_log_path, paths_message)

  hdf5_paths <- assemble_hdf5_shards(
    shard_paths = shard_paths,
    manifest = manifest,
    output_dir = options$output_dir,
    pipeline_version = options$pipeline_version,
    tca_version = options$tca_version,
    write_tsv = options$write_tsv
  )
  tsv_paths <- attr(hdf5_paths, "tsv_paths")
  expected_shard_ids <- sort(unique(manifest$shard_id))
  statistics <- initialize_reconstruction_stats(colnames(X))
  qc_points <- list()
  for (index in seq_along(shard_paths)) {
    shard <- read_tensor_shard(shard_paths[[index]])
    shard_id <- shard$shard_id
    if (!identical(shard_id, expected_shard_ids[[index]])) {
      stop("Tensor shard order changed during reconstruction", call. = FALSE)
    }
    genes <- manifest$gene_name[manifest$shard_id == shard_id]
    C2_shard <- if (is.null(C2) || ncol(C2) == 0L) NULL else C2
    deltas_shard <- if (is.null(C2_shard)) {
      NULL
    } else {
      model$deltas_hat[genes, colnames(C2_shard), drop = FALSE]
    }
    reconstructed <- reconstruct_tensor_shard(
      tensor = shard$tensor,
      weights = weights,
      C2 = C2_shard,
      deltas_hat = deltas_shard
    )
    observed <- X[genes, , drop = FALSE]
    statistics <- update_reconstruction_stats(
      statistics,
      observed,
      reconstructed
    )
    qc_points[[index]] <- sample_qc_values(observed, reconstructed)
  }
  reconstruction_by_sample <- finalize_reconstruction_stats(statistics)
  qc_points <- dplyr::bind_rows(qc_points)
  qc_paths <- write_qc_reports(
    weights = weights,
    reconstruction_by_sample = reconstruction_by_sample,
    observed = qc_points$observed,
    reconstructed = qc_points$reconstructed,
    gene_count = nrow(X),
    output_dir = options$output_dir
  )
  output_inventory <- tibble::tibble(
    logical_name = paste0(slugify_cell_group(names(hdf5_paths)), "_expression"),
    path = unname(hdf5_paths),
    n_genes = nrow(X),
    n_samples = ncol(X),
    scale = "log2_cpm",
    cell_group = names(hdf5_paths)
  )
  if (length(tsv_paths) > 0L) {
    output_inventory <- dplyr::bind_rows(
      output_inventory,
      tibble::tibble(
        logical_name = paste0(
          slugify_cell_group(names(tsv_paths)),
          "_expression_tsv"
        ),
        path = unname(tsv_paths),
        n_genes = nrow(X),
        n_samples = ncol(X),
        scale = "log2_cpm",
        cell_group = names(tsv_paths)
      )
    )
  }
  inventory_path <- file.path(options$output_dir, "assembled_outputs.tsv")
  readr::write_tsv(output_inventory, inventory_path, na = "")
  complete_message <- sprintf(
    paste0(
      "stage=assemble event=stage_complete output_dimensions=genes:%d ",
      "samples:%d groups:%d output_paths=%s scale=log2_cpm"
    ),
    nrow(X),
    ncol(X),
    length(hdf5_paths),
    paste(c(hdf5_paths, tsv_paths, qc_paths, inventory_path), collapse = ",")
  )
  message(sprintf("%s utc_complete=%s", complete_message, tensor_utc_time()))
  append_tensor_log(assemble_log_path, complete_message)
}

tryCatch(
  run_assemble_tca_outputs(),
  error = function(error) {
    error_message <- sprintf(
      "stage=assemble status=failed utc_time=%s message=%s",
      tensor_utc_time(),
      conditionMessage(error)
    )
    message(error_message)
    if (!is.null(assemble_log_path)) {
      append_tensor_log(assemble_log_path, error_message)
    }
    quit(status = 1L)
  }
)
