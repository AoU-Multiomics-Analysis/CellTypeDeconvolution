#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "bootstrap.R"))

export_log_path <- NULL

run_export_tca_beds <- function() {
  option_list <- list(
    optparse::make_option("--expression", type = "character",
      help = "Coordinate-preserving BED of positive linear CPM values."),
    optparse::make_option("--tca-expression", dest = "tca_expression", type = "character",
      help = "Filtered gene-by-sample log2(CPM) TSV with gene_id first column."),
    optparse::make_option("--model", type = "character",
      help = "Cohort-wide fitted TCA model RDS."),
    optparse::make_option("--weights", type = "character",
      help = "Sample-by-group TCA weight TSV with sample_id first column."),
    optparse::make_option("--covariates", type = "character", default = NULL,
      help = "Optional sample-by-covariate TSV with sample_id first column."),
    optparse::make_option("--num-cores", dest = "num_cores", type = "integer", default = 8L,
      help = "Number of CPU cores for tensor extraction."),
    optparse::make_option("--output-dir", dest = "output_dir", type = "character",
      help = "Output directory for cell-type BED and QC files."),
    optparse::make_option("--log-file", dest = "log_file", type = "character", default = NULL,
      help = "Export detail log path.")
  )
  options <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
  required_options <- c(
    "expression", "tca_expression", "model", "weights", "output_dir"
  )
  missing_options <- required_options[vapply(
    options[required_options],
    function(value) is.null(value) || !nzchar(value),
    logical(1)
  )]
  if (length(missing_options) > 0L) {
    stop(sprintf("Missing required options: %s", paste(missing_options, collapse = ", ")),
      call. = FALSE
    )
  }
  dir.create(options$output_dir, recursive = TRUE, showWarnings = FALSE)
  export_log_path <<- if (is.null(options$log_file) || !nzchar(options$log_file)) {
    file.path(options$output_dir, "export_tca_beds.log")
  } else {
    options$log_file
  }
  append_tensor_log(export_log_path, sprintf(
    "stage=export_tca_beds event=stage_start output_dir=%s",
    normalizePath(options$output_dir)
  ))

  X <- read_numeric_matrix(options$tca_expression, "gene_id")
  expression <- read_expression_bed(options$expression)
  coordinates <- expression$coordinates
  weights <- read_numeric_matrix(options$weights, "sample_id")
  model <- readRDS(options$model)
  C2 <- if (is.null(options$covariates) || !nzchar(options$covariates)) {
    model$C2
  } else {
    read_numeric_matrix(options$covariates, "sample_id")
  }
  if (!is.null(C2) && ncol(C2) == 0L) {
    C2 <- NULL
  }
  if (!identical(colnames(X), rownames(weights))) {
    stop("Weight sample order must match expression sample order exactly", call. = FALSE)
  }
  if (!identical(colnames(weights), colnames(model$W))) {
    stop("Weight group order must match the TCA model exactly", call. = FALSE)
  }
  if (!is.null(C2) && !identical(rownames(C2), colnames(X))) {
    stop("Covariate sample order must match expression sample order exactly", call. = FALSE)
  }
  excluded_constant_genes <- count_excluded_constant_genes(
    coordinates,
    rownames(X)
  )
  coordinate_index <- match(rownames(X), coordinates$gene_id)
  if (anyNA(coordinate_index)) {
    stop("Every modeled gene_id must have BED coordinates", call. = FALSE)
  }
  coordinates <- coordinates[coordinate_index, , drop = FALSE]
  if (!identical(coordinates$gene_id, rownames(X))) {
    stop("Coordinate gene order must match TCA expression exactly", call. = FALSE)
  }
  dimensions_message <- sprintf(
    "stage=export_tca_beds input_dimensions=genes:%d samples:%d groups:%d covariates:%d scale=log2_cpm",
    nrow(X), ncol(X), ncol(weights), if (is.null(C2)) 0L else ncol(C2)
  )
  paths_message <- sprintf(
    "stage=export_tca_beds input_paths=%s output_dir=%s",
    paste(c(options$tca_expression, options$expression, options$model, options$weights,
      options$covariates), collapse = ","), options$output_dir
  )
  message(dimensions_message)
  message(paths_message)
  append_tensor_log(export_log_path, dimensions_message)
  append_tensor_log(export_log_path, paths_message)

  tensor <- extract_full_tensor(X, model, options$num_cores, export_log_path)
  bed_result <- write_cell_type_beds(tensor, coordinates, options$output_dir)
  deltas_hat <- if (is.null(C2)) {
    NULL
  } else {
    model$deltas_hat[rownames(X), colnames(C2), drop = FALSE]
  }
  reconstructed <- reconstruct_tensor(
    tensor,
    weights,
    C2,
    deltas_hat
  )
  statistics <- initialize_reconstruction_stats(colnames(X)) |>
    update_reconstruction_stats(X, reconstructed)
  reconstruction_by_sample <- finalize_reconstruction_stats(statistics)
  qc_points <- sample_qc_values(X, reconstructed)
  qc_paths <- write_qc_reports(
    weights = weights,
    reconstruction_by_sample = reconstruction_by_sample,
    observed = qc_points$observed,
    reconstructed = qc_points$reconstructed,
    gene_count = nrow(X),
    excluded_constant_gene_count = excluded_constant_genes,
    output_dir = options$output_dir
  )
  inventory_path <- file.path(options$output_dir, "cell_type_bed_inventory.tsv")
  readr::write_tsv(bed_result$inventory, inventory_path, na = "")
  output_paths <- c(
    unname(bed_result$paths), inventory_path, unname(qc_paths)
  )
  complete_message <- sprintf(
    "stage=export_tca_beds event=stage_complete output_dimensions=genes:%d samples:%d groups:%d excluded_constant_genes:%d output_paths=%s scale=log2_cpm",
    nrow(X), ncol(X), ncol(weights), excluded_constant_genes,
    paste(output_paths, collapse = ",")
  )
  message(sprintf("%s utc_complete=%s", complete_message, tensor_utc_time()))
  append_tensor_log(export_log_path, complete_message)
}

tryCatch(
  run_export_tca_beds(),
  error = function(error) {
    error_message <- sprintf(
      "stage=export_tca_beds status=failed utc_time=%s message=%s",
      tensor_utc_time(), conditionMessage(error)
    )
    message(error_message)
    if (!is.null(export_log_path)) {
      append_tensor_log(export_log_path, error_message)
    }
    quit(status = 1L)
  }
)
