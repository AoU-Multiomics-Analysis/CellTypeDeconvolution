#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) < 1L || length(arguments) > 2L) {
  stop(
    "Usage: assert_outputs.R OUTPUTS_JSON [FIXTURE_DIRECTORY]",
    call. = FALSE
  )
}
if (!requireNamespace("hdf5r", quietly = TRUE)) {
  stop("The hdf5r package is required for smoke assertions", call. = FALSE)
}

outputs_path <- arguments[[1L]]
fixture_directory <- if (length(arguments) == 2L) {
  arguments[[2L]]
} else {
  "tests/fixtures"
}
outputs <- jsonlite::read_json(outputs_path, simplifyVector = TRUE)
estimated_proportions_key <-
  "cell_type_deconvolution.estimated_proportions"
expected_proportion_mode <- if (is.null(outputs[[estimated_proportions_key]])) {
  "precomputed"
} else {
  "dtangle"
}

output_value <- function(name) {
  key <- paste0("cell_type_deconvolution.", name)
  value <- outputs[[key]]
  if (is.null(value) || length(value) == 0L) {
    stop(sprintf("Required workflow output is absent: %s", key), call. = FALSE)
  }
  value
}

read_matrix_table <- function(name, id_column) {
  table <- readr::read_tsv(
    output_value(name),
    col_types = readr::cols(.default = readr::col_character()),
    name_repair = "minimal",
    show_col_types = FALSE,
    progress = FALSE
  )
  if (!identical(names(table)[[1L]], id_column)) {
    stop(
      sprintf("%s must start with %s", name, id_column),
      call. = FALSE
    )
  }
  values <- table[-1L] |>
    dplyr::mutate(dplyr::across(dplyr::everything(), readr::parse_double)) |>
    as.matrix()
  rownames(values) <- table[[id_column]]
  values
}

expected_samples <- readLines(
  file.path(fixture_directory, "expected_samples.txt"),
  warn = FALSE
)
expected_full_mapped_genes <- readLines(
  file.path(fixture_directory, "expected_genes.txt"),
  warn = FALSE
)
expected_groups <- readLines(
  file.path(fixture_directory, "expected_groups.txt"),
  warn = FALSE
)
signature_genes <- readr::read_tsv(
  file.path(fixture_directory, "synthetic_signature.tsv"),
  col_types = readr::cols(.default = readr::col_character()),
  show_col_types = FALSE,
  progress = FALSE
)$gene_symbol
non_signature_genes <- setdiff(expected_full_mapped_genes, signature_genes)
stopifnot(length(non_signature_genes) == 6L)

proportions <- read_matrix_table("proportions_lm22", "sample_id")
combined <- read_matrix_table("proportions_combined", "sample_id")
tca_weights <- read_matrix_table("tca_weights", "sample_id")
stopifnot(all(is.finite(proportions)), all(proportions >= 0))
stopifnot(max(abs(rowSums(proportions) - 1)) < 1e-8)
stopifnot(all(is.finite(combined)), all(combined >= 0))
stopifnot(max(abs(rowSums(combined) - 1)) < 1e-8)
stopifnot(all(tca_weights > 0), all(is.finite(tca_weights)))
stopifnot(max(abs(rowSums(tca_weights) - 1)) < 1e-8)
stopifnot(identical(expected_samples, rownames(proportions)))
stopifnot(identical(expected_samples, rownames(combined)))
stopifnot(identical(expected_samples, rownames(tca_weights)))
stopifnot(identical(expected_groups, colnames(combined)))
stopifnot(identical(expected_groups, colnames(tca_weights)))

tca_expression <- read_matrix_table("tca_expression", "gene_name")
observed_genes <- rownames(tca_expression)
observed_samples <- colnames(tca_expression)
stopifnot(identical(expected_samples, observed_samples))
stopifnot(identical(expected_full_mapped_genes, observed_genes))
stopifnot(any(non_signature_genes %in% observed_genes))

shard_manifest <- readr::read_tsv(
  output_value("gene_shard_manifest"),
  show_col_types = FALSE,
  progress = FALSE
)
stopifnot(identical(expected_full_mapped_genes, shard_manifest$gene_name))
stopifnot(identical(seq_along(expected_full_mapped_genes), shard_manifest$gene_index))

reconstruction <- readr::read_tsv(
  output_value("reconstruction_by_sample"),
  show_col_types = FALSE,
  progress = FALSE
)
stopifnot(identical(expected_samples, reconstruction$sample_id))
stopifnot(all(reconstruction$gene_count == length(expected_full_mapped_genes)))
stopifnot(all(is.finite(reconstruction$correlation)))
stopifnot(all(is.finite(reconstruction$rmse)))

group_hdf5 <- output_value("group_hdf5")
group_tsv <- output_value("group_tsv")
stopifnot(length(group_hdf5) == length(expected_groups))
stopifnot(length(group_tsv) == length(expected_groups))

hdf5_groups <- purrr::map_chr(group_hdf5, function(path) {
  h5 <- hdf5r::H5File$new(path, mode = "r")
  on.exit(h5$close_all(), add = TRUE)
  stopifnot(all(c("expression", "gene_name", "sample_id") %in% names(h5)))
  genes <- as.character(h5[["gene_name"]][])
  samples <- as.character(h5[["sample_id"]][])
  dimensions <- h5[["expression"]]$dims
  cell_group <- as.character(h5$attr_open("cell_group")$read())
  scale <- as.character(h5$attr_open("scale")$read())
  stopifnot(
    identical(dimensions, c(
      length(expected_full_mapped_genes),
      length(expected_samples)
    )),
    identical(genes, expected_full_mapped_genes),
    identical(samples, expected_samples),
    identical(scale, "log2_cpm"),
    cell_group %in% expected_groups,
    any(non_signature_genes %in% genes)
  )
  cell_group
})
stopifnot(setequal(hdf5_groups, expected_groups))
stopifnot(anyDuplicated(hdf5_groups) == 0L)

purrr::walk(group_tsv, function(path) {
  table <- readr::read_tsv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    name_repair = "minimal",
    show_col_types = FALSE,
    progress = FALSE
  )
  stopifnot(
    identical(table$gene_name, expected_full_mapped_genes),
    identical(names(table)[-1L], expected_samples),
    any(non_signature_genes %in% table$gene_name)
  )
})

inventory <- readr::read_tsv(
  output_value("output_inventory"),
  show_col_types = FALSE,
  progress = FALSE
)
stopifnot(
  nrow(inventory) == 2L * length(expected_groups),
  all(inventory$n_genes == length(expected_full_mapped_genes)),
  all(inventory$n_samples == length(expected_samples)),
  all(inventory$scale == "log2_cpm"),
  setequal(inventory$cell_group, expected_groups)
)

manifest <- jsonlite::read_json(
  output_value("output_manifest"),
  simplifyVector = FALSE
)
parameters <- manifest$parameters
required_parameter_names <- c(
  "proportion_mode", "min_lm22_overlap", "dtangle_marker_fraction",
  "dtangle_marker_method", "dtangle_quantile_normalize",
  "group_mean_threshold", "zero_floor", "tca_shard_size",
  "tca_max_iters", "random_seed", "write_tsv",
  "hdf5_gene_chunk_max", "hdf5_sample_chunk_max",
  "hdf5_gzip_level", "scale"
)
numeric_parameter <- function(name) {
  as.numeric(parameters[[name]])
}
stopifnot(
  setequal(names(parameters), required_parameter_names),
  identical(parameters$proportion_mode, expected_proportion_mode),
  numeric_parameter("min_lm22_overlap") == 0.80,
  numeric_parameter("dtangle_marker_fraction") == 0.10,
  identical(parameters$dtangle_marker_method, "ratio"),
  identical(parameters$dtangle_quantile_normalize, FALSE),
  numeric_parameter("group_mean_threshold") == 0.0001,
  numeric_parameter("zero_floor") == 0.000001,
  numeric_parameter("tca_shard_size") == 24,
  numeric_parameter("tca_max_iters") == 10,
  numeric_parameter("random_seed") == 20260901,
  identical(parameters$write_tsv, TRUE),
  numeric_parameter("hdf5_gene_chunk_max") == 500,
  numeric_parameter("hdf5_sample_chunk_max") == 256,
  numeric_parameter("hdf5_gzip_level") == 6,
  identical(parameters$scale, "log2_cpm"),
  identical(manifest$tca_version, "1.2.1")
)
manifest_outputs <- manifest$outputs
stopifnot(length(manifest_outputs) == nrow(inventory))
purrr::walk(manifest_outputs, function(entry) {
  stopifnot(
    identical(entry$scale, "log2_cpm"),
    identical(unlist(entry$dimensions, use.names = FALSE), c(
      length(expected_full_mapped_genes),
      length(expected_samples)
    )),
    entry$cell_group %in% expected_groups
  )
})

required_file_outputs <- c(
  "prepared_cpm", "prepared_log2_cpm", "mapping_report",
  "prepare_excluded_genes", "prepare_log", "cell_group_filter_report",
  "proportions_log", "tca_model", "tca_model_log", "tca_excluded_genes",
  "fit_tca_log", "qc_summary", "qc_plots", "assembly_log",
  "effective_parameters_file", "manifest_log"
)
purrr::walk(required_file_outputs, function(name) {
  stopifnot(file.exists(output_value(name)))
})
purrr::walk(output_value("tensor_shards"), ~ stopifnot(file.exists(.x)))
purrr::walk(output_value("tensor_shard_logs"), ~ stopifnot(file.exists(.x)))

message(sprintf(
  paste0(
    "Smoke assertions passed: %d genes, %d samples, %d groups, ",
    "scale=log2_cpm"
  ),
  length(expected_full_mapped_genes),
  length(expected_samples),
  length(expected_groups)
))
