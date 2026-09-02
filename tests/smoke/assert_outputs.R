#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) < 1L || length(arguments) > 2L) {
  stop(
    "Usage: assert_outputs.R OUTPUTS_JSON [FIXTURE_DIRECTORY]",
    call. = FALSE
  )
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

tca_expression <- read_matrix_table("tca_expression", "gene_id")
observed_genes <- rownames(tca_expression)
observed_samples <- colnames(tca_expression)
stopifnot(identical(expected_samples, observed_samples))
stopifnot(identical(expected_full_mapped_genes, observed_genes))

reconstruction <- readr::read_tsv(
  output_value("reconstruction_by_sample"),
  show_col_types = FALSE,
  progress = FALSE
)
stopifnot(identical(expected_samples, reconstruction$sample_id))
stopifnot(all(reconstruction$gene_count == length(expected_full_mapped_genes)))
stopifnot(all(is.finite(reconstruction$correlation)))
stopifnot(all(is.finite(reconstruction$rmse)))

cell_type_beds <- output_value("cell_type_beds")
stopifnot(length(cell_type_beds) == length(expected_groups))
purrr::walk(cell_type_beds, function(path) {
  table <- readr::read_tsv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    name_repair = "minimal",
    show_col_types = FALSE,
    progress = FALSE
  )
  stopifnot(
    identical(table$gene_id, expected_full_mapped_genes),
    identical(names(table)[-(1:4)], expected_samples)
  )
})

inventory <- readr::read_tsv(
  output_value("cell_type_bed_inventory"),
  show_col_types = FALSE,
  progress = FALSE
)
stopifnot(
  nrow(inventory) == length(expected_groups),
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
  "group_mean_threshold", "zero_floor", "tca_max_iters",
  "random_seed", "scale"
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
  numeric_parameter("tca_max_iters") == 10,
  numeric_parameter("random_seed") == 20260901,
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
  "prepared_tca_cpm", "prepared_tca_log2_cpm",
  "prepared_dtangle_cpm", "prepared_dtangle_log2_cpm",
  "prepared_coordinates", "mapping_report",
  "prepare_excluded_genes", "prepare_log", "cell_group_filter_report",
  "proportions_log", "tca_model", "tca_model_log", "tca_excluded_genes",
  "fit_tca_log", "qc_summary", "qc_plots", "export_log",
  "export_detail_log", "output_manifest", "manifest_log",
  "effective_parameters_file"
)
purrr::walk(required_file_outputs, function(name) {
  stopifnot(file.exists(output_value(name)))
})

message(sprintf(
  paste0(
    "Smoke assertions passed: %d genes, %d samples, %d groups, ",
    "scale=log2_cpm"
  ),
  length(expected_full_mapped_genes),
  length(expected_samples),
  length(expected_groups)
))
