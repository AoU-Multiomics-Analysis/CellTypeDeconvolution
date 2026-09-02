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
expected_gene_ids <- readLines(
  file.path(fixture_directory, "expected_gene_ids.txt"),
  warn = FALSE
)
expected_coordinates <- readr::read_tsv(
  file.path(fixture_directory, "expected_coordinates.tsv"),
  col_types = readr::cols(
    `#chr` = readr::col_character(),
    start = readr::col_integer(),
    end = readr::col_integer(),
    gene_id = readr::col_character()
  ),
  show_col_types = FALSE,
  progress = FALSE
)
expected_groups <- readLines(
  file.path(fixture_directory, "expected_groups.txt"),
  warn = FALSE
)
proportions <- read_matrix_table("proportions_lm22", "sample_id")
combined <- read_matrix_table("proportions_combined", "sample_id")
tca_weights <- read_matrix_table("tca_weights", "sample_id")
expected_lm22_types <- c(
  "B cells naive", "B cells memory", "Plasma cells", "T cells CD8",
  "T cells CD4 naive", "T cells CD4 memory resting",
  "T cells CD4 memory activated", "T cells follicular helper",
  "T cells regulatory (Tregs)", "T cells gamma delta",
  "NK cells resting", "NK cells activated", "Monocytes",
  "Macrophages M0", "Macrophages M1", "Macrophages M2",
  "Dendritic cells resting", "Dendritic cells activated",
  "Mast cells resting", "Mast cells activated", "Eosinophils", "Neutrophils"
)
stopifnot(all(is.finite(proportions)), all(proportions >= 0))
stopifnot(identical(colnames(proportions), expected_lm22_types))
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
stopifnot(identical(expected_gene_ids, observed_genes))
stopifnot(all(c("ENSGSYN000001", "ENSGDUP000001") %in% observed_genes))
stopifnot(any(grepl("^ENSGEXTRA", observed_genes)))

if (identical(expected_proportion_mode, "dtangle")) {
  shared_bulk <- read_matrix_table("dtangle_shared_bulk", "gene_symbol")
  stopifnot(sum(rownames(shared_bulk) == "SYN_GENE_001") == 1L)
}

reconstruction <- readr::read_tsv(
  output_value("reconstruction_by_sample"),
  show_col_types = FALSE,
  progress = FALSE
)
stopifnot(identical(expected_samples, reconstruction$sample_id))
stopifnot(all(reconstruction$gene_count == length(expected_gene_ids)))
stopifnot(all(is.finite(reconstruction$correlation)))
stopifnot(all(is.finite(reconstruction$rmse)))

cell_type_beds <- output_value("cell_type_beds")
cell_type_bed_paths <- normalizePath(cell_type_beds)
stopifnot(
  length(cell_type_bed_paths) == length(expected_groups),
  all(grepl("[.]bed[.]gz$", cell_type_bed_paths)),
  anyDuplicated(cell_type_bed_paths) == 0L
)
purrr::walk(cell_type_bed_paths, function(path) {
  bed <- readr::read_tsv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    name_repair = "minimal",
    show_col_types = FALSE,
    progress = FALSE
  )
  stopifnot(
    identical(names(bed)[1:4], c("#chr", "start", "end", "gene_id")),
    identical(bed$gene_id, expected_gene_ids),
    identical(names(bed)[-(1:4)], expected_samples),
    identical(bed[["#chr"]], expected_coordinates[["#chr"]]),
    identical(readr::parse_integer(bed$start), expected_coordinates$start),
    identical(readr::parse_integer(bed$end), expected_coordinates$end)
  )
  values <- bed[-(1:4)] |>
    dplyr::mutate(dplyr::across(dplyr::everything(), readr::parse_double)) |>
    as.matrix()
  stopifnot(all(is.finite(values)))
})

inventory <- readr::read_tsv(
  output_value("cell_type_bed_inventory"),
  show_col_types = FALSE,
  progress = FALSE
)
inventory_basenames <- inventory$path
stopifnot(
  nrow(inventory) == length(expected_groups),
  all(inventory$n_genes == length(expected_gene_ids)),
  all(inventory$n_samples == length(expected_samples)),
  all(inventory$scale == "log2_cpm"),
  setequal(inventory$cell_group, expected_groups),
  all(nzchar(inventory$slug)),
  anyDuplicated(inventory$slug) == 0L,
  all(grepl("^[^/]+[.]bed[.]gz$", inventory_basenames)),
  anyDuplicated(inventory_basenames) == 0L,
  setequal(basename(cell_type_bed_paths), inventory_basenames)
)
public_inventory <- readr::read_tsv(
  output_value("output_inventory"),
  show_col_types = FALSE,
  progress = FALSE
)
stopifnot(
  identical(names(public_inventory), names(inventory)),
  isTRUE(all.equal(
    as.data.frame(public_inventory),
    as.data.frame(inventory),
    check.attributes = FALSE
  ))
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
    identical(entry$path, entry$file_name),
    grepl("^[^/]+[.]bed[.]gz$", entry$path),
    identical(entry$scale, "log2_cpm"),
    identical(unlist(entry$dimensions, use.names = FALSE), c(
      length(expected_gene_ids),
      length(expected_samples)
    )),
    entry$cell_group %in% expected_groups
  )
})

qc_summary <- readr::read_tsv(
  output_value("qc_summary"),
  show_col_types = FALSE,
  progress = FALSE
)
required_qc_metrics <- c(
  "gene_count", "sample_count", "cell_group_count",
  "excluded_constant_gene_count",
  "correlation_min", "correlation_median", "correlation_mean",
  "correlation_max", "rmse_min", "rmse_median", "rmse_mean",
  "rmse_max", "lm22_gene_count", "lm22_cell_type_count",
  "lm22_value_validation", "lm22_value_min", "lm22_value_max",
  "input_proportion_max_row_sum_error",
  "combined_proportion_max_row_sum_error",
  "adjusted_weight_max_row_sum_error",
  "normalization_adjustment_max_abs", "zero_values_adjusted",
  "tca_internal_iterations", "tca_max_internal_iterations",
  "tca_convergence", "tca_tau_hat"
)
stopifnot(setequal(required_qc_metrics, qc_summary$metric))
qc_status <- stats::setNames(qc_summary$status, qc_summary$metric)
if (identical(expected_proportion_mode, "dtangle")) {
  stopifnot(identical(qc_status[["lm22_value_validation"]], "passed"))
} else {
  stopifnot(identical(
    qc_status[["lm22_value_validation"]],
    "not_applicable_precomputed_mode"
  ))
}
stopifnot(qc_status[["tca_convergence"]] %in% c(
  "converged", "max_iterations_reached"
))

required_file_outputs <- c(
  "proportion_mode_validation_log", "proportions_lm22",
  "proportions_combined", "tca_weights", "cell_group_filter_report",
  "proportions_log", "tca_model", "tca_model_log", "tca_expression",
  "tca_excluded_genes", "fit_tca_log", "cell_type_bed_inventory",
  "reconstruction_by_sample", "qc_summary", "qc_plots", "export_log",
  "export_detail_log", "output_manifest", "output_inventory",
  "manifest_log", "effective_parameters_file"
)
if (identical(expected_proportion_mode, "dtangle")) {
  required_file_outputs <- c(required_file_outputs, c(
    "estimated_proportions", "dtangle_markers", "dtangle_metadata",
    "dtangle_overlap_report", "transformed_lm22", "dtangle_shared_bulk",
    "dtangle_log"
  ))
}
purrr::walk(required_file_outputs, function(name) {
  stopifnot(file.exists(output_value(name)))
})

message(sprintf(
  paste0(
    "Smoke assertions passed: %d genes, %d samples, %d groups, ",
    "scale=log2_cpm"
  ),
  length(expected_gene_ids),
  length(expected_samples),
  length(expected_groups)
))
