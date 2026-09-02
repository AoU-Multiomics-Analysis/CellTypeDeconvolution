validate_reconstruction_inputs <- function(tensor, weights, C2, deltas_hat) {
  validate_tensor_list(tensor)
  validate_tensor_matrix(weights, "Reconstruction weights")
  source_names <- names(tensor)
  gene_names <- rownames(tensor[[1L]])
  sample_ids <- colnames(tensor[[1L]])
  if (!identical(colnames(weights), source_names)) {
    stop("Reconstruction weight source order must match tensor source order",
      call. = FALSE
    )
  }
  if (!identical(rownames(weights), sample_ids)) {
    stop("Reconstruction weight sample order must match tensor sample order",
      call. = FALSE
    )
  }
  if (is.null(C2)) {
    if (!is.null(deltas_hat) && length(deltas_hat) > 0L) {
      stop("deltas_hat requires C2", call. = FALSE)
    }
    return(invisible(TRUE))
  }
  if (!is.matrix(C2) || !is.numeric(C2) || any(!is.finite(C2)) ||
      is.null(rownames(C2)) || is.null(colnames(C2))) {
    stop("C2 must be a finite matrix with identifiers", call. = FALSE)
  }
  if (!identical(rownames(C2), sample_ids)) {
    stop("C2 sample order must match tensor sample order", call. = FALSE)
  }
  if (!is.matrix(deltas_hat) || !is.numeric(deltas_hat) ||
      any(!is.finite(deltas_hat)) ||
      !identical(rownames(deltas_hat), gene_names) ||
      !identical(colnames(deltas_hat), colnames(C2))) {
    stop("deltas_hat dimensions and order must match genes and C2",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

reconstruct_tensor <- function(
    tensor,
    weights,
    C2 = NULL,
    deltas_hat = NULL) {
  validate_reconstruction_inputs(tensor, weights, C2, deltas_hat)
  source_reconstruction <- matrix(
    0,
    nrow = nrow(tensor[[1L]]),
    ncol = ncol(tensor[[1L]]),
    dimnames = dimnames(tensor[[1L]])
  )
  for (source_index in seq_along(tensor)) {
    for (sample_index in seq_len(ncol(source_reconstruction))) {
      source_reconstruction[, sample_index] <-
        source_reconstruction[, sample_index] +
        tensor[[source_index]][, sample_index] *
        weights[sample_index, source_index]
    }
  }
  if (!is.null(C2) && ncol(C2) > 0L) {
    source_reconstruction <- source_reconstruction + t(C2 %*% t(deltas_hat))
  }
  if (any(!is.finite(source_reconstruction))) {
    stop("Reconstructed values must be finite", call. = FALSE)
  }
  source_reconstruction
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

initialize_reconstruction_stats <- function(sample_ids) {
  if (!is.character(sample_ids) || length(sample_ids) == 0L ||
      anyNA(sample_ids) || any(!nzchar(sample_ids)) ||
      anyDuplicated(sample_ids) > 0L) {
    stop("sample_ids must be a non-empty vector of unique identifiers",
      call. = FALSE
    )
  }
  tibble::tibble(
    sample_id = sample_ids,
    n = rep(0, length(sample_ids)),
    sum_x = rep(0, length(sample_ids)),
    sum_y = rep(0, length(sample_ids)),
    sum_x2 = rep(0, length(sample_ids)),
    sum_y2 = rep(0, length(sample_ids)),
    sum_xy = rep(0, length(sample_ids)),
    sum_squared_error = rep(0, length(sample_ids))
  )
}

update_reconstruction_stats <- function(stats, observed, reconstructed) {
  validate_tensor_matrix(observed, "Observed expression")
  validate_tensor_matrix(reconstructed, "Reconstructed expression")
  if (!identical(dim(observed), dim(reconstructed))) {
    stop("Observed and reconstructed dimensions must match", call. = FALSE)
  }
  if (!identical(dimnames(observed), dimnames(reconstructed))) {
    stop("Observed and reconstructed order must match", call. = FALSE)
  }
  if (!identical(stats$sample_id, colnames(observed))) {
    stop("Reconstruction statistic sample order does not match", call. = FALSE)
  }
  stats$n <- stats$n + nrow(observed)
  stats$sum_x <- unname(stats$sum_x + colSums(observed))
  stats$sum_y <- unname(stats$sum_y + colSums(reconstructed))
  stats$sum_x2 <- unname(stats$sum_x2 + colSums(observed^2))
  stats$sum_y2 <- unname(stats$sum_y2 + colSums(reconstructed^2))
  stats$sum_xy <- unname(stats$sum_xy + colSums(observed * reconstructed))
  stats$sum_squared_error <- stats$sum_squared_error +
    unname(colSums((observed - reconstructed)^2))
  stats
}

finalize_reconstruction_stats <- function(stats) {
  numerator <- stats$n * stats$sum_xy - stats$sum_x * stats$sum_y
  denominator <- sqrt(
    (stats$n * stats$sum_x2 - stats$sum_x^2) *
      (stats$n * stats$sum_y2 - stats$sum_y^2)
  )
  if (any(stats$n < 2) || any(!is.finite(denominator)) || any(denominator <= 0)) {
    stop("Per-sample reconstruction correlation is not finite", call. = FALSE)
  }
  metrics <- tibble::tibble(
    sample_id = stats$sample_id,
    gene_count = as.integer(stats$n),
    correlation = numerator / denominator,
    rmse = sqrt(stats$sum_squared_error / stats$n)
  )
  if (any(!is.finite(metrics$correlation)) || any(!is.finite(metrics$rmse))) {
    stop("Per-sample reconstruction metrics must be finite", call. = FALSE)
  }
  metrics
}

make_qc_plots <- function(weights, observed, reconstructed) {
  validate_tensor_matrix(weights, "QC weights")
  if (!is.numeric(observed) || !is.numeric(reconstructed) ||
      length(observed) == 0L || length(observed) != length(reconstructed) ||
      any(!is.finite(observed)) || any(!is.finite(reconstructed))) {
    stop("QC observed and reconstructed values must be finite and aligned",
      call. = FALSE
    )
  }
  proportion_data <- tibble::as_tibble(
    weights,
    rownames = "sample_id",
    .name_repair = "minimal"
  ) |>
    tidyr::pivot_longer(
      cols = -"sample_id",
      names_to = "cell_group",
      values_to = "proportion"
    )
  reconstruction_data <- tibble::tibble(
    observed = observed,
    reconstructed = reconstructed,
    residual = observed - reconstructed
  )
  minimal_theme <- ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank()
    )

  proportion_plot <- ggplot2::ggplot(
    proportion_data,
    ggplot2::aes(x = .data$cell_group, y = .data$proportion)
  ) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.12, alpha = 0.25, size = 0.5) +
    ggplot2::labs(x = "Cell group", y = "Proportion") +
    minimal_theme +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    ))
  reconstruction_plot <- ggplot2::ggplot(
    reconstruction_data,
    ggplot2::aes(x = .data$observed, y = .data$reconstructed)
  ) +
    ggplot2::geom_point(alpha = 0.25, size = 0.6) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::labs(
      x = "Observed log2(CPM)",
      y = "Reconstructed log2(CPM)"
    ) +
    minimal_theme
  residual_plot <- ggplot2::ggplot(
    reconstruction_data,
    ggplot2::aes(x = .data$residual)
  ) +
    ggplot2::geom_histogram(bins = 50L) +
    ggplot2::labs(x = "Residual log2(CPM)", y = "Count") +
    minimal_theme

  list(
    proportions = proportion_plot,
    reconstruction = reconstruction_plot,
    residuals = residual_plot
  )
}

make_qc_summary <- function(
    reconstruction_by_sample,
    gene_count,
    group_count,
    excluded_constant_gene_count = 0L) {
  tibble::tibble(
    metric = c(
      "gene_count", "sample_count", "cell_group_count",
      "excluded_constant_gene_count",
      "correlation_min", "correlation_median", "correlation_mean",
      "correlation_max", "rmse_min", "rmse_median", "rmse_mean",
      "rmse_max"
    ),
    value = c(
      gene_count,
      nrow(reconstruction_by_sample),
      group_count,
      excluded_constant_gene_count,
      min(reconstruction_by_sample$correlation),
      stats::median(reconstruction_by_sample$correlation),
      mean(reconstruction_by_sample$correlation),
      max(reconstruction_by_sample$correlation),
      min(reconstruction_by_sample$rmse),
      stats::median(reconstruction_by_sample$rmse),
      mean(reconstruction_by_sample$rmse),
      max(reconstruction_by_sample$rmse)
    )
  )
}

write_qc_reports <- function(
    weights,
    reconstruction_by_sample,
    observed,
    reconstructed,
    gene_count,
    excluded_constant_gene_count = 0L,
    output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_paths <- c(
    reconstruction_by_sample = file.path(
      output_dir, "reconstruction_by_sample.tsv"
    ),
    qc_summary = file.path(output_dir, "qc_summary.tsv"),
    qc_plots = file.path(output_dir, "qc_plots.pdf")
  )
  summary <- make_qc_summary(
    reconstruction_by_sample,
    gene_count,
    ncol(weights),
    excluded_constant_gene_count
  )
  plots <- make_qc_plots(weights, observed, reconstructed)
  readr::write_tsv(
    reconstruction_by_sample,
    output_paths[["reconstruction_by_sample"]],
    na = ""
  )
  readr::write_tsv(summary, output_paths[["qc_summary"]], na = "")
  grDevices::pdf(output_paths[["qc_plots"]], width = 8, height = 6)
  on.exit(grDevices::dev.off(), add = TRUE)
  purrr::walk(plots, print)
  unname(grDevices::dev.off())
  on.exit(NULL, add = FALSE)
  output_paths
}

parse_tca_convergence <- function(tca_log_lines) {
  if (!is.character(tca_log_lines) || length(tca_log_lines) == 0L) {
    stop("TCA model log must contain text lines", call. = FALSE)
  }
  iteration_matches <- stringr::str_match(
    tca_log_lines,
    "Iteration ([0-9]+) out of ([0-9]+) internal iterations"
  )
  observed <- !is.na(iteration_matches[, 1L])
  if (!any(observed)) {
    stop("TCA model log does not contain internal iteration records", call. = FALSE)
  }
  iterations <- as.integer(iteration_matches[observed, 2L])
  maximum_iterations <- as.integer(iteration_matches[observed, 3L])
  if (anyNA(iterations) || anyNA(maximum_iterations)) {
    stop("TCA model log contains invalid iteration records", call. = FALSE)
  }
  converged <- any(grepl(
    "Internal loop converged.",
    tca_log_lines,
    fixed = TRUE
  ))
  list(
    iterations = max(iterations),
    maximum_iterations = max(maximum_iterations),
    converged = converged,
    status = if (converged) "converged" else "max_iterations_reached"
  )
}

build_pipeline_qc_summary <- function(
    export_summary,
    original_proportions,
    combined_proportions,
    tca_weights,
    filter_report,
    tca_model,
    tca_log_lines,
    dtangle_metadata = NULL) {
  if (!inherits(export_summary, "data.frame") ||
      !all(c("metric", "value") %in% names(export_summary)) ||
      nrow(export_summary) == 0L || anyNA(export_summary$metric) ||
      any(!nzchar(export_summary$metric)) ||
      anyDuplicated(export_summary$metric) > 0L) {
    stop("export_summary must contain unique metric and value columns", call. = FALSE)
  }
  matrices <- list(
    original_proportions = original_proportions,
    combined_proportions = combined_proportions,
    tca_weights = tca_weights
  )
  valid_matrices <- purrr::map_lgl(matrices, function(value) {
    is.matrix(value) && is.numeric(value) && nrow(value) > 0L &&
      ncol(value) > 0L && all(is.finite(value)) &&
      !is.null(rownames(value)) && !is.null(colnames(value))
  })
  if (!all(valid_matrices)) {
    stop("Proportion QC inputs must be finite non-empty matrices", call. = FALSE)
  }
  if (!identical(rownames(original_proportions), rownames(combined_proportions)) ||
      !identical(rownames(combined_proportions), rownames(tca_weights))) {
    stop("Proportion QC sample order must match exactly", call. = FALSE)
  }
  required_filter_columns <- c(
    "cell_group", "retained", "zero_count_before", "zero_floor"
  )
  if (!inherits(filter_report, "data.frame") ||
      !all(required_filter_columns %in% names(filter_report))) {
    stop("filter_report is missing required adjustment columns", call. = FALSE)
  }
  retained_report <- filter_report |>
    dplyr::filter(.data$retained) |>
    dplyr::select(dplyr::all_of(required_filter_columns))
  if (!identical(retained_report$cell_group, colnames(tca_weights)) ||
      !all(retained_report$cell_group %in% colnames(combined_proportions))) {
    stop("Retained group order must match TCA weight order", call. = FALSE)
  }
  adjusted_before_normalization <- combined_proportions[
    , retained_report$cell_group, drop = FALSE
  ]
  for (group_index in seq_len(ncol(adjusted_before_normalization))) {
    zero_values <- adjusted_before_normalization[, group_index] == 0
    adjusted_before_normalization[zero_values, group_index] <-
      retained_report$zero_floor[[group_index]]
  }
  if (!is.list(tca_model) || !is.numeric(tca_model$tau_hat) ||
      length(tca_model$tau_hat) != 1L || !is.finite(tca_model$tau_hat)) {
    stop("TCA model must contain one finite tau_hat value", call. = FALSE)
  }
  convergence <- parse_tca_convergence(tca_log_lines)
  base_summary <- tibble::as_tibble(export_summary) |>
    dplyr::transmute(
      metric = as.character(.data$metric),
      value = as.numeric(.data$value),
      status = "observed"
    )
  proportion_summary <- tibble::tibble(
    metric = c(
      "input_proportion_max_row_sum_error",
      "combined_proportion_max_row_sum_error",
      "adjusted_weight_max_row_sum_error",
      "normalization_adjustment_max_abs",
      "zero_values_adjusted",
      "tca_internal_iterations",
      "tca_max_internal_iterations",
      "tca_convergence",
      "tca_tau_hat"
    ),
    value = c(
      max(abs(rowSums(original_proportions) - 1)),
      max(abs(rowSums(combined_proportions) - 1)),
      max(abs(rowSums(tca_weights) - 1)),
      max(abs(tca_weights - adjusted_before_normalization)),
      sum(retained_report$zero_count_before),
      convergence$iterations,
      convergence$maximum_iterations,
      as.numeric(convergence$converged),
      tca_model$tau_hat
    ),
    status = c(
      rep("passed", 7L),
      convergence$status,
      "fitted"
    )
  )
  lm22_status <- if (is.null(dtangle_metadata)) {
    tibble::tibble(
      metric = c(
        "lm22_gene_count", "lm22_cell_type_count",
        "lm22_value_min", "lm22_value_max", "lm22_value_validation"
      ),
      value = rep(NA_real_, 5L),
      status = "not_applicable_precomputed_mode"
    )
  } else {
    required_lm22_fields <- c(
      "gene_count", "cell_type_count", "value_min", "value_max",
      "validation_status"
    )
    if (!is.list(dtangle_metadata) ||
        !is.list(dtangle_metadata$lm22_qc) ||
        !all(required_lm22_fields %in% names(dtangle_metadata$lm22_qc)) ||
        !identical(dtangle_metadata$lm22_qc$validation_status, "passed")) {
      stop("dtangle metadata is missing passed LM22 QC", call. = FALSE)
    }
    qc <- dtangle_metadata$lm22_qc
    tibble::tibble(
      metric = c(
        "lm22_gene_count", "lm22_cell_type_count",
        "lm22_value_min", "lm22_value_max", "lm22_value_validation"
      ),
      value = c(
        as.numeric(qc$gene_count),
        as.numeric(qc$cell_type_count),
        as.numeric(qc$value_min),
        as.numeric(qc$value_max),
        1
      ),
      status = "passed"
    )
  }
  dplyr::bind_rows(base_summary, lm22_status, proportion_summary)
}

validate_manifest_outputs <- function(outputs) {
  required_columns <- c(
    "logical_name", "path", "n_genes", "n_samples", "scale",
    "cell_group"
  )
  if (!inherits(outputs, "data.frame") ||
      !all(required_columns %in% names(outputs))) {
    stop(
      paste0(
        "outputs must be a data frame with columns: ",
        paste(required_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  outputs <- tibble::as_tibble(outputs)
  character_fields <- c("logical_name", "path", "scale", "cell_group")
  valid_character_fields <- purrr::map_lgl(
    outputs[character_fields],
    ~ is.character(.x) && !anyNA(.x) && all(nzchar(trimws(.x)))
  )
  if (!all(valid_character_fields)) {
    stop(
      "Manifest logical names, paths, scales, and cell groups must be non-empty",
      call. = FALSE
    )
  }
  if (nrow(outputs) == 0L || anyDuplicated(outputs$logical_name) > 0L ||
      anyDuplicated(outputs$path) > 0L ||
      anyDuplicated(outputs$cell_group) > 0L ||
      any(!file.exists(outputs$path))) {
    stop("Manifest outputs must be unique existing files", call. = FALSE)
  }
  dimensions <- c(outputs$n_genes, outputs$n_samples)
  if (!is.numeric(dimensions) || anyNA(dimensions) ||
      any(!is.finite(dimensions)) || any(dimensions < 1) ||
      any(dimensions != floor(dimensions))) {
    stop(
      "Manifest output dimensions must be positive integer values",
      call. = FALSE
    )
  }
  outputs$n_genes <- as.integer(outputs$n_genes)
  outputs$n_samples <- as.integer(outputs$n_samples)
  outputs
}

validate_container_image <- function(container_image) {
  local_smoke_image <- "celltype-deconvolution:test"
  default_image <- "ghcr.io/aou-multiomics-analysis/celltypedeconvolution:latest"
  digest_pattern <- "^[^[:space:]@]+@sha256:[0-9a-f]{64}$"
  valid <- is.character(container_image) && length(container_image) == 1L &&
    !is.na(container_image) &&
    (identical(container_image, local_smoke_image) ||
      identical(container_image, default_image) ||
      grepl(digest_pattern, container_image))
  if (!valid) {
    stop(
      paste0(
        "container_image must use an immutable SHA-256 digest or equal an ",
        "approved local smoke or GitHub default image"
      ),
      call. = FALSE
    )
  }
  container_image
}

build_output_manifest <- function(
    outputs,
    tca_version,
    parameters,
    container_image) {
  outputs <- validate_manifest_outputs(outputs)
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("The digest package is required for SHA-256 checksums", call. = FALSE)
  }
  if (!is.list(parameters)) {
    stop("parameters must be a named list", call. = FALSE)
  }
  container_image <- validate_container_image(container_image)
  output_entries <- purrr::pmap(outputs, function(
      logical_name,
      path,
      n_genes,
      n_samples,
      scale,
      cell_group,
      ...) {
    list(
      logical_name = as.character(logical_name),
      file_name = basename(path),
      path = basename(path),
      sha256 = digest::digest(
        file = path,
        algo = "sha256",
        serialize = FALSE
      ),
      dimensions = c(as.integer(n_genes), as.integer(n_samples)),
      scale = as.character(scale),
      cell_group = as.character(cell_group)
    )
  })
  list(
    schema_version = "1.0",
    created_utc = format(
      Sys.time(),
      "%Y-%m-%dT%H:%M:%SZ",
      tz = "UTC"
    ),
    tca_version = tca_version,
    software_versions = list(
      R = as.character(getRversion()),
      TCA = tca_version
    ),
    parameters = parameters,
    container_image = container_image,
    outputs = output_entries
  )
}

write_output_manifest <- function(path, manifest) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("The jsonlite package is required for manifest output", call. = FALSE)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(
    manifest,
    path,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )
  invisible(path)
}
