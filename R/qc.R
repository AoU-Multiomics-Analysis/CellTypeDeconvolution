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
  source_reconstruction <- Reduce(
    `+`,
    purrr::map2(
      tensor,
      seq_len(ncol(weights)),
      ~ sweep(.x, 2L, weights[, .y], "*")
    )
  )
  covariate_term <- if (is.null(C2) || ncol(C2) == 0L) {
    matrix(
      0,
      nrow = nrow(source_reconstruction),
      ncol = ncol(source_reconstruction),
      dimnames = dimnames(source_reconstruction)
    )
  } else {
    t(C2 %*% t(deltas_hat))
  }
  reconstructed <- source_reconstruction + covariate_term
  if (any(!is.finite(reconstructed))) {
    stop("Reconstructed values must be finite", call. = FALSE)
  }
  reconstructed
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

make_qc_summary <- function(reconstruction_by_sample, gene_count, group_count) {
  tibble::tibble(
    metric = c(
      "gene_count", "sample_count", "cell_group_count",
      "correlation_min", "correlation_median", "correlation_mean",
      "correlation_max", "rmse_min", "rmse_median", "rmse_mean",
      "rmse_max"
    ),
    value = c(
      gene_count,
      nrow(reconstruction_by_sample),
      group_count,
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
    ncol(weights)
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
  if (nrow(outputs) == 0L || anyDuplicated(outputs$logical_name) > 0L ||
      any(!file.exists(outputs$path))) {
    stop("Manifest outputs must be unique existing files", call. = FALSE)
  }
  if (anyNA(outputs$n_genes) || anyNA(outputs$n_samples) ||
      any(outputs$n_genes < 1L) || any(outputs$n_samples < 1L)) {
    stop("Manifest output dimensions must be positive", call. = FALSE)
  }
  outputs
}

build_output_manifest <- function(
    outputs,
    pipeline_version,
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
      path = normalizePath(path),
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
    pipeline_version = pipeline_version,
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
