validate_lm22 <- function(lm22_linear) {
  if (!is.matrix(lm22_linear) || !is.numeric(lm22_linear)) {
    stop("LM22 must be a numeric matrix", call. = FALSE)
  }
  if (nrow(lm22_linear) == 0L) {
    stop("LM22 must contain at least one gene", call. = FALSE)
  }

  required_cell_types <- lm22_cell_types()
  observed_cell_types <- colnames(lm22_linear)
  if (is.null(observed_cell_types) || length(observed_cell_types) != 22L ||
      anyDuplicated(observed_cell_types) > 0L ||
      !setequal(observed_cell_types, required_cell_types)) {
    stop("LM22 must contain exactly the 22 standard LM22 columns", call. = FALSE)
  }

  gene_symbols <- rownames(lm22_linear)
  if (is.null(gene_symbols)) {
    stop("LM22 gene symbols must be present", call. = FALSE)
  }
  gene_symbols <- trimws(gene_symbols)
  if (anyNA(gene_symbols) || any(!nzchar(gene_symbols))) {
    stop("LM22 gene symbols must be non-missing and non-empty", call. = FALSE)
  }
  if (anyDuplicated(gene_symbols) > 0L) {
    stop("LM22 gene symbols must be unique after trimming", call. = FALSE)
  }
  if (any(!is.finite(lm22_linear))) {
    stop("LM22 values must be finite", call. = FALSE)
  }
  if (any(lm22_linear <= 0)) {
    stop("LM22 values must be positive", call. = FALSE)
  }

  invisible(TRUE)
}

standardize_lm22 <- function(lm22_linear) {
  validate_lm22(lm22_linear)
  lm22_linear <- lm22_linear[, lm22_cell_types(), drop = FALSE]
  rownames(lm22_linear) <- trimws(rownames(lm22_linear))
  lm22_linear
}

transform_lm22 <- function(lm22_linear) {
  log2(standardize_lm22(lm22_linear))
}

validate_bulk_log <- function(bulk_log) {
  if (!is.matrix(bulk_log) || !is.numeric(bulk_log)) {
    stop("Bulk log expression must be a numeric matrix", call. = FALSE)
  }
  if (nrow(bulk_log) == 0L || ncol(bulk_log) == 0L) {
    stop("Bulk log expression must contain at least one gene and one sample", call. = FALSE)
  }
  if (is.null(rownames(bulk_log)) || is.null(colnames(bulk_log))) {
    stop("Bulk log expression must have gene and sample identifiers", call. = FALSE)
  }

  gene_symbols <- trimws(rownames(bulk_log))
  if (anyNA(gene_symbols) || any(!nzchar(gene_symbols)) ||
      anyDuplicated(gene_symbols) > 0L) {
    stop("Bulk log expression gene symbols must be unique after trimming", call. = FALSE)
  }
  if (anyNA(colnames(bulk_log)) || any(!nzchar(colnames(bulk_log))) ||
      anyDuplicated(colnames(bulk_log)) > 0L) {
    stop("Bulk log expression sample identifiers must be unique and non-empty", call. = FALSE)
  }
  if (any(!is.finite(bulk_log))) {
    stop("Bulk log expression values must be finite", call. = FALSE)
  }
  if (any(bulk_log < 0)) {
    stop("Bulk log expression values must be nonnegative", call. = FALSE)
  }

  bulk_log
}

standardize_bulk_log <- function(bulk_log) {
  bulk_log <- validate_bulk_log(bulk_log)
  rownames(bulk_log) <- trimws(rownames(bulk_log))
  bulk_log
}

prepare_dtangle_inputs <- function(
    bulk_log,
    lm22_linear,
    min_overlap = pipeline_defaults()$min_lm22_overlap,
    quantile_normalize = FALSE) {
  if (!is.numeric(min_overlap) || length(min_overlap) != 1L ||
      !is.finite(min_overlap) || min_overlap <= 0 || min_overlap > 1) {
    stop("min_overlap must be a finite value in (0, 1]", call. = FALSE)
  }
  if (!is.logical(quantile_normalize) || length(quantile_normalize) != 1L ||
      is.na(quantile_normalize)) {
    stop("quantile_normalize must be one non-missing logical value", call. = FALSE)
  }

  transformed_lm22 <- transform_lm22(lm22_linear)
  bulk_log <- standardize_bulk_log(bulk_log)
  shared_genes <- rownames(transformed_lm22)[
    rownames(transformed_lm22) %in% rownames(bulk_log)
  ]
  overlap_fraction <- length(shared_genes) / nrow(transformed_lm22)
  if (overlap_fraction < min_overlap) {
    stop(
      sprintf("LM22 overlap %.3f is below %.3f.", overlap_fraction, min_overlap),
      call. = FALSE
    )
  }

  shared_lm22 <- transformed_lm22[shared_genes, , drop = FALSE]
  shared_bulk <- bulk_log[match(shared_genes, rownames(bulk_log)), , drop = FALSE]
  if (isTRUE(quantile_normalize)) {
    joined_profiles <- cbind(shared_lm22, shared_bulk)
    normalized_profiles <- limma::normalizeBetweenArrays(joined_profiles)
    shared_lm22 <- normalized_profiles[, colnames(shared_lm22), drop = FALSE]
    shared_bulk <- normalized_profiles[, colnames(shared_bulk), drop = FALSE]
  }

  list(
    Y = t(shared_bulk),
    references = t(shared_lm22),
    transformed_lm22 = transformed_lm22,
    shared_bulk = shared_bulk,
    overlap_report = tibble::tibble(
      reference_gene_count = nrow(transformed_lm22),
      overlap_count = length(shared_genes),
      overlap_fraction = overlap_fraction,
      min_overlap = min_overlap,
      quantile_normalize = quantile_normalize
    )
  )
}

estimate_dtangle <- function(
    inputs,
    marker_fraction = pipeline_defaults()$marker_fraction) {
  if (!is.list(inputs) || !all(c("Y", "references", "overlap_report") %in% names(inputs))) {
    stop("inputs must be returned by prepare_dtangle_inputs", call. = FALSE)
  }
  if (!is.numeric(marker_fraction) || length(marker_fraction) != 1L ||
      !is.finite(marker_fraction) || marker_fraction <= 0 || marker_fraction > 1) {
    stop("marker_fraction must be a finite value in (0, 1]", call. = FALSE)
  }
  if (!requireNamespace("dtangle", quietly = TRUE)) {
    stop("The dtangle package is required for proportion estimation", call. = FALSE)
  }

  fit <- dtangle::dtangle(
    Y = inputs$Y,
    references = inputs$references,
    n_markers = marker_fraction,
    data_type = "rna-seq",
    marker_method = "ratio"
  )
  selected_marker_counts <- lengths(fit$markers)
  if (length(selected_marker_counts) != nrow(inputs$references) ||
      any(selected_marker_counts == 0L)) {
    stop("Every LM22 cell type must have at least one selected marker", call. = FALSE)
  }

  proportions <- fit$estimates
  if (is.null(dim(proportions))) {
    proportions <- matrix(
      proportions,
      nrow = nrow(inputs$Y),
      dimnames = list(rownames(inputs$Y), rownames(inputs$references))
    )
  }
  colnames(proportions) <- rownames(inputs$references)
  if (any(!is.finite(proportions)) || any(proportions < 0)) {
    stop("dtangle estimates must be finite and nonnegative", call. = FALSE)
  }
  if (any(abs(rowSums(proportions) - 1) > 1e-8)) {
    stop("dtangle estimate rows must sum to one within 1e-8", call. = FALSE)
  }

  markers <- purrr::map2_dfr(
    fit$markers,
    rownames(inputs$references),
    function(indices, cell_type) {
      tibble::tibble(
        cell_type = cell_type,
        marker_rank = seq_along(indices),
        gene_symbol = colnames(inputs$Y)[indices]
      )
    }
  )
  overlap <- inputs$overlap_report[1, , drop = FALSE]
  metadata <- list(
    dtangle_version = as.character(utils::packageVersion("dtangle")),
    gamma = unname(fit$gamma),
    marker_method = "ratio",
    marker_fraction = marker_fraction,
    overlap_count = overlap$overlap_count,
    overlap_fraction = overlap$overlap_fraction,
    quantile_normalize = overlap$quantile_normalize,
    sample_count = nrow(inputs$Y),
    reference_dimensions = list(
      cell_types = nrow(inputs$references),
      genes = ncol(inputs$references)
    )
  )

  list(proportions = proportions, markers = markers, metadata = metadata)
}
