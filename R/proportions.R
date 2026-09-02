validate_lm22_proportions <- function(proportions) {
  if (!is.matrix(proportions) || !is.numeric(proportions)) {
    stop("LM22 proportions must be a numeric matrix", call. = FALSE)
  }
  if (nrow(proportions) == 0L) {
    stop("LM22 proportions must contain at least one sample", call. = FALSE)
  }

  sample_ids <- rownames(proportions)
  if (is.null(sample_ids) || anyNA(sample_ids) || any(!nzchar(sample_ids)) ||
      anyDuplicated(sample_ids) > 0L) {
    stop("LM22 proportions must have unique, non-empty sample identifiers", call. = FALSE)
  }

  observed_columns <- colnames(proportions)
  required_columns <- lm22_cell_types()
  if (is.null(observed_columns) || length(observed_columns) != 22L ||
      anyDuplicated(observed_columns) > 0L ||
      !setequal(observed_columns, required_columns)) {
    stop("LM22 proportions must contain exactly the 22 standard LM22 columns", call. = FALSE)
  }
  if (any(!is.finite(proportions))) {
    stop("LM22 proportions must be finite", call. = FALSE)
  }
  if (any(proportions < 0)) {
    stop("LM22 proportions must be nonnegative", call. = FALSE)
  }

  row_sum_difference <- abs(rowSums(proportions) - 1)
  if (any(row_sum_difference > 1e-6 + .Machine$double.eps)) {
    stop("LM22 proportion rows must sum to one within 1e-6", call. = FALSE)
  }

  invisible(TRUE)
}

validate_group_parameters <- function(mean_threshold, zero_floor) {
  if (!is.numeric(mean_threshold) || length(mean_threshold) != 1L ||
      !is.finite(mean_threshold) || mean_threshold < 0) {
    stop("mean_threshold must be one finite nonnegative value", call. = FALSE)
  }
  if (!is.numeric(zero_floor) || length(zero_floor) != 1L ||
      !is.finite(zero_floor) || zero_floor <= 0) {
    stop("zero_floor must be one finite positive value", call. = FALSE)
  }

  invisible(TRUE)
}

validate_combined_proportions <- function(combined) {
  if (!is.matrix(combined) || !is.numeric(combined) || nrow(combined) == 0L ||
      ncol(combined) == 0L) {
    stop("Combined proportions must be a non-empty numeric matrix", call. = FALSE)
  }
  if (is.null(rownames(combined)) || is.null(colnames(combined)) ||
      anyNA(rownames(combined)) || any(!nzchar(rownames(combined))) ||
      anyDuplicated(rownames(combined)) > 0L || anyNA(colnames(combined)) ||
      any(!nzchar(colnames(combined))) || anyDuplicated(colnames(combined)) > 0L) {
    stop("Combined proportions must have unique, non-empty sample and group identifiers", call. = FALSE)
  }
  if (any(!is.finite(combined))) {
    stop("Combined proportions must be finite", call. = FALSE)
  }
  if (any(combined < 0)) {
    stop("Combined proportions must be nonnegative", call. = FALSE)
  }

  invisible(TRUE)
}

combine_lm22_proportions <- function(proportions) {
  validate_lm22_proportions(proportions)
  group_map <- lm22_group_map()
  combined <- vapply(
    group_map,
    function(cell_types) rowSums(proportions[, cell_types, drop = FALSE]),
    numeric(nrow(proportions))
  )
  if (nrow(proportions) == 1L) {
    combined <- matrix(
      combined,
      nrow = 1L,
      dimnames = list(rownames(proportions), names(group_map))
    )
  } else {
    rownames(combined) <- rownames(proportions)
  }
  combined
}

filter_and_adjust_groups <- function(
    combined,
    mean_threshold = pipeline_defaults()$group_mean_threshold,
    zero_floor = pipeline_defaults()$zero_floor) {
  validate_combined_proportions(combined)
  validate_group_parameters(mean_threshold, zero_floor)

  cohort_means <- colMeans(combined)
  retained <- cohort_means >= mean_threshold
  report <- tibble::tibble(
    cell_group = colnames(combined),
    cohort_mean = unname(cohort_means),
    threshold = mean_threshold,
    retained = retained,
    filter_reason = ifelse(retained, "retained", "below_threshold"),
    zero_count_before = colSums(combined == 0),
    zero_floor = zero_floor
  )
  if (sum(retained) < 2L) {
    stop("At least two cell groups must be retained for TCA", call. = FALSE)
  }

  adjusted <- combined[, retained, drop = FALSE]
  adjusted[adjusted == 0] <- zero_floor
  weights <- sweep(adjusted, 1L, rowSums(adjusted), "/")

  list(weights = weights, report = report)
}

process_proportions <- function(
    proportions,
    mean_threshold = pipeline_defaults()$group_mean_threshold,
    zero_floor = pipeline_defaults()$zero_floor) {
  validate_lm22_proportions(proportions)
  validate_group_parameters(mean_threshold, zero_floor)

  combined <- combine_lm22_proportions(proportions)
  filtered <- filter_and_adjust_groups(combined, mean_threshold, zero_floor)

  list(
    original = proportions,
    combined = combined,
    tca_weights = filtered$weights,
    report = filtered$report
  )
}
