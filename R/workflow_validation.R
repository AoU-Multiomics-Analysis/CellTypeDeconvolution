validate_proportion_mode <- function(
    lm22_defined,
    precomputed_proportions_defined) {
  values <- list(
    lm22_defined = lm22_defined,
    precomputed_proportions_defined = precomputed_proportions_defined
  )
  valid <- purrr::map_lgl(values, function(value) {
    is.logical(value) && length(value) == 1L && !is.na(value)
  })
  if (!all(valid)) {
    stop("Proportion-mode presence flags must be true or false", call. = FALSE)
  }
  if (identical(lm22_defined, precomputed_proportions_defined)) {
    stop(
      "exactly one of lm22 and precomputed_proportions must be provided",
      call. = FALSE
    )
  }
  if (lm22_defined) {
    list(selected_mode = "dtangle", estimate_proportions = TRUE)
  } else {
    list(selected_mode = "precomputed", estimate_proportions = FALSE)
  }
}
