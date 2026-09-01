read_numeric_matrix <- function(path, id_column) {
  if (!is.character(id_column) || length(id_column) != 1L || is.na(id_column)) {
    stop("id_column must be one non-missing character value", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop(sprintf("Matrix file does not exist: %s", path), call. = FALSE)
  }

  table <- readr::read_tsv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    name_repair = "minimal",
    progress = FALSE,
    show_col_types = FALSE
  )
  column_names <- names(table)
  if (length(column_names) == 0L || !identical(column_names[[1L]], id_column)) {
    stop(sprintf("The first column must be '%s'", id_column), call. = FALSE)
  }
  if (anyDuplicated(column_names) > 0L) {
    stop("Matrix columns must have unique names", call. = FALSE)
  }

  identifiers <- table[[id_column]]
  if (anyNA(identifiers) || any(!nzchar(identifiers))) {
    stop(sprintf("ID column '%s' contains missing or empty identifiers", id_column), call. = FALSE)
  }
  if (anyDuplicated(identifiers) > 0L) {
    duplicate_id <- identifiers[[which(duplicated(identifiers))[1L]]]
    stop(sprintf("ID column '%s' contains duplicate identifier '%s'", id_column, duplicate_id), call. = FALSE)
  }

  value_names <- setdiff(column_names, id_column)
  if (length(value_names) == 0L) {
    stop("Matrix must contain at least one numeric column", call. = FALSE)
  }

  values <- lapply(value_names, function(value_name) {
    raw_values <- table[[value_name]]
    parsed_values <- suppressWarnings(readr::parse_double(raw_values))
    invalid_numeric <- is.na(parsed_values)
    if (any(invalid_numeric)) {
      stop(sprintf("Column '%s' contains missing or nonnumeric values", value_name), call. = FALSE)
    }
    if (any(!is.finite(parsed_values))) {
      stop(sprintf("Column '%s' contains nonfinite values", value_name), call. = FALSE)
    }
    parsed_values
  })

  matrix_values <- do.call(cbind, values)
  dimnames(matrix_values) <- list(identifiers, value_names)
  matrix_values
}

write_numeric_matrix <- function(x, path, id_column) {
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("x must be a numeric matrix", call. = FALSE)
  }
  if (!is.character(id_column) || length(id_column) != 1L || is.na(id_column)) {
    stop("id_column must be one non-missing character value", call. = FALSE)
  }
  if (is.null(rownames(x)) || anyNA(rownames(x)) || any(!nzchar(rownames(x)))) {
    stop("x must have non-missing, non-empty row names", call. = FALSE)
  }
  if (anyDuplicated(rownames(x)) > 0L) {
    stop("x must have unique row names", call. = FALSE)
  }
  if (is.null(colnames(x)) || anyNA(colnames(x)) || any(!nzchar(colnames(x)))) {
    stop("x must have non-missing, non-empty column names", call. = FALSE)
  }
  if (anyDuplicated(colnames(x)) > 0L || id_column %in% colnames(x)) {
    stop("x must have unique column names that do not include id_column", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop("x contains nonfinite values", call. = FALSE)
  }

  output <- tibble::rownames_to_column(as.data.frame(x, check.names = FALSE), var = id_column)
  readr::write_tsv(output, path, na = "")
  invisible(path)
}

assert_identical_ids <- function(expected, observed, label) {
  if (identical(expected, observed)) {
    return(invisible(TRUE))
  }

  total_ids <- max(length(expected), length(observed))
  first_difference <- which(vapply(seq_len(total_ids), function(index) {
    if (index > length(expected) || index > length(observed)) {
      return(TRUE)
    }
    !identical(expected[[index]], observed[[index]])
  }, logical(1)))[[1L]]
  expected_value <- if (first_difference <= length(expected)) expected[[first_difference]] else "<missing>"
  observed_value <- if (first_difference <= length(observed)) observed[[first_difference]] else "<missing>"
  stop(
    sprintf(
      "%s identifiers differ at position %d: expected '%s'; observed '%s'",
      label, first_difference, expected_value, observed_value
    ),
    call. = FALSE
  )
}
