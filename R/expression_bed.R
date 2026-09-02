expression_bed_columns <- function() c("#chr", "start", "end", "gene_id")

read_expression_bed <- function(path) {
  table <- readr::read_tsv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    name_repair = "minimal",
    show_col_types = FALSE,
    progress = FALSE
  )
  required <- expression_bed_columns()
  if (ncol(table) <= length(required) ||
      !identical(names(table)[seq_along(required)], required)) {
    stop(
      "Expression BED must start with #chr, start, end, and gene_id",
      call. = FALSE
    )
  }
  coordinates <- table |>
    dplyr::transmute(
      `#chr` = trimws(.data[["#chr"]]),
      start = readr::parse_integer(.data$start),
      end = readr::parse_integer(.data$end),
      gene_id = trimws(.data$gene_id)
    )
  if (anyNA(coordinates$start) || anyNA(coordinates$end) ||
      any(coordinates$start < 0L) ||
      any(coordinates$start >= coordinates$end)) {
    stop("BED start must be a non-negative integer less than end", call. = FALSE)
  }
  if (anyNA(coordinates[["#chr"]]) || anyNA(coordinates$gene_id) ||
      any(!nzchar(coordinates[["#chr"]])) ||
      any(!nzchar(coordinates$gene_id)) ||
      anyDuplicated(coordinates$gene_id) > 0L) {
    stop(
      paste0(
        "BED chromosome and gene_id values must not be missing or empty. ",
        "gene_id values must be unique."
      ),
      call. = FALSE
    )
  }
  sample_ids <- names(table)[-(seq_along(required))]
  if (any(!nzchar(sample_ids)) || anyDuplicated(sample_ids) > 0L) {
    stop("BED sample identifiers must be non-empty and unique", call. = FALSE)
  }
  cpm <- table |>
    dplyr::select(dplyr::all_of(sample_ids)) |>
    dplyr::mutate(dplyr::across(
      dplyr::everything(),
      ~ suppressWarnings(readr::parse_double(.x, na = character()))
    )) |>
    as.matrix()
  rownames(cpm) <- coordinates$gene_id
  cpm <- validate_cpm_matrix(cpm)
  if (any(cpm <= 0)) {
    stop("CPM values must be strictly positive", call. = FALSE)
  }
  list(coordinates = coordinates, cpm = cpm)
}

write_expression_bed <- function(path, coordinates, matrix) {
  required <- expression_bed_columns()
  if (!inherits(coordinates, "data.frame") ||
      !identical(names(coordinates), required)) {
    stop("coordinates must have the exact BED columns", call. = FALSE)
  }
  if (!is.matrix(matrix) || !is.numeric(matrix) ||
      !identical(rownames(matrix), coordinates$gene_id) ||
      any(!is.finite(matrix))) {
    stop("BED matrix genes and finite values must match coordinates", call. = FALSE)
  }
  output <- dplyr::bind_cols(
    tibble::as_tibble(coordinates),
    tibble::as_tibble(matrix, .name_repair = "minimal")
  )
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(output, path, na = "")
  invisible(path)
}
