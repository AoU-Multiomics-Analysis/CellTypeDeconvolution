validate_expression_matrix <- function(expression) {
  if (!is.matrix(expression) || !is.numeric(expression)) {
    stop("expression must be a numeric matrix", call. = FALSE)
  }
  if (nrow(expression) == 0L || ncol(expression) == 0L) {
    stop("expression must contain at least one gene and one sample", call. = FALSE)
  }
  if (is.null(rownames(expression)) || anyNA(rownames(expression)) ||
      any(!nzchar(rownames(expression)))) {
    stop("gene identifiers must be non-missing and non-empty", call. = FALSE)
  }
  if (anyDuplicated(rownames(expression)) > 0L) {
    stop("gene identifiers must be unique", call. = FALSE)
  }
  if (is.null(colnames(expression)) || anyNA(colnames(expression)) ||
      any(!nzchar(colnames(expression)))) {
    stop("sample identifiers must be non-missing and non-empty", call. = FALSE)
  }
  if (anyDuplicated(colnames(expression)) > 0L) {
    stop("sample identifiers must be unique", call. = FALSE)
  }
  if (any(!is.finite(expression))) {
    stop("expression values must be finite", call. = FALSE)
  }
  if (any(expression < 0)) {
    stop("expression values must be nonnegative", call. = FALSE)
  }

  invisible(TRUE)
}

validate_annotation <- function(annotation) {
  if (!inherits(annotation, "data.frame")) {
    stop("annotation must be a data frame", call. = FALSE)
  }
  required_columns <- c("gene_id", "gene_symbol")
  missing_columns <- setdiff(required_columns, names(annotation))
  if (length(missing_columns) > 0L) {
    stop(
      sprintf(
        "annotation is missing required columns: %s",
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  annotation <- tibble::as_tibble(annotation)
  annotation$gene_id <- as.character(annotation$gene_id)
  annotation$gene_symbol <- trimws(as.character(annotation$gene_symbol))
  if (anyNA(annotation$gene_id) || any(!nzchar(annotation$gene_id))) {
    stop("annotation gene_id values must be non-missing and non-empty", call. = FALSE)
  }
  if (anyDuplicated(annotation$gene_id) > 0L) {
    stop("annotation gene_id values must be unique", call. = FALSE)
  }

  annotation
}

read_optional_annotation <- function(path) {
  if (is.null(path) || identical(path, "")) {
    return(NULL)
  }
  if (!is.character(path) || length(path) != 1L || is.na(path)) {
    stop("annotation path must be one non-missing character value", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop(sprintf("Annotation file does not exist: %s", path), call. = FALSE)
  }

  annotation <- readr::read_tsv(
    path,
    name_repair = "minimal",
    progress = FALSE,
    show_col_types = FALSE
  )
  validate_annotation(annotation)
}

counts_to_tpm <- function(counts, gene_length_bp) {
  validate_expression_matrix(counts)
  if (!identical(rownames(counts), names(gene_length_bp))) {
    stop("gene_length_bp names must exactly match count gene identifiers", call. = FALSE)
  }
  if (any(!is.finite(gene_length_bp)) || any(gene_length_bp <= 0)) {
    stop("Count input requires positive gene_length_bp values", call. = FALSE)
  }

  rates <- counts / (gene_length_bp / 1000)
  totals <- colSums(rates)
  if (any(totals <= 0)) {
    stop("Each sample must have a positive total expression rate.", call. = FALSE)
  }
  tpm <- sweep(rates, 2, totals, "/") * 1e6
  dimnames(tpm) <- dimnames(counts)
  tpm
}

collapse_to_symbols <- function(values, annotation) {
  validate_expression_matrix(values)
  annotation <- validate_annotation(annotation)
  selected_annotation <- annotation |>
    dplyr::filter(.data$gene_id %in% rownames(values)) |>
    dplyr::filter(!is.na(.data$gene_symbol), nzchar(.data$gene_symbol)) |>
    dplyr::select("gene_id", "gene_symbol")

  if (nrow(selected_annotation) == 0L) {
    return(matrix(
      numeric(),
      nrow = 0L,
      ncol = ncol(values),
      dimnames = list(character(), colnames(values))
    ))
  }

  collapsed <- tibble::as_tibble(values, rownames = "gene_id") |>
    tidyr::pivot_longer(
      cols = -"gene_id",
      names_to = "sample_id",
      values_to = "value"
    ) |>
    dplyr::inner_join(selected_annotation, by = "gene_id") |>
    dplyr::group_by(.data$gene_symbol, .data$sample_id) |>
    dplyr::summarise(value = sum(.data$value), .groups = "drop") |>
    dplyr::arrange(.data$gene_symbol) |>
    tidyr::pivot_wider(names_from = "sample_id", values_from = "value") |>
    dplyr::arrange(.data$gene_symbol) |>
    dplyr::select("gene_symbol", dplyr::all_of(colnames(values)))

  output <- as.matrix(dplyr::select(collapsed, -"gene_symbol"))
  storage.mode(output) <- "double"
  rownames(output) <- collapsed$gene_symbol
  output
}

normalise_tpm <- function(values) {
  totals <- colSums(values)
  if (any(totals <= 0)) {
    stop("Each sample must have a positive total TPM.", call. = FALSE)
  }
  output <- sweep(values, 2, totals, "/") * 1e6
  dimnames(output) <- dimnames(values)
  output
}

annotation_for_expression <- function(expression, annotation, expression_type) {
  if (is.null(annotation)) {
    if (identical(expression_type, "counts")) {
      stop("Count input requires annotation with gene lengths", call. = FALSE)
    }
    return(tibble::tibble(
      gene_id = rownames(expression),
      gene_symbol = rownames(expression),
      mapping_action = "input_identifier_retained"
    ))
  }

  annotation <- validate_annotation(annotation)
  if (identical(expression_type, "counts")) {
    if (!("gene_length_bp" %in% names(annotation))) {
      stop("Count input requires a gene_length_bp annotation column", call. = FALSE)
    }
    annotation$gene_length_bp <- suppressWarnings(as.numeric(annotation$gene_length_bp))
    matching_annotation <- annotation[match(rownames(expression), annotation$gene_id), , drop = FALSE]
    missing_genes <- rownames(expression)[is.na(matching_annotation$gene_id)]
    if (length(missing_genes) > 0L) {
      stop(
        sprintf(
          "gene_length_bp is missing for gene_id: %s",
          paste(missing_genes, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    if (any(!is.finite(matching_annotation$gene_length_bp)) ||
        any(matching_annotation$gene_length_bp <= 0)) {
      stop("Count input requires positive gene_length_bp values", call. = FALSE)
    }
  }

  annotation[match(rownames(expression), annotation$gene_id), , drop = FALSE]
}

make_mapping_report <- function(expression, annotation) {
  if (is.null(annotation)) {
    return(tibble::tibble(
      gene_id = rownames(expression),
      gene_symbol = rownames(expression),
      mapping_action = "input_identifier_retained"
    ))
  }

  symbols <- annotation$gene_symbol
  duplicate_symbols <- !is.na(symbols) & nzchar(symbols) &
    (duplicated(symbols) | duplicated(symbols, fromLast = TRUE))
  tibble::tibble(
    gene_id = annotation$gene_id,
    gene_symbol = dplyr::na_if(annotation$gene_symbol, ""),
    mapping_action = dplyr::case_when(
      is.na(symbols) | !nzchar(symbols) ~ "unmapped",
      duplicate_symbols ~ "duplicate_symbol_collapsed",
      TRUE ~ "mapped"
    )
  ) |>
    dplyr::arrange(.data$gene_id)
}

make_excluded_report <- function(mapping_report, tpm) {
  unmapped <- mapping_report |>
    dplyr::filter(.data$mapping_action == "unmapped") |>
    dplyr::transmute(
      gene_id = .data$gene_id,
      gene_symbol = .data$gene_symbol,
      reason = "unmapped_gene_symbol"
    )
  constant_rows <- apply(tpm, 1L, function(values) {
    length(unique(values)) == 1L
  })
  constants <- tibble::tibble(
    gene_id = rownames(tpm)[constant_rows],
    gene_symbol = rownames(tpm)[constant_rows],
    reason = "constant_expression"
  )

  dplyr::bind_rows(unmapped, constants) |>
    dplyr::arrange(.data$gene_id, .data$reason)
}

prepare_expression <- function(expression, expression_type, annotation = NULL) {
  validate_expression_matrix(expression)
  if (!is.character(expression_type) || length(expression_type) != 1L ||
      is.na(expression_type) || !(expression_type %in% c("counts", "tpm"))) {
    stop("expression_type must be 'counts' or 'tpm'", call. = FALSE)
  }

  expression_annotation <- annotation_for_expression(
    expression,
    annotation,
    expression_type
  )
  if (identical(expression_type, "counts")) {
    gene_lengths <- expression_annotation$gene_length_bp
    names(gene_lengths) <- expression_annotation$gene_id
    tpm <- counts_to_tpm(expression, gene_lengths)
  } else {
    tpm <- expression
  }

  mapping_report <- make_mapping_report(
    expression,
    if (is.null(annotation)) NULL else expression_annotation
  )
  if (is.null(annotation)) {
    tpm <- tpm[sort(rownames(tpm)), , drop = FALSE]
  } else {
    tpm <- collapse_to_symbols(tpm, expression_annotation)
  }
  tpm <- normalise_tpm(tpm)
  excluded_genes <- make_excluded_report(mapping_report, tpm)

  list(
    tpm = tpm,
    log_expression = log2(tpm + 1),
    mapping_report = mapping_report,
    excluded_genes = excluded_genes
  )
}
