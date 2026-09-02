validate_cpm_matrix <- function(cpm) {
  if (!is.matrix(cpm) || !is.numeric(cpm)) {
    stop("cpm must be a numeric matrix", call. = FALSE)
  }
  if (nrow(cpm) == 0L || ncol(cpm) == 0L) {
    stop("cpm must contain at least one gene and one sample", call. = FALSE)
  }
  if (is.null(rownames(cpm)) || anyNA(rownames(cpm))) {
    stop("CPM gene identifiers must be non-missing", call. = FALSE)
  }

  gene_ids <- trimws(rownames(cpm))
  if (any(!nzchar(gene_ids))) {
    stop("CPM gene identifiers must be non-empty", call. = FALSE)
  }
  if (anyDuplicated(gene_ids) > 0L) {
    stop("CPM gene identifiers must be unique after trimming", call. = FALSE)
  }
  if (is.null(colnames(cpm)) || anyNA(colnames(cpm)) ||
      any(!nzchar(colnames(cpm)))) {
    stop("CPM sample identifiers must be non-missing and non-empty", call. = FALSE)
  }
  if (anyDuplicated(colnames(cpm)) > 0L) {
    stop("CPM sample identifiers must be unique", call. = FALSE)
  }
  if (any(!is.finite(cpm))) {
    stop("CPM values must be finite", call. = FALSE)
  }

  rownames(cpm) <- gene_ids
  cpm
}

extract_gtf_attribute <- function(attributes, key) {
  if (!is.character(key) || length(key) != 1L || is.na(key) || !nzchar(key)) {
    stop("key must be one non-missing, non-empty character value", call. = FALSE)
  }

  pattern <- paste0("(?:^|;)[[:space:]]*", key,
    "[[:space:]]+\\\"([^\\\"]+)\\\"")
  stringr::str_match(attributes, pattern)[, 2L]
}

validate_gtf_gene_annotation <- function(annotation) {
  if (!inherits(annotation, "data.frame")) {
    stop("GTF annotation must be a data frame", call. = FALSE)
  }

  required_columns <- c("gene_id", "gene_name", "gene_type")
  missing_columns <- setdiff(required_columns, names(annotation))
  if (length(missing_columns) > 0L) {
    stop(
      sprintf(
        "GTF annotation is missing required columns: %s",
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  annotation <- tibble::as_tibble(annotation) |>
    dplyr::transmute(
      gene_id = trimws(as.character(.data$gene_id)),
      gene_name = trimws(as.character(.data$gene_name)),
      gene_type = trimws(as.character(.data$gene_type))
    )
  if (nrow(annotation) == 0L) {
    stop("GTF annotation must contain at least one gene record", call. = FALSE)
  }
  if (anyNA(annotation$gene_id) || any(!nzchar(annotation$gene_id))) {
    stop("GTF annotation gene_id values must be non-missing and non-empty", call. = FALSE)
  }
  if (anyDuplicated(annotation$gene_id) > 0L) {
    stop("GTF annotation gene_id values must be unique", call. = FALSE)
  }
  if (!any(!is.na(annotation$gene_name) & nzchar(annotation$gene_name))) {
    stop("GTF annotation must contain at least one usable gene_name", call. = FALSE)
  }

  annotation
}

read_gtf_gene_annotation <- function(path, chunk_size = 100000L) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("GTF path must be one non-missing, non-empty character value", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop(sprintf("GTF file does not exist: %s", path), call. = FALSE)
  }
  if (!is.numeric(chunk_size) || length(chunk_size) != 1L || is.na(chunk_size) ||
      !is.finite(chunk_size) || chunk_size < 1L || chunk_size != as.integer(chunk_size)) {
    stop("chunk_size must be a positive integer", call. = FALSE)
  }

  connection <- if (grepl("[.]gz$", path, ignore.case = TRUE)) {
    gzfile(path, "rt")
  } else {
    file(path, "rt")
  }
  on.exit(close(connection), add = TRUE)

  gene_rows <- list()
  chunk_index <- 0L
  repeat {
    lines <- readLines(connection, n = as.integer(chunk_size), warn = FALSE)
    if (length(lines) == 0L) {
      break
    }

    records <- lines[!startsWith(trimws(lines), "#")]
    if (length(records) == 0L) {
      next
    }
    fields_by_record <- strsplit(records, "\t", fixed = TRUE)
    field_counts <- lengths(fields_by_record)
    if (any(field_counts != 9L)) {
      stop("Each non-comment GTF record must contain nine tab-separated fields", call. = FALSE)
    }
    fields <- do.call(rbind, fields_by_record)
    is_gene <- fields[, 3L] == "gene"
    if (!any(is_gene)) {
      next
    }

    gene_fields <- fields[is_gene, , drop = FALSE]
    chunk_index <- chunk_index + 1L
    gene_rows[[chunk_index]] <- tibble::tibble(
      gene_id = extract_gtf_attribute(gene_fields[, 9L], "gene_id"),
      gene_name = extract_gtf_attribute(gene_fields[, 9L], "gene_name"),
      gene_type = dplyr::coalesce(
        extract_gtf_attribute(gene_fields[, 9L], "gene_type"),
        extract_gtf_attribute(gene_fields[, 9L], "gene_biotype")
      )
    )
  }

  if (length(gene_rows) == 0L) {
    stop("GTF annotation must contain at least one gene record", call. = FALSE)
  }
  validate_gtf_gene_annotation(dplyr::bind_rows(gene_rows))
}

collapse_cpm_to_gene_names <- function(cpm, annotation) {
  cpm <- validate_cpm_matrix(cpm)
  annotation <- validate_gtf_gene_annotation(annotation)
  annotation_index <- match(rownames(cpm), annotation$gene_id)
  mapped_gene_names <- annotation$gene_name[annotation_index]
  usable <- !is.na(mapped_gene_names) & nzchar(mapped_gene_names)
  if (!any(usable)) {
    return(matrix(
      numeric(),
      nrow = 0L,
      ncol = ncol(cpm),
      dimnames = list(character(), colnames(cpm))
    ))
  }

  mapped_gene_names <- mapped_gene_names[usable]
  gene_name_order <- unique(mapped_gene_names)
  group_index <- match(mapped_gene_names, gene_name_order)
  output <- rowsum(
    cpm[usable, , drop = FALSE],
    group = group_index,
    reorder = FALSE
  )
  storage.mode(output) <- "double"
  rownames(output) <- gene_name_order
  colnames(output) <- colnames(cpm)
  output
}

make_cpm_mapping_report <- function(cpm, annotation) {
  matching_indices <- match(rownames(cpm), annotation$gene_id)
  gene_names <- annotation$gene_name[matching_indices]
  usable_gene_names <- !is.na(gene_names) & nzchar(gene_names)
  duplicate_gene_names <- usable_gene_names &
    (duplicated(gene_names) | duplicated(gene_names, fromLast = TRUE))

  tibble::tibble(
    gene_id = rownames(cpm),
    gene_name = gene_names,
    mapping_action = dplyr::case_when(
      is.na(matching_indices) ~ "missing_gtf_gene_id",
      !usable_gene_names ~ "missing_gene_name",
      duplicate_gene_names ~ "duplicate_gene_name_aggregated_for_dtangle",
      TRUE ~ "mapped"
    )
  )
}

make_tca_expression <- function(expression) {
  if (!is.list(expression) || !all(c("coordinates", "cpm") %in% names(expression))) {
    stop("expression must contain coordinates and cpm", call. = FALSE)
  }
  cpm <- validate_cpm_matrix(expression$cpm)
  if (!identical(rownames(cpm), expression$coordinates$gene_id)) {
    stop("Expression genes and BED coordinates must match exactly", call. = FALSE)
  }
  log2(cpm)
}

make_dtangle_expression <- function(expression, annotation) {
  if (!is.list(expression) || !all(c("coordinates", "cpm") %in% names(expression))) {
    stop("expression must contain coordinates and cpm", call. = FALSE)
  }
  cpm <- validate_cpm_matrix(expression$cpm)
  annotation <- validate_gtf_gene_annotation(annotation)
  mapping_report <- make_cpm_mapping_report(cpm, annotation)
  symbol_cpm <- collapse_cpm_to_gene_names(cpm, annotation)
  if (nrow(symbol_cpm) == 0L || any(symbol_cpm <= 0)) {
    stop("No positive CPM genes map to a usable gene symbol", call. = FALSE)
  }
  list(log_expression = log2(symbol_cpm), mapping_report = mapping_report)
}

make_excluded_genes <- function(mapping_report) {
  mapping_report |>
    dplyr::filter(!(.data$mapping_action %in% c(
      "mapped", "duplicate_gene_name_aggregated_for_dtangle"
    ))) |>
    dplyr::transmute(
      gene_id = .data$gene_id,
      gene_name = .data$gene_name,
      reason = .data$mapping_action
    )
}
