tensor_utc_time <- function() {
  format(Sys.time(), tz = "UTC", usetz = TRUE)
}

append_tensor_log <- function(log_file, message_text) {
  if (is.null(log_file)) {
    return(invisible(NULL))
  }
  if (!is.character(log_file) || length(log_file) != 1L ||
      is.na(log_file) || !nzchar(log_file)) {
    stop("log_file must be NULL or one non-empty path", call. = FALSE)
  }
  cat(
    sprintf("utc_time=%s %s\n", tensor_utc_time(), message_text),
    file = log_file,
    append = TRUE
  )
  invisible(log_file)
}

validate_tensor_positive_integer <- function(value, label) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < 1 || value != as.integer(value)) {
    stop(sprintf("%s must be one positive integer", label), call. = FALSE)
  }
  as.integer(value)
}

validate_tensor_tca_version <- function() {
  required_version <- "1.2.1"
  if (!requireNamespace("TCA", quietly = TRUE)) {
    stop("The TCA package is required for tensor extraction", call. = FALSE)
  }
  observed_version <- as.character(utils::packageVersion("TCA"))
  if (!identical(observed_version, required_version)) {
    stop(
      sprintf(
        "TCA version %s is required; found %s",
        required_version,
        observed_version
      ),
      call. = FALSE
    )
  }
  observed_version
}

validate_tensor_matrix <- function(x, label) {
  if (!is.matrix(x) || !is.numeric(x) || nrow(x) == 0L || ncol(x) == 0L) {
    stop(sprintf("%s must be a non-empty numeric matrix", label), call. = FALSE)
  }
  row_ids <- rownames(x)
  column_ids <- colnames(x)
  if (is.null(row_ids) || anyNA(row_ids) || any(!nzchar(row_ids)) ||
      anyDuplicated(row_ids) > 0L) {
    stop(sprintf("%s row identifiers must be unique and non-empty", label),
      call. = FALSE
    )
  }
  if (is.null(column_ids) || anyNA(column_ids) || any(!nzchar(column_ids)) ||
      anyDuplicated(column_ids) > 0L) {
    stop(sprintf("%s column identifiers must be unique and non-empty", label),
      call. = FALSE
    )
  }
  if (any(!is.finite(x))) {
    stop(sprintf("%s values must be finite", label), call. = FALSE)
  }
  invisible(TRUE)
}

validate_tensor_list <- function(tensor) {
  if (!is.list(tensor) || length(tensor) == 0L) {
    stop("tensor must be a non-empty list", call. = FALSE)
  }
  source_names <- names(tensor)
  if (is.null(source_names) || anyNA(source_names) ||
      any(!nzchar(source_names)) || anyDuplicated(source_names) > 0L) {
    stop("tensor source names must be unique and non-empty", call. = FALSE)
  }
  purrr::iwalk(tensor, function(source_matrix, source_name) {
    validate_tensor_matrix(source_matrix, sprintf("Tensor source '%s'", source_name))
  })
  reference_dimensions <- dim(tensor[[1L]])
  reference_dimnames <- dimnames(tensor[[1L]])
  compatible <- purrr::map_lgl(tensor, function(source_matrix) {
    identical(dim(source_matrix), reference_dimensions) &&
      identical(dimnames(source_matrix), reference_dimnames)
  })
  if (!all(compatible)) {
    stop("All tensor sources must have identical dimensions and order", call. = FALSE)
  }
  slugs <- slugify_cell_group(source_names)
  if (any(!nzchar(slugs)) || anyDuplicated(slugs) > 0L) {
    stop("Tensor source names must have unique non-empty slugs", call. = FALSE)
  }
  invisible(TRUE)
}

validate_tensor_contract <- function(
    tensor,
    gene_ids,
    sample_ids,
    cell_groups) {
  if (!is.character(cell_groups) || length(cell_groups) == 0L ||
      anyNA(cell_groups) || any(!nzchar(cell_groups)) ||
      anyDuplicated(cell_groups) > 0L) {
    stop("Tensor cell groups must be unique and non-empty", call. = FALSE)
  }
  if (!is.list(tensor) || is.null(names(tensor)) ||
      anyNA(names(tensor)) || any(!nzchar(names(tensor))) ||
      !identical(names(tensor), cell_groups)) {
    stop("Tensor cell groups must match the TCA weights exactly", call. = FALSE)
  }
  slugs <- slugify_cell_group(cell_groups)
  if (any(!nzchar(slugs)) || anyDuplicated(slugs) > 0L) {
    stop("Tensor cell groups must have unique non-empty slugs", call. = FALSE)
  }
  valid <- purrr::map_lgl(tensor, function(source_matrix) {
    is.matrix(source_matrix) && is.numeric(source_matrix) &&
      identical(rownames(source_matrix), gene_ids) &&
      identical(colnames(source_matrix), sample_ids) &&
      identical(dim(source_matrix), c(length(gene_ids), length(sample_ids))) &&
      all(is.finite(source_matrix))
  })
  if (!all(valid)) {
    stop("Tensor genes, samples, dimensions, and values must match", call. = FALSE)
  }
  invisible(TRUE)
}

count_excluded_constant_genes <- function(coordinates, modeled_gene_ids) {
  if (!inherits(coordinates, "data.frame") ||
      !"gene_id" %in% names(coordinates) ||
      anyNA(coordinates$gene_id) || any(!nzchar(coordinates$gene_id)) ||
      anyDuplicated(coordinates$gene_id) > 0L) {
    stop("coordinates must contain unique non-empty gene_id values", call. = FALSE)
  }
  if (!is.character(modeled_gene_ids) || length(modeled_gene_ids) == 0L ||
      anyNA(modeled_gene_ids) || any(!nzchar(modeled_gene_ids)) ||
      anyDuplicated(modeled_gene_ids) > 0L ||
      !all(modeled_gene_ids %in% coordinates$gene_id)) {
    stop(
      "Every modeled gene must occur once in the prepared coordinates",
      call. = FALSE
    )
  }
  as.integer(nrow(coordinates) - length(modeled_gene_ids))
}

extract_full_tensor <- function(X, model, num_cores = 1L, log_file = NULL) {
  validate_tensor_matrix(X, "TCA expression")
  num_cores <- validate_tensor_positive_integer(num_cores, "num_cores")
  if (!is.list(model) || !is.matrix(model$W) ||
      !identical(colnames(X), rownames(model$W))) {
    stop("Model sample order must match expression sample order exactly", call. = FALSE)
  }
  tca_version <- validate_tensor_tca_version()
  append_tensor_log(
    log_file,
    sprintf(
      "stage=tensor_extract event=extract_start genes=%d samples=%d sources=%d num_cores=%d scale=log2_cpm tca_version=%s",
      nrow(X), ncol(X), ncol(model$W), num_cores, tca_version
    )
  )
  tensor <- tryCatch(
    TCA::tensor(
      X = X,
      tca.mdl = model,
      scale = FALSE,
      parallel = num_cores > 1L,
      num_cores = num_cores,
      log_file = log_file,
      verbose = TRUE
    ),
    error = function(error) {
      append_tensor_log(
        log_file,
        sprintf("stage=tensor_extract event=extract_error message=%s", conditionMessage(error))
      )
      stop(error)
    }
  )
  if (!is.list(tensor) || length(tensor) != ncol(model$W)) {
    stop("TCA tensor source dimension does not match the model", call. = FALSE)
  }
  names(tensor) <- colnames(model$W)
  tensor <- purrr::map(tensor, function(source_matrix) {
    if (!is.matrix(source_matrix) || !identical(dim(source_matrix), dim(X))) {
      stop("TCA tensor returned an invalid gene-by-sample dimension", call. = FALSE)
    }
    dimnames(source_matrix) <- dimnames(X)
    source_matrix
  })
  validate_tensor_contract(tensor, rownames(X), colnames(X), colnames(model$W))
  append_tensor_log(
    log_file,
    "stage=tensor_extract event=extract_complete scale=log2_cpm"
  )
  tensor
}

write_cell_type_beds <- function(tensor, coordinates, output_dir) {
  if (!inherits(coordinates, "data.frame") ||
      !identical(names(coordinates), c("#chr", "start", "end", "gene_id"))) {
    stop("coordinates must have the exact BED columns", call. = FALSE)
  }
  if (!is.character(output_dir) || length(output_dir) != 1L ||
      is.na(output_dir) || !nzchar(output_dir)) {
    stop("output_dir must be one non-empty path", call. = FALSE)
  }
  validate_tensor_contract(
    tensor,
    coordinates$gene_id,
    colnames(tensor[[1L]]),
    names(tensor)
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- stats::setNames(
    file.path(output_dir, paste0(slugify_cell_group(names(tensor)), ".bed.gz")),
    names(tensor)
  )
  purrr::iwalk(paths, function(path, cell_group) {
    write_expression_bed(path, coordinates, tensor[[cell_group]])
  })
  inventory <- tibble::tibble(
    logical_name = paste0(slugify_cell_group(names(paths)), "_expression"),
    path = basename(unname(paths)),
    n_genes = nrow(tensor[[1L]]),
    n_samples = ncol(tensor[[1L]]),
    scale = "log2_cpm",
    cell_group = names(paths),
    slug = slugify_cell_group(names(paths))
  )
  list(paths = paths, inventory = inventory)
}
