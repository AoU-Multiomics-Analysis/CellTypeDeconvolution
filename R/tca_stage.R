validate_matrix_identifiers <- function(x, row_label, column_label) {
  row_ids <- rownames(x)
  column_ids <- colnames(x)
  if (is.null(row_ids) || anyNA(row_ids) || any(!nzchar(row_ids)) ||
      anyDuplicated(row_ids) > 0L) {
    stop(
      sprintf("%s identifiers must be unique and non-empty", row_label),
      call. = FALSE
    )
  }
  if (is.null(column_ids) || anyNA(column_ids) || any(!nzchar(column_ids)) ||
      anyDuplicated(column_ids) > 0L) {
    stop(
      sprintf("%s identifiers must be unique and non-empty", column_label),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

validate_tca_inputs <- function(X, W, C2 = NULL) {
  if (!is.matrix(X) || !is.numeric(X) || nrow(X) == 0L || ncol(X) == 0L) {
    stop("TCA expression must be a non-empty numeric matrix", call. = FALSE)
  }
  validate_matrix_identifiers(X, "Expression gene", "Expression sample")
  if (any(!is.finite(X))) {
    stop("TCA expression values must be finite", call. = FALSE)
  }

  if (!is.matrix(W) || !is.numeric(W) || nrow(W) == 0L || ncol(W) < 2L) {
    stop("TCA weights must be a numeric matrix with at least two groups", call. = FALSE)
  }
  validate_matrix_identifiers(W, "Weight sample", "Weight group")
  if (!identical(colnames(X), rownames(W))) {
    stop(
      "TCA weight sample order must match expression sample order exactly",
      call. = FALSE
    )
  }
  if (any(!is.finite(W))) {
    stop("TCA weights must be finite", call. = FALSE)
  }
  if (any(W <= 0)) {
    stop("TCA weights must be strictly positive", call. = FALSE)
  }
  if (any(abs(rowSums(W) - 1) > 1e-8)) {
    stop("TCA weight rows must sum to one within 1e-8", call. = FALSE)
  }

  if (!is.null(C2)) {
    if (!is.matrix(C2) || !is.numeric(C2) || nrow(C2) == 0L || ncol(C2) == 0L) {
      stop("TCA covariates must be a non-empty numeric matrix", call. = FALSE)
    }
    validate_matrix_identifiers(C2, "Covariate sample", "Covariate")
    if (!identical(colnames(X), rownames(C2))) {
      stop(
        "TCA covariate sample order must match expression sample order exactly",
        call. = FALSE
      )
    }
    if (anyNA(C2)) {
      stop("TCA input contains missing covariates", call. = FALSE)
    }
    if (any(!is.finite(C2))) {
      stop("TCA covariates must be finite", call. = FALSE)
    }
    constant_columns <- apply(
      C2,
      2L,
      function(values) length(unique(values)) == 1L
    )
    intercept_names <- tolower(colnames(C2)) %in% c("intercept", "(intercept)")
    if (any(constant_columns) || any(intercept_names)) {
      stop("TCA covariates must not contain an intercept column", call. = FALSE)
    }
  }

  invisible(TRUE)
}

remove_constant_features <- function(X) {
  if (!is.matrix(X) || !is.numeric(X) || nrow(X) == 0L || ncol(X) == 0L) {
    stop("TCA expression must be a non-empty numeric matrix", call. = FALSE)
  }
  gene_ids <- rownames(X)
  if (is.null(gene_ids) || anyNA(gene_ids) || any(!nzchar(gene_ids)) ||
      anyDuplicated(gene_ids) > 0L) {
    stop("Expression gene identifiers must be unique and non-empty", call. = FALSE)
  }
  constant <- apply(
    X,
    1L,
    function(values) length(unique(values)) == 1L
  )

  list(
    matrix = X[!constant, , drop = FALSE],
    report = tibble::tibble(
      gene_name = rownames(X)[constant],
      reason = rep("constant_expression", sum(constant))
    )
  )
}

validate_positive_integer <- function(value, label) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < 1 || value != as.integer(value)) {
    stop(sprintf("%s must be one positive integer", label), call. = FALSE)
  }
  as.integer(value)
}

make_gene_shard_manifest <- function(
    genes,
    shard_size = pipeline_defaults()$tensor_shard_size) {
  if (!is.character(genes) || length(genes) == 0L || anyNA(genes) ||
      any(!nzchar(genes)) || anyDuplicated(genes) > 0L) {
    stop("genes must be a non-empty vector of unique identifiers", call. = FALSE)
  }
  shard_size <- validate_positive_integer(shard_size, "shard_size")
  gene_index <- seq_along(genes)
  shard_id <- as.integer((gene_index - 1L) %/% shard_size + 1L)
  index_within_shard <- as.integer((gene_index - 1L) %% shard_size + 1L)

  tibble::tibble(
    gene_index = as.integer(gene_index),
    gene_name = genes,
    shard_id = shard_id,
    shard_name = sprintf("shard_%05d", shard_id),
    index_within_shard = index_within_shard
  )
}

write_gene_shards <- function(
    genes,
    shard_size = pipeline_defaults()$tensor_shard_size,
    output_dir) {
  if (!is.character(output_dir) || length(output_dir) != 1L ||
      is.na(output_dir) || !nzchar(output_dir)) {
    stop("output_dir must be one non-empty path", call. = FALSE)
  }
  manifest <- make_gene_shard_manifest(genes, shard_size)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  shard_names <- unique(manifest$shard_name)
  purrr::walk(shard_names, function(shard_name) {
    shard_genes <- manifest$gene_name[manifest$shard_name == shard_name]
    writeLines(shard_genes, file.path(output_dir, paste0(shard_name, ".txt")))
  })

  manifest
}

tca_utc_time <- function() {
  format(Sys.time(), tz = "UTC", usetz = TRUE)
}

append_tca_log <- function(log_file, message_text) {
  if (!is.character(log_file) || length(log_file) != 1L ||
      is.na(log_file) || !nzchar(log_file)) {
    stop("log_file must be one non-empty path", call. = FALSE)
  }
  cat(
    sprintf("utc_time=%s %s\n", tca_utc_time(), message_text),
    file = log_file,
    append = TRUE
  )
  invisible(log_file)
}

validate_tca_version <- function() {
  required_version <- "1.2.1"
  if (!requireNamespace("TCA", quietly = TRUE)) {
    stop("The TCA package is required for model fitting", call. = FALSE)
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

fit_tca_stage <- function(
    X,
    W,
    C2 = NULL,
    num_cores = 1L,
    max_iters = 10L,
    random_seed = 20260901L,
    log_file) {
  validate_tca_inputs(X, W, C2)
  num_cores <- validate_positive_integer(num_cores, "num_cores")
  max_iters <- validate_positive_integer(max_iters, "max_iters")
  random_seed <- validate_positive_integer(random_seed, "random_seed")
  filtered <- remove_constant_features(X)
  if (nrow(filtered$matrix) == 0L) {
    stop("No variable genes remain after constant-gene removal", call. = FALSE)
  }
  X <- filtered$matrix
  tca_version <- validate_tca_version()

  append_tca_log(
    log_file,
    sprintf(
      paste0(
        "stage=tca event=fit_start scale=log2_cpm genes=%d samples=%d ",
        "groups=%d covariates=%d excluded_constant_genes=%d ",
        "num_cores=%d max_iters=%d random_seed=%d tca_version=%s"
      ),
      nrow(X),
      ncol(X),
      ncol(W),
      if (is.null(C2)) 0L else ncol(C2),
      nrow(filtered$report),
      num_cores,
      max_iters,
      random_seed,
      tca_version
    )
  )

  set.seed(random_seed)
  model <- tryCatch(
    TCA::tca(
      X = X,
      W = W,
      C2 = C2,
      refit_W = FALSE,
      vars.mle = FALSE,
      constrain_mu = FALSE,
      parallel = num_cores > 1L,
      num_cores = num_cores,
      max_iters = max_iters,
      log_file = log_file,
      verbose = TRUE
    ),
    error = function(error) {
      append_tca_log(
        log_file,
        sprintf("stage=tca event=fit_error message=%s", conditionMessage(error))
      )
      stop(error)
    }
  )
  append_tca_log(log_file, "stage=tca event=fit_complete")

  list(
    model = model,
    X = X,
    excluded_genes = filtered$report
  )
}
