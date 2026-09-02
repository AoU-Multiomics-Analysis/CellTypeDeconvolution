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

extract_tensor_shard <- function(
    X,
    model,
    genes,
    num_cores = 1L,
    log_file = NULL) {
  validate_tensor_matrix(X, "TCA expression")
  num_cores <- validate_tensor_positive_integer(num_cores, "num_cores")
  if (!is.list(model) || !is.matrix(model$W) ||
      is.null(rownames(model$W)) || is.null(colnames(model$W))) {
    stop("model must be a fitted TCA model with named weights", call. = FALSE)
  }
  if (!identical(colnames(X), rownames(model$W))) {
    stop("Model sample order must match expression sample order exactly", call. = FALSE)
  }
  if (!is.character(genes) || length(genes) == 0L || anyNA(genes) ||
      any(!nzchar(genes)) || anyDuplicated(genes) > 0L) {
    stop("genes must be a non-empty vector of unique identifiers", call. = FALSE)
  }
  if (!all(genes %in% rownames(X))) {
    stop("Every shard gene must be present in expression", call. = FALSE)
  }
  model_genes <- rownames(model$mus_hat)
  if (is.null(model_genes) || !all(genes %in% model_genes)) {
    stop("Every shard gene must be present in the fitted model", call. = FALSE)
  }
  tca_version <- validate_tensor_tca_version()
  append_tensor_log(
    log_file,
    sprintf(
      paste0(
        "stage=tensor_extract event=extract_start genes=%d samples=%d ",
        "sources=%d num_cores=%d scale=log2_cpm tca_version=%s"
      ),
      length(genes),
      ncol(X),
      ncol(model$W),
      num_cores,
      tca_version
    )
  )

  tensor <- tryCatch({
    subset_model <- TCA::tcasub(
      model,
      genes,
      log_file = log_file,
      verbose = TRUE
    )
    TCA::tensor(
      X = X[genes, , drop = FALSE],
      tca.mdl = subset_model,
      scale = FALSE,
      parallel = num_cores > 1L,
      num_cores = num_cores,
      log_file = log_file,
      verbose = TRUE
    )
  }, error = function(error) {
    append_tensor_log(
      log_file,
      sprintf(
        "stage=tensor_extract event=extract_error message=%s",
        conditionMessage(error)
      )
    )
    stop(error)
  })

  if (!is.list(tensor) || length(tensor) != ncol(model$W)) {
    stop("TCA tensor source dimension does not match the model", call. = FALSE)
  }
  names(tensor) <- colnames(model$W)
  tensor <- purrr::map(tensor, function(source_matrix) {
    if (!is.matrix(source_matrix) ||
        !identical(dim(source_matrix), c(length(genes), ncol(X)))) {
      stop("TCA tensor returned an invalid feature-by-sample dimension",
        call. = FALSE
      )
    }
    rownames(source_matrix) <- genes
    colnames(source_matrix) <- colnames(X)
    source_matrix
  })
  validate_tensor_list(tensor)
  append_tensor_log(
    log_file,
    "stage=tensor_extract event=extract_complete scale=log2_cpm"
  )
  tensor
}

write_tensor_shard <- function(path, tensor, shard_id) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("path must be one non-empty file path", call. = FALSE)
  }
  shard_id <- validate_tensor_positive_integer(shard_id, "shard_id")
  validate_tensor_list(tensor)
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("The hdf5r package is required for tensor shard output", call. = FALSE)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  h5 <- hdf5r::H5File$new(path, mode = "w")
  on.exit(h5$close_all(), add = TRUE)

  gene_names <- rownames(tensor[[1L]])
  sample_ids <- colnames(tensor[[1L]])
  source_names <- names(tensor)
  h5[["gene_name"]] <- enc2utf8(gene_names)
  h5[["sample_id"]] <- enc2utf8(sample_ids)
  hdf5r::h5attr(h5, "shard_id") <- shard_id
  hdf5r::h5attr(h5, "source_names") <- enc2utf8(source_names)
  hdf5r::h5attr(h5, "scale") <- "log2_cpm"
  sources_group <- h5$create_group("sources")

  purrr::iwalk(tensor, function(source_matrix, source_name) {
    source_group <- sources_group$create_group(slugify_cell_group(source_name))
    source_group$create_dataset(
      "expression",
      robj = source_matrix,
      chunk_dims = c(
        min(500L, nrow(source_matrix)),
        min(256L, ncol(source_matrix))
      ),
      gzip_level = 6L
    )
    hdf5r::h5attr(source_group, "cell_group") <- enc2utf8(source_name)
  })
  invisible(path)
}

read_tensor_shard <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Tensor shard does not exist: %s", path), call. = FALSE)
  }
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("The hdf5r package is required for tensor shard input", call. = FALSE)
  }
  h5 <- hdf5r::H5File$new(path, mode = "r")
  on.exit(h5$close_all(), add = TRUE)
  required_datasets <- c("gene_name", "sample_id", "sources")
  if (!all(required_datasets %in% names(h5))) {
    stop("Tensor shard has an invalid dimension or missing dataset", call. = FALSE)
  }
  required_attributes <- c("shard_id", "source_names", "scale")
  if (!all(required_attributes %in% hdf5r::h5attr_names(h5))) {
    stop("Tensor shard is missing required attributes", call. = FALSE)
  }
  scale <- h5$attr_open("scale")$read()
  if (!identical(scale, "log2_cpm")) {
    stop("Tensor shard scale must be log2_cpm", call. = FALSE)
  }
  shard_id <- as.integer(h5$attr_open("shard_id")$read())
  source_names <- as.character(h5$attr_open("source_names")$read())
  gene_names <- as.character(h5[["gene_name"]][])
  sample_ids <- as.character(h5[["sample_id"]][])
  source_slugs <- slugify_cell_group(source_names)
  sources_group <- h5[["sources"]]
  if (!setequal(names(sources_group), source_slugs)) {
    stop("Tensor shard sources do not match source_names", call. = FALSE)
  }
  tensor <- purrr::map2(source_names, source_slugs, function(source_name, slug) {
    source_path <- paste0("sources/", slug)
    if (!("expression" %in% names(h5[[source_path]]))) {
      stop("Tensor shard has an invalid source dimension", call. = FALSE)
    }
    source_matrix <- h5[[paste0(source_path, "/expression")]]$read(
      drop = FALSE
    )
    if (!identical(dim(source_matrix), c(length(gene_names), length(sample_ids)))) {
      stop("Tensor shard has an invalid source dimension", call. = FALSE)
    }
    rownames(source_matrix) <- gene_names
    colnames(source_matrix) <- sample_ids
    source_matrix
  })
  names(tensor) <- source_names
  validate_tensor_list(tensor)
  list(shard_id = shard_id, tensor = tensor)
}

validate_shard_manifest <- function(manifest) {
  if (!inherits(manifest, "data.frame")) {
    stop("manifest must be a data frame", call. = FALSE)
  }
  required_columns <- c(
    "gene_index", "gene_name", "shard_id", "shard_name",
    "index_within_shard"
  )
  missing_columns <- setdiff(required_columns, names(manifest))
  if (length(missing_columns) > 0L) {
    stop(
      sprintf(
        "manifest is missing required columns: %s",
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  manifest <- tibble::as_tibble(manifest)
  if (nrow(manifest) == 0L) {
    stop("manifest must contain at least one gene", call. = FALSE)
  }
  if (anyNA(manifest$gene_name) || any(!nzchar(manifest$gene_name))) {
    stop("manifest gene identifiers must be non-missing and non-empty",
      call. = FALSE
    )
  }
  if (anyDuplicated(manifest$gene_name) > 0L) {
    stop("manifest contains duplicate gene identifiers", call. = FALSE)
  }
  expected_gene_index <- as.integer(seq_len(nrow(manifest)))
  if (!identical(as.integer(manifest$gene_index), expected_gene_index)) {
    stop("manifest gene order must follow gene_index exactly", call. = FALSE)
  }
  shard_ids <- as.integer(manifest$shard_id)
  if (anyNA(shard_ids) || any(shard_ids < 1L) ||
      !identical(sort(unique(shard_ids)), seq_len(max(shard_ids)))) {
    stop("manifest shard identifiers must be consecutive positive integers",
      call. = FALSE
    )
  }
  expected_shard_names <- sprintf("shard_%05d", shard_ids)
  if (!identical(as.character(manifest$shard_name), expected_shard_names)) {
    stop("manifest shard_name values do not match shard identifiers",
      call. = FALSE
    )
  }
  expected_within <- ave(
    seq_len(nrow(manifest)),
    shard_ids,
    FUN = seq_along
  ) |> as.integer()
  if (!identical(as.integer(manifest$index_within_shard), expected_within)) {
    stop("manifest index_within_shard values are misordered", call. = FALSE)
  }
  manifest$shard_id <- shard_ids
  manifest
}

write_tensor_tsv_rows <- function(connection, genes, source_matrix, header) {
  rows <- tibble::as_tibble(source_matrix, .name_repair = "minimal") |>
    dplyr::mutate(gene_name = genes, .before = 1L)
  utils::write.table(
    rows,
    file = connection,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = header,
    append = !header
  )
  invisible(connection)
}

assemble_hdf5_shards <- function(
    shard_paths,
    manifest,
    output_dir,
    pipeline_version,
    tca_version,
    write_tsv = FALSE) {
  manifest <- validate_shard_manifest(manifest)
  if (!is.character(shard_paths)) {
    stop("shard_paths must be a character vector", call. = FALSE)
  }
  expected_shard_ids <- sort(unique(manifest$shard_id))
  if (length(shard_paths) != length(expected_shard_ids) ||
      any(!file.exists(shard_paths))) {
    stop("missing tensor shard files", call. = FALSE)
  }
  if (!is.character(output_dir) || length(output_dir) != 1L ||
      is.na(output_dir) || !nzchar(output_dir)) {
    stop("output_dir must be one non-empty path", call. = FALSE)
  }
  if (!is.character(pipeline_version) || length(pipeline_version) != 1L ||
      is.na(pipeline_version) || !nzchar(pipeline_version)) {
    stop("pipeline_version must be one non-empty value", call. = FALSE)
  }
  if (!is.character(tca_version) || length(tca_version) != 1L ||
      is.na(tca_version) || !nzchar(tca_version)) {
    stop("tca_version must be one non-empty value", call. = FALSE)
  }
  if (!is.logical(write_tsv) || length(write_tsv) != 1L || is.na(write_tsv)) {
    stop("write_tsv must be one non-missing logical value", call. = FALSE)
  }

  observed_shard_ids <- integer(length(shard_paths))
  source_names <- NULL
  sample_ids <- NULL
  for (path_index in seq_along(shard_paths)) {
    shard <- read_tensor_shard(shard_paths[[path_index]])
    shard_id <- shard$shard_id
    observed_shard_ids[[path_index]] <- shard_id
    if (!(shard_id %in% expected_shard_ids)) {
      stop("missing or unexpected tensor shard identifier", call. = FALSE)
    }
    if (is.null(source_names)) {
      source_names <- names(shard$tensor)
      sample_ids <- colnames(shard$tensor[[1L]])
    } else {
      if (!identical(names(shard$tensor), source_names)) {
        stop("Tensor shard source order does not match", call. = FALSE)
      }
      if (!identical(colnames(shard$tensor[[1L]]), sample_ids)) {
        stop("Tensor shard sample order does not match", call. = FALSE)
      }
    }
    expected_genes <- manifest |>
      dplyr::filter(.data$shard_id == .env$shard_id) |>
      dplyr::arrange(.data$index_within_shard) |>
      dplyr::pull("gene_name")
    observed_genes <- rownames(shard$tensor[[1L]])
    if (length(observed_genes) != length(expected_genes)) {
      stop("Tensor shard dimension does not match the manifest", call. = FALSE)
    }
    if (!identical(observed_genes, expected_genes)) {
      stop("Tensor shard gene order does not match the manifest", call. = FALSE)
    }
  }
  if (anyDuplicated(observed_shard_ids) > 0L) {
    stop("duplicate tensor shard identifiers", call. = FALSE)
  }
  if (!setequal(observed_shard_ids, expected_shard_ids)) {
    stop("missing tensor shard identifiers", call. = FALSE)
  }
  if (!identical(observed_shard_ids, expected_shard_ids)) {
    stop("Tensor shard order does not match the manifest", call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  hdf5_paths <- stats::setNames(
    file.path(output_dir, paste0(slugify_cell_group(source_names), ".h5")),
    source_names
  )
  h5_files <- purrr::map(hdf5_paths, function(path) {
    hdf5r::H5File$new(path, mode = "w")
  })
  on.exit(purrr::walk(h5_files, ~ .x$close_all()), add = TRUE)
  expression_datasets <- purrr::map2(
    h5_files,
    source_names,
    function(h5, source_name) {
      h5[["gene_name"]] <- enc2utf8(manifest$gene_name)
      h5[["sample_id"]] <- enc2utf8(sample_ids)
      hdf5r::h5attr(h5, "cell_group") <- enc2utf8(source_name)
      hdf5r::h5attr(h5, "scale") <- "log2_cpm"
      hdf5r::h5attr(h5, "pipeline_version") <- enc2utf8(pipeline_version)
      hdf5r::h5attr(h5, "tca_version") <- enc2utf8(tca_version)
      h5$create_dataset(
        "expression",
        dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE,
        dims = c(nrow(manifest), length(sample_ids)),
        chunk_dims = c(
          min(500L, nrow(manifest)),
          min(256L, length(sample_ids))
        ),
        gzip_level = 6L
      )
    }
  )

  tsv_paths <- stats::setNames(character(), character())
  tsv_connections <- NULL
  if (write_tsv) {
    tsv_paths <- stats::setNames(
      file.path(
        output_dir,
        paste0(slugify_cell_group(source_names), ".tsv.gz")
      ),
      source_names
    )
    tsv_connections <- purrr::map(tsv_paths, ~ gzfile(.x, open = "wt"))
    on.exit(purrr::walk(tsv_connections, close), add = TRUE)
  }

  for (path_index in seq_along(shard_paths)) {
    shard <- read_tensor_shard(shard_paths[[path_index]])
    shard_id <- expected_shard_ids[[path_index]]
    positions <- which(manifest$shard_id == shard_id)
    purrr::iwalk(shard$tensor, function(source_matrix, source_name) {
      expression_datasets[[source_name]][positions, ] <- source_matrix
      if (write_tsv) {
        write_tensor_tsv_rows(
          tsv_connections[[source_name]],
          rownames(source_matrix),
          source_matrix,
          header = shard_id == expected_shard_ids[[1L]]
        )
      }
    })
  }
  purrr::walk(h5_files, ~ .x$flush())
  attr(hdf5_paths, "tsv_paths") <- tsv_paths
  hdf5_paths
}
