#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "bootstrap.R"))

extract_log_path <- NULL

run_extract_tca_shard <- function() {
  message(sprintf("stage=extract utc_start=%s", tensor_utc_time()))
  option_list <- list(
    optparse::make_option(
      "--expression-log",
      dest = "expression_log",
      type = "character",
      help = "Gene-by-sample log2(CPM) TSV with gene_name first."
    ),
    optparse::make_option(
      "--model",
      type = "character",
      help = "Cohort-wide fitted TCA model RDS."
    ),
    optparse::make_option(
      "--genes",
      type = "character",
      help = "Ordered newline-delimited gene shard."
    ),
    optparse::make_option(
      "--shard-id",
      dest = "shard_id",
      type = "integer",
      help = "Positive shard identifier."
    ),
    optparse::make_option(
      "--num-cores",
      dest = "num_cores",
      type = "integer",
      default = 1L,
      help = "Number of cores for TCA tensor extraction."
    ),
    optparse::make_option(
      "--output",
      type = "character",
      help = "Output shard HDF5 path."
    ),
    optparse::make_option(
      "--log-file",
      dest = "log_file",
      type = "character",
      default = NULL,
      help = "Extraction log path. Defaults to <output>.log."
    )
  )
  options <- optparse::parse_args(optparse::OptionParser(
    option_list = option_list
  ))
  required_options <- c(
    "expression_log", "model", "genes", "shard_id", "output"
  )
  missing_options <- required_options[vapply(
    options[required_options],
    function(value) {
      is.null(value) || length(value) != 1L || is.na(value) ||
        (is.character(value) && !nzchar(value))
    },
    logical(1)
  )]
  if (length(missing_options) > 0L) {
    stop(
      sprintf(
        "Missing required options: %s",
        paste(missing_options, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  dir.create(dirname(options$output), recursive = TRUE, showWarnings = FALSE)
  extract_log_path <<- if (is.null(options$log_file)) {
    paste0(options$output, ".log")
  } else {
    options$log_file
  }
  X <- read_numeric_matrix(options$expression_log, "gene_name")
  model <- readRDS(options$model)
  genes <- readLines(options$genes, warn = FALSE)
  if (length(genes) == 0L || anyNA(genes) || any(!nzchar(genes))) {
    stop("The gene shard must contain non-empty gene identifiers", call. = FALSE)
  }
  dimensions_message <- sprintf(
    paste0(
      "stage=extract input_dimensions=genes:%d samples:%d sources:%d ",
      "shard_genes:%d shard_id:%d scale=log2_cpm"
    ),
    nrow(X),
    ncol(X),
    ncol(model$W),
    length(genes),
    options$shard_id
  )
  paths_message <- sprintf(
    "stage=extract input_paths=%s output_paths=%s,%s",
    paste(
      c(options$expression_log, options$model, options$genes),
      collapse = ","
    ),
    options$output,
    extract_log_path
  )
  message(dimensions_message)
  message(paths_message)
  append_tensor_log(extract_log_path, dimensions_message)
  append_tensor_log(extract_log_path, paths_message)

  tensor <- extract_tensor_shard(
    X = X,
    model = model,
    genes = genes,
    num_cores = options$num_cores,
    log_file = extract_log_path
  )
  write_tensor_shard(options$output, tensor, options$shard_id)
  complete_message <- sprintf(
    paste0(
      "stage=extract event=stage_complete output_dimensions=genes:%d ",
      "samples:%d sources:%d output=%s scale=log2_cpm"
    ),
    length(genes),
    ncol(X),
    length(tensor),
    options$output
  )
  message(sprintf("%s utc_complete=%s", complete_message, tensor_utc_time()))
  append_tensor_log(extract_log_path, complete_message)
}

tryCatch(
  run_extract_tca_shard(),
  error = function(error) {
    error_message <- sprintf(
      "stage=extract status=failed utc_time=%s message=%s",
      tensor_utc_time(),
      conditionMessage(error)
    )
    message(error_message)
    if (!is.null(extract_log_path)) {
      append_tensor_log(extract_log_path, error_message)
    }
    quit(status = 1L)
  }
)
