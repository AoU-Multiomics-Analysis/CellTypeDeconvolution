#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "bootstrap.R"))

manifest_log_path <- NULL

run_build_manifest <- function() {
  message(sprintf("stage=manifest utc_start=%s", tensor_utc_time()))
  option_list <- list(
    optparse::make_option(
      "--outputs",
      type = "character",
      help = "TSV output inventory from assembly."
    ),
    optparse::make_option(
      "--pipeline-version",
      dest = "pipeline_version",
      type = "character",
      help = "Pipeline version."
    ),
    optparse::make_option(
      "--tca-version",
      dest = "tca_version",
      type = "character",
      default = "1.2.1",
      help = "TCA package version."
    ),
    optparse::make_option(
      "--parameters-json",
      dest = "parameters_json",
      type = "character",
      default = NULL,
      help = "Optional JSON object with scientific and storage parameters."
    ),
    optparse::make_option(
      "--container-image",
      dest = "container_image",
      type = "character",
      help = "Container image name and immutable digest."
    ),
    optparse::make_option(
      "--output",
      type = "character",
      help = "Output manifest JSON path."
    ),
    optparse::make_option(
      "--log-file",
      dest = "log_file",
      type = "character",
      default = NULL,
      help = "Manifest log path."
    )
  )
  options <- optparse::parse_args(optparse::OptionParser(
    option_list = option_list
  ))
  required_options <- c(
    "outputs", "pipeline_version", "container_image", "output"
  )
  missing_options <- required_options[vapply(
    options[required_options],
    function(value) is.null(value) || !nzchar(value),
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
  manifest_log_path <<- if (is.null(options$log_file)) {
    paste0(options$output, ".log")
  } else {
    options$log_file
  }
  outputs <- readr::read_tsv(
    options$outputs,
    col_types = readr::cols(
      logical_name = readr::col_character(),
      path = readr::col_character(),
      n_genes = readr::col_integer(),
      n_samples = readr::col_integer(),
      scale = readr::col_character(),
      cell_group = readr::col_character()
    ),
    show_col_types = FALSE,
    progress = FALSE
  )
  parameters <- if (is.null(options$parameters_json) ||
      !nzchar(options$parameters_json)) {
    list(
      tensor_shard_size = pipeline_defaults()$tensor_shard_size,
      hdf5_gene_chunk_max = 500L,
      hdf5_sample_chunk_max = 256L,
      hdf5_gzip_level = 6L,
      scale = "log2_cpm"
    )
  } else {
    jsonlite::read_json(options$parameters_json, simplifyVector = FALSE)
  }
  dimensions_message <- sprintf(
    "stage=manifest input_dimensions=outputs:%d",
    nrow(outputs)
  )
  paths_message <- sprintf(
    "stage=manifest input_paths=%s output_paths=%s",
    options$outputs,
    options$output
  )
  message(dimensions_message)
  message(paths_message)
  append_tensor_log(manifest_log_path, dimensions_message)
  append_tensor_log(manifest_log_path, paths_message)

  manifest <- build_output_manifest(
    outputs = outputs,
    pipeline_version = options$pipeline_version,
    tca_version = options$tca_version,
    parameters = parameters,
    container_image = options$container_image
  )
  write_output_manifest(options$output, manifest)
  complete_message <- sprintf(
    "stage=manifest event=stage_complete outputs=%d output_paths=%s",
    length(manifest$outputs),
    options$output
  )
  message(sprintf("%s utc_complete=%s", complete_message, tensor_utc_time()))
  append_tensor_log(manifest_log_path, complete_message)
}

tryCatch(
  run_build_manifest(),
  error = function(error) {
    error_message <- sprintf(
      "stage=manifest status=failed utc_time=%s message=%s",
      tensor_utc_time(),
      conditionMessage(error)
    )
    message(error_message)
    if (!is.null(manifest_log_path)) {
      append_tensor_log(manifest_log_path, error_message)
    }
    quit(status = 1L)
  }
)
