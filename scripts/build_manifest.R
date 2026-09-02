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
      help = "TSV output inventory from direct BED export."
    ),
    optparse::make_option(
      "--export-qc-summary",
      dest = "export_qc_summary",
      type = "character",
      help = "Direct BED export QC summary TSV."
    ),
    optparse::make_option(
      "--original-proportions",
      dest = "original_proportions",
      type = "character",
      help = "Sample-by-LM22 proportion TSV."
    ),
    optparse::make_option(
      "--combined-proportions",
      dest = "combined_proportions",
      type = "character",
      help = "Sample-by-combined-group proportion TSV."
    ),
    optparse::make_option(
      "--tca-weights",
      dest = "tca_weights",
      type = "character",
      help = "Sample-by-retained-group TCA weight TSV."
    ),
    optparse::make_option(
      "--filter-report",
      dest = "filter_report",
      type = "character",
      help = "Cell-group filter and zero-adjustment report TSV."
    ),
    optparse::make_option(
      "--model",
      type = "character",
      help = "Fitted TCA model RDS."
    ),
    optparse::make_option(
      "--model-log",
      dest = "model_log",
      type = "character",
      help = "TCA model log with internal iteration records."
    ),
    optparse::make_option(
      "--dtangle-metadata",
      dest = "dtangle_metadata",
      type = "character",
      default = NULL,
      help = "Optional dtangle metadata JSON with LM22 QC."
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
      "--qc-output",
      dest = "qc_output",
      type = "character",
      help = "Final pipeline QC summary TSV path."
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
    "outputs", "export_qc_summary",
    "original_proportions", "combined_proportions", "tca_weights",
    "filter_report", "model", "model_log",
    "container_image", "output", "qc_output"
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
  purrr::walk(
    unique(c(dirname(options$output), dirname(options$qc_output))),
    ~ dir.create(.x, recursive = TRUE, showWarnings = FALSE)
  )
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
  export_qc_summary <- readr::read_tsv(
    options$export_qc_summary,
    show_col_types = FALSE,
    progress = FALSE
  )
  original_proportions <- read_numeric_matrix(
    options$original_proportions,
    "sample_id"
  )
  combined_proportions <- read_numeric_matrix(
    options$combined_proportions,
    "sample_id"
  )
  tca_weights <- read_numeric_matrix(options$tca_weights, "sample_id")
  filter_report <- readr::read_tsv(
    options$filter_report,
    col_types = readr::cols(
      cell_group = readr::col_character(),
      cohort_mean = readr::col_double(),
      threshold = readr::col_double(),
      retained = readr::col_logical(),
      filter_reason = readr::col_character(),
      zero_count_before = readr::col_integer(),
      zero_floor = readr::col_double()
    ),
    show_col_types = FALSE,
    progress = FALSE
  )
  tca_model <- readRDS(options$model)
  tca_log_lines <- readLines(options$model_log, warn = FALSE)
  dtangle_metadata <- if (is.null(options$dtangle_metadata) ||
      !nzchar(options$dtangle_metadata)) {
    NULL
  } else {
    jsonlite::read_json(options$dtangle_metadata, simplifyVector = FALSE)
  }
  parameters <- if (is.null(options$parameters_json) ||
      !nzchar(options$parameters_json)) {
    list(scale = "log2_cpm")
  } else {
    jsonlite::read_json(options$parameters_json, simplifyVector = FALSE)
  }
  dimensions_message <- sprintf(
    paste0(
      "stage=manifest input_dimensions=outputs:%d samples:%d ",
      "lm22_types:%d retained_groups:%d"
    ),
    nrow(outputs),
    nrow(original_proportions),
    ncol(original_proportions),
    ncol(tca_weights)
  )
  paths_message <- sprintf(
    "stage=manifest input_paths=%s output_paths=%s",
    options$outputs,
    paste(c(options$output, options$qc_output), collapse = ",")
  )
  message(dimensions_message)
  message(paths_message)
  append_tensor_log(manifest_log_path, dimensions_message)
  append_tensor_log(manifest_log_path, paths_message)

  manifest <- build_output_manifest(
    outputs = outputs,
    tca_version = options$tca_version,
    parameters = parameters,
    container_image = options$container_image
  )
  qc_summary <- build_pipeline_qc_summary(
    export_summary = export_qc_summary,
    original_proportions = original_proportions,
    combined_proportions = combined_proportions,
    tca_weights = tca_weights,
    filter_report = filter_report,
    tca_model = tca_model,
    tca_log_lines = tca_log_lines,
    dtangle_metadata = dtangle_metadata
  )
  write_output_manifest(options$output, manifest)
  readr::write_tsv(qc_summary, options$qc_output, na = "")
  complete_message <- sprintf(
    paste0(
      "stage=manifest event=stage_complete outputs=%d qc_metrics=%d ",
      "output_paths=%s"
    ),
    length(manifest$outputs),
    nrow(qc_summary),
    paste(c(options$output, options$qc_output), collapse = ",")
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
