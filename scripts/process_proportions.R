#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "bootstrap.R"))

utc_time <- function() format(Sys.time(), tz = "UTC", usetz = TRUE)

run_proportion_stage <- function() {
  option_list <- list(
    optparse::make_option(
      "--proportions",
      type = "character",
      help = "Sample-by-22 LM22 proportion matrix TSV with sample_id first column."
    ),
    optparse::make_option(
      "--mean-threshold",
      dest = "mean_threshold",
      type = "double",
      default = pipeline_defaults()$group_mean_threshold,
      help = "Retain cell groups whose cohort mean is at least this value."
    ),
    optparse::make_option(
      "--zero-floor",
      dest = "zero_floor",
      type = "double",
      default = pipeline_defaults()$zero_floor,
      help = "Replace exact zero retained proportions with this positive value."
    ),
    optparse::make_option(
      "--output-dir",
      dest = "output_dir",
      type = "character",
      help = "Directory for proportion processing outputs."
    )
  )
  options <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
  required_options <- c("proportions", "output_dir")
  missing_options <- required_options[vapply(
    options[required_options],
    function(value) is.null(value) || !nzchar(value),
    logical(1)
  )]
  if (length(missing_options) > 0L) {
    stop(
      sprintf("Missing required options: %s", paste(missing_options, collapse = ", ")),
      call. = FALSE
    )
  }

  message(sprintf("stage=proportions utc_start=%s", utc_time()))
  proportions <- read_numeric_matrix(options$proportions, "sample_id")
  message(sprintf(
    "stage=proportions input_dimensions=samples:%d lm22_cell_types:%d",
    nrow(proportions), ncol(proportions)
  ))
  message(sprintf(
    "stage=proportions settings=mean_threshold:%g zero_floor:%g",
    options$mean_threshold, options$zero_floor
  ))
  result <- process_proportions(
    proportions = proportions,
    mean_threshold = options$mean_threshold,
    zero_floor = options$zero_floor
  )
  message(sprintf(
    "stage=proportions output_dimensions=samples:%d combined_groups:%d retained_groups:%d",
    nrow(result$original), ncol(result$combined), ncol(result$tca_weights)
  ))

  dir.create(options$output_dir, recursive = TRUE, showWarnings = FALSE)
  output_paths <- list(
    original = file.path(options$output_dir, "proportions_lm22.tsv"),
    combined = file.path(options$output_dir, "proportions_combined.tsv"),
    tca_weights = file.path(options$output_dir, "proportions_tca_weights.tsv"),
    report = file.path(options$output_dir, "cell_group_filter_report.tsv")
  )
  message(sprintf(
    "stage=proportions output_paths=%s",
    paste(unlist(output_paths, use.names = FALSE), collapse = ",")
  ))
  write_numeric_matrix(result$original, output_paths$original, "sample_id")
  write_numeric_matrix(result$combined, output_paths$combined, "sample_id")
  write_numeric_matrix(result$tca_weights, output_paths$tca_weights, "sample_id")
  readr::write_tsv(result$report, output_paths$report, na = "")
  message(sprintf("stage=proportions utc_complete=%s", utc_time()))
}

tryCatch(
  run_proportion_stage(),
  error = function(error) {
    message(sprintf(
      "stage=proportions status=failed utc_time=%s message=%s",
      utc_time(), conditionMessage(error)
    ))
    quit(status = 1L)
  }
)
