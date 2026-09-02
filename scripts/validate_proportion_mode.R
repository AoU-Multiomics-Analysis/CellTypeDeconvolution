#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "bootstrap.R"))

parse_presence_flag <- function(value, label) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !value %in% c("true", "false")) {
    stop(sprintf("%s must be true or false", label), call. = FALSE)
  }
  identical(value, "true")
}

run_validation <- function() {
  option_list <- list(
    optparse::make_option(
      "--lm22-defined",
      dest = "lm22_defined",
      type = "character"
    ),
    optparse::make_option(
      "--precomputed-defined",
      dest = "precomputed_defined",
      type = "character"
    ),
    optparse::make_option(
      "--output-dir",
      dest = "output_dir",
      type = "character"
    )
  )
  options <- optparse::parse_args(optparse::OptionParser(
    option_list = option_list
  ))
  required <- c("lm22_defined", "precomputed_defined", "output_dir")
  missing <- required[vapply(
    options[required],
    function(value) is.null(value) || !nzchar(value),
    logical(1)
  )]
  if (length(missing) > 0L) {
    stop(
      sprintf("Missing required options: %s", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  message(sprintf("stage=validate_proportion_mode utc_start=%s", tensor_utc_time()))
  result <- validate_proportion_mode(
    parse_presence_flag(options$lm22_defined, "lm22_defined"),
    parse_presence_flag(
      options$precomputed_defined,
      "precomputed_proportions_defined"
    )
  )
  dir.create(options$output_dir, recursive = TRUE, showWarnings = FALSE)
  mode_path <- file.path(options$output_dir, "selected_mode.txt")
  estimate_path <- file.path(options$output_dir, "estimate_proportions.txt")
  writeLines(result$selected_mode, mode_path)
  writeLines(tolower(as.character(result$estimate_proportions)), estimate_path)
  message(sprintf(
    paste0(
      "stage=validate_proportion_mode selected_mode=%s ",
      "estimate_proportions=%s output_paths=%s,%s utc_complete=%s"
    ),
    result$selected_mode,
    tolower(as.character(result$estimate_proportions)),
    mode_path,
    estimate_path,
    tensor_utc_time()
  ))
}

tryCatch(
  run_validation(),
  error = function(error) {
    message(sprintf(
      "stage=validate_proportion_mode status=failed utc_time=%s message=%s",
      tensor_utc_time(),
      conditionMessage(error)
    ))
    quit(status = 1L)
  }
)
