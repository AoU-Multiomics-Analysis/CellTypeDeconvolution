#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]))
source(file.path(dirname(script_path), "bootstrap.R"))

utc_time <- function() format(Sys.time(), tz = "UTC", usetz = TRUE)

run_prepare_expression <- function() {
  option_list <- list(
    optparse::make_option(
      "--expression",
      type = "character",
      help = "BED or BED.GZ positive linear CPM matrix with coordinate columns."
    ),
    optparse::make_option(
      "--gtf",
      type = "character",
      help = "GTF or GTF.GZ file that maps gene_id values to gene_name values."
    ),
    optparse::make_option(
      "--output-dir",
      dest = "output_dir",
      type = "character",
      help = "Directory for prepared CPM outputs."
    )
  )
  options <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
  required_options <- c("expression", "gtf", "output_dir")
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

  message(sprintf("stage=prepare_expression utc_start=%s", utc_time()))
  expression <- read_expression_bed(options$expression)
  message(sprintf(
    "stage=prepare_expression input_rows=%d samples=%d",
    nrow(expression$cpm),
    ncol(expression$cpm)
  ))
  annotation <- read_gtf_gene_annotation(options$gtf)
  message(sprintf(
    "stage=prepare_expression gtf_gene_count=%d",
    nrow(annotation)
  ))
  result <- prepare_expression_bed(expression = expression, annotation = annotation)
  message(sprintf(
    "stage=prepare_expression mapped_tca_rows=%d dtangle_symbols=%d samples=%d",
    nrow(result$tca_cpm),
    nrow(result$dtangle_cpm),
    ncol(result$tca_cpm)
  ))

  dir.create(options$output_dir, recursive = TRUE, showWarnings = FALSE)
  output_paths <- list(
    tca_cpm = file.path(options$output_dir, "prepared_tca_cpm.tsv.gz"),
    tca_log = file.path(options$output_dir, "prepared_tca_log2_cpm.tsv.gz"),
    dtangle_cpm = file.path(options$output_dir, "prepared_dtangle_cpm.tsv.gz"),
    dtangle_log = file.path(options$output_dir, "prepared_dtangle_log2_cpm.tsv.gz"),
    coordinates = file.path(options$output_dir, "prepared_coordinates.tsv"),
    mapping_report = file.path(options$output_dir, "gene_mapping_report.tsv"),
    excluded_genes = file.path(options$output_dir, "excluded_genes.tsv")
  )
  purrr::iwalk(output_paths, function(path, name) {
    message(sprintf("stage=prepare_expression output_%s=%s", name, path))
  })
  write_numeric_matrix(result$tca_cpm, output_paths$tca_cpm, "gene_id")
  write_numeric_matrix(result$tca_log_expression, output_paths$tca_log, "gene_id")
  write_numeric_matrix(result$dtangle_cpm, output_paths$dtangle_cpm, "gene_name")
  write_numeric_matrix(result$dtangle_log_expression, output_paths$dtangle_log, "gene_name")
  readr::write_tsv(result$coordinates, output_paths$coordinates, na = "")
  readr::write_tsv(result$mapping_report, output_paths$mapping_report, na = "")
  readr::write_tsv(result$excluded_genes, output_paths$excluded_genes, na = "")
  message(sprintf("stage=prepare_expression utc_complete=%s", utc_time()))
}

tryCatch(
  run_prepare_expression(),
  error = function(error) {
    message(sprintf(
      "stage=prepare_expression status=failed utc_time=%s message=%s",
      utc_time(), conditionMessage(error)
    ))
    quit(status = 1L)
  }
)
