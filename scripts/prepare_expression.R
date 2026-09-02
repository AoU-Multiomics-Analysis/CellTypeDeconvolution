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
      help = "Gene-by-sample positive linear CPM matrix TSV with a gene_id first column."
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
  cpm <- read_numeric_matrix(options$expression, "gene_id")
  message(sprintf(
    "stage=prepare_expression input_dimensions=genes:%d samples:%d",
    nrow(cpm),
    ncol(cpm)
  ))
  annotation <- read_gtf_gene_annotation(options$gtf)
  message(sprintf(
    "stage=prepare_expression gtf_gene_count=%d",
    nrow(annotation)
  ))
  result <- prepare_expression(cpm = cpm, annotation = annotation)
  mapped_gene_name_count <- sum(result$mapping_report$mapping_action %in% c(
    "mapped", "duplicate_gene_name_collapsed"
  ))
  message(sprintf(
    "stage=prepare_expression mapped_gene_name_count=%d collapsed_gene_name_count=%d",
    mapped_gene_name_count,
    nrow(result$cpm)
  ))

  dir.create(options$output_dir, recursive = TRUE, showWarnings = FALSE)
  output_paths <- list(
    cpm = file.path(options$output_dir, "prepared_cpm.tsv.gz"),
    log_expression = file.path(options$output_dir, "prepared_log2_cpm.tsv.gz"),
    mapping_report = file.path(options$output_dir, "gene_mapping_report.tsv"),
    excluded_genes = file.path(options$output_dir, "excluded_genes.tsv")
  )
  message(sprintf(
    "stage=prepare_expression output_paths=%s",
    paste(unlist(output_paths, use.names = FALSE), collapse = ",")
  ))
  write_numeric_matrix(result$cpm, output_paths$cpm, "gene_name")
  write_numeric_matrix(result$log_expression, output_paths$log_expression, "gene_name")
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
