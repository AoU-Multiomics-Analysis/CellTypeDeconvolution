#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]))
source(file.path(dirname(script_path), "bootstrap.R"))

option_list <- list(
  optparse::make_option(
    "--expression",
    type = "character",
    help = "Gene-by-sample expression matrix TSV with a gene_id first column."
  ),
  optparse::make_option(
    "--expression-type",
    dest = "expression_type",
    type = "character",
    help = "Expression scale: counts or tpm."
  ),
  optparse::make_option(
    "--annotation",
    type = "character",
    default = "",
    help = "Optional gene annotation TSV."
  ),
  optparse::make_option(
    "--output-dir",
    dest = "output_dir",
    type = "character",
    help = "Directory for prepared expression outputs."
  )
)
options <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
required_options <- c("expression", "expression_type", "output_dir")
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

utc_time <- function() format(Sys.time(), tz = "UTC", usetz = TRUE)
message(sprintf("stage=prepare_expression utc_start=%s", utc_time()))
expression <- read_numeric_matrix(options$expression, "gene_id")
message(sprintf(
  "stage=prepare_expression input_dimensions=genes:%d samples:%d",
  nrow(expression),
  ncol(expression)
))
result <- prepare_expression(
  expression = expression,
  expression_type = options$expression_type,
  annotation = read_optional_annotation(options$annotation)
)

dir.create(options$output_dir, recursive = TRUE, showWarnings = FALSE)
output_paths <- list(
  tpm = file.path(options$output_dir, "prepared_tpm.tsv.gz"),
  log_expression = file.path(options$output_dir, "prepared_log2_tpm_plus_1.tsv.gz"),
  mapping_report = file.path(options$output_dir, "gene_mapping_report.tsv"),
  excluded_genes = file.path(options$output_dir, "excluded_genes.tsv")
)
message(sprintf(
  "stage=prepare_expression output_paths=%s",
  paste(unlist(output_paths, use.names = FALSE), collapse = ",")
))
write_numeric_matrix(result$tpm, output_paths$tpm, "gene_id")
write_numeric_matrix(result$log_expression, output_paths$log_expression, "gene_id")
readr::write_tsv(result$mapping_report, output_paths$mapping_report, na = "")
readr::write_tsv(result$excluded_genes, output_paths$excluded_genes, na = "")
message(sprintf("stage=prepare_expression utc_complete=%s", utc_time()))
