#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]))
source(file.path(dirname(script_path), "bootstrap.R"))

utc_time <- function() format(Sys.time(), tz = "UTC", usetz = TRUE)

run_dtangle_stage <- function() {
  option_list <- list(
    optparse::make_option(
      "--expression",
      type = "character",
      help = "Coordinate-preserving BED of positive linear CPM values."
    ),
    optparse::make_option(
      "--gtf",
      type = "character",
      help = "GTF used to map expression gene_id values to LM22 gene symbols."
    ),
    optparse::make_option(
      "--lm22",
      type = "character",
      help = "Standard positive linear LM22 matrix TSV with gene_symbol first column."
    ),
    optparse::make_option(
      "--min-overlap",
      dest = "min_overlap",
      type = "double",
      default = pipeline_defaults()$min_lm22_overlap,
      help = "Minimum fraction of LM22 genes shared with bulk expression."
    ),
    optparse::make_option(
      "--marker-fraction",
      dest = "marker_fraction",
      type = "double",
      default = pipeline_defaults()$marker_fraction,
      help = "Fraction of reference genes used as dtangle markers."
    ),
    optparse::make_option(
      "--marker-method",
      dest = "marker_method",
      type = "character",
      default = "ratio",
      help = "dtangle marker method; only ratio is supported."
    ),
    optparse::make_option(
      "--quantile-normalize",
      dest = "quantile_normalize",
      action = "store_true",
      default = FALSE,
      help = "Jointly quantile normalize LM22 and bulk profiles for dtangle only."
    ),
    optparse::make_option(
      "--output-dir",
      dest = "output_dir",
      type = "character",
      help = "Directory for dtangle outputs."
    )
  )
  options <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
  required_options <- c("expression", "gtf", "lm22", "output_dir")
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

  message(sprintf("stage=dtangle utc_start=%s", utc_time()))
  expression <- read_expression_bed(options$expression)
  annotation <- read_gtf_gene_annotation(options$gtf)
  dtangle_expression <- make_dtangle_expression(expression, annotation)
  bulk_log <- dtangle_expression$log_expression
  lm22_linear <- read_numeric_matrix(options$lm22, "gene_symbol")
  message(sprintf(
    "stage=dtangle bulk_dimensions=genes:%d samples:%d lm22_dimensions=genes:%d cell_types:%d",
    nrow(bulk_log), ncol(bulk_log), nrow(lm22_linear), ncol(lm22_linear)
  ))
  message(sprintf(
    paste0(
      "stage=dtangle settings=min_overlap:%.3f marker_fraction:%.3f ",
      "marker_method:%s quantile_normalize:%s"
    ),
    options$min_overlap,
    options$marker_fraction,
    options$marker_method,
    options$quantile_normalize
  ))
  inputs <- prepare_dtangle_inputs(
    bulk_log = bulk_log,
    lm22_linear = lm22_linear,
    min_overlap = options$min_overlap,
    quantile_normalize = options$quantile_normalize
  )
  message(sprintf(
    "stage=dtangle overlap=shared_genes:%d overlap_fraction:%.3f",
    inputs$overlap_count,
    inputs$overlap_fraction
  ))
  fit <- estimate_dtangle(
    inputs,
    marker_fraction = options$marker_fraction,
    marker_method = options$marker_method
  )
  marker_counts <- unlist(fit$metadata$marker_counts, use.names = TRUE)
  message(sprintf("stage=dtangle package=dtangle version=%s", fit$metadata$dtangle_version))
  message(sprintf(
    "stage=dtangle markers_per_cell_type=%s",
    paste(sprintf("%s:%d", names(marker_counts), marker_counts), collapse = ",")
  ))

  dir.create(options$output_dir, recursive = TRUE, showWarnings = FALSE)
  output_paths <- list(
    proportions = file.path(options$output_dir, "dtangle_proportions.tsv"),
    markers = file.path(options$output_dir, "dtangle_markers.tsv"),
    metadata = file.path(options$output_dir, "dtangle_metadata.json"),
    overlap = file.path(options$output_dir, "dtangle_overlap.tsv"),
    lm22_log = file.path(options$output_dir, "dtangle_lm22_log.tsv.gz"),
    shared_bulk = file.path(options$output_dir, "dtangle_shared_bulk.tsv.gz")
  )
  message(sprintf(
    "stage=dtangle output_paths=%s",
    paste(unlist(output_paths, use.names = FALSE), collapse = ",")
  ))
  write_numeric_matrix(fit$proportions, output_paths$proportions, "sample_id")
  readr::write_tsv(fit$markers, output_paths$markers, na = "")
  jsonlite::write_json(fit$metadata, output_paths$metadata, auto_unbox = TRUE, pretty = TRUE)
  readr::write_tsv(inputs$overlap_report, output_paths$overlap, na = "")
  write_numeric_matrix(inputs$transformed_lm22, output_paths$lm22_log, "gene_symbol")
  write_numeric_matrix(inputs$shared_bulk, output_paths$shared_bulk, "gene_symbol")
  message(sprintf("stage=dtangle utc_complete=%s", utc_time()))
}

tryCatch(
  run_dtangle_stage(),
  error = function(error) {
    message(sprintf(
      "stage=dtangle status=failed utc_time=%s message=%s",
      utc_time(), conditionMessage(error)
    ))
    quit(status = 1L)
  }
)
