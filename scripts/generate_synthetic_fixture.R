#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1L || !nzchar(arguments[[1L]])) {
  stop(
    "Usage: generate_synthetic_fixture.R OUTPUT_DIRECTORY",
    call. = FALSE
  )
}

output_directory <- arguments[[1L]]
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

seed <- 20260901L
set.seed(seed)

cell_types <- c(
  "B cells naive", "B cells memory", "Plasma cells", "T cells CD8",
  "T cells CD4 naive", "T cells CD4 memory resting",
  "T cells CD4 memory activated", "T cells follicular helper",
  "T cells regulatory (Tregs)", "T cells gamma delta",
  "NK cells resting", "NK cells activated", "Monocytes",
  "Macrophages M0", "Macrophages M1", "Macrophages M2",
  "Dendritic cells resting", "Dendritic cells activated",
  "Mast cells resting", "Mast cells activated", "Eosinophils", "Neutrophils"
)
cell_groups <- c(
  "B cells", "CD4 T cells", "CD8 T cells", "Gamma-delta T cells",
  "NK cells", "Monocyte/myeloid", "Neutrophils", "Eosinophils",
  "Dendritic cells", "Mast cells"
)
sample_ids <- sprintf("SAMPLE_%02d", seq_len(12L))
signature_gene_names <- sprintf("SYN_GENE_%03d", seq_len(66L))
signature_gene_ids <- sprintf("ENSGSYN%06d", seq_len(66L))
extra_gene_names <- sprintf("EXTRA_GENE_%02d", seq_len(6L))
extra_gene_ids <- sprintf("ENSGEXTRA%06d", seq_len(6L))

signature <- matrix(
  stats::runif(66L * 22L, min = 0.25, max = 1.25),
  nrow = 66L,
  ncol = 22L,
  dimnames = list(signature_gene_names, cell_types)
)
for (cell_type_index in seq_along(cell_types)) {
  marker_indices <- (3L * cell_type_index - 2L):(3L * cell_type_index)
  signature[marker_indices, cell_type_index] <-
    signature[marker_indices, cell_type_index] + 60
}

weights <- matrix(
  stats::rgamma(length(sample_ids) * length(cell_types), shape = 2, rate = 1),
  nrow = length(sample_ids),
  ncol = length(cell_types),
  dimnames = list(sample_ids, cell_types)
)
weights <- sweep(weights, 1L, rowSums(weights), "/")

extra_gene_expression <- matrix(
  stats::runif(length(extra_gene_names) * length(sample_ids), min = 0.5, max = 4),
  nrow = length(extra_gene_names),
  ncol = length(sample_ids),
  dimnames = list(extra_gene_names, sample_ids)
)
bulk_linear <- rbind(signature %*% t(weights), extra_gene_expression)
bulk_cpm <- sweep(bulk_linear, 2L, colSums(bulk_linear), "/") * 1e6
stopifnot(
  all(bulk_cpm > 0),
  max(abs(colSums(bulk_cpm) - 1e6)) < 1e-6
)

write_matrix <- function(matrix_value, path, id_column) {
  matrix_value |>
    as.data.frame(check.names = FALSE) |>
    tibble::rownames_to_column(var = id_column) |>
    readr::write_tsv(path, na = "")
}

expression_cpm <- bulk_cpm
rownames(expression_cpm) <- c(signature_gene_ids, extra_gene_ids)
write_matrix(
  expression_cpm,
  file.path(output_directory, "synthetic_cpm.tsv"),
  "gene_id"
)
write_matrix(
  signature,
  file.path(output_directory, "synthetic_signature.tsv"),
  "gene_symbol"
)
write_matrix(
  weights,
  file.path(output_directory, "precomputed_proportions.tsv"),
  "sample_id"
)

batch_indicator <- tibble::tibble(
  sample_id = sample_ids,
  batch_indicator = rep(c(0, 1), each = 6L)
)
readr::write_tsv(
  batch_indicator,
  file.path(output_directory, "batch_indicator.tsv"),
  na = ""
)

all_gene_ids <- c(signature_gene_ids, extra_gene_ids)
all_gene_names <- c(signature_gene_names, extra_gene_names)
gene_types <- ifelse(
  seq_along(all_gene_ids) %% 5L == 0L,
  "lncRNA",
  "protein_coding"
)
gtf_lines <- purrr::pmap_chr(
  list(
    gene_id = all_gene_ids,
    gene_name = all_gene_names,
    gene_type = gene_types,
    index = seq_along(all_gene_ids)
  ),
  function(gene_id, gene_name, gene_type, index) {
    start <- index * 1000L
    end <- start + 499L
    attributes <- sprintf(
      'gene_id "%s"; gene_name "%s"; gene_type "%s";',
      gene_id,
      gene_name,
      gene_type
    )
    paste(
      "chr1", "synthetic", "gene", start, end, ".", "+", ".",
      attributes,
      sep = "\t"
    )
  }
)
writeLines(gtf_lines, file.path(output_directory, "synthetic.gtf"))

writeLines(sample_ids, file.path(output_directory, "expected_samples.txt"))
writeLines(all_gene_ids, file.path(output_directory, "expected_genes.txt"))
writeLines(cell_groups, file.path(output_directory, "expected_groups.txt"))

common_inputs <- list(
  "cell_type_deconvolution.expression" = "tests/fixtures/synthetic_cpm.tsv",
  "cell_type_deconvolution.gtf" = "tests/fixtures/synthetic.gtf",
  "cell_type_deconvolution.covariates" = "tests/fixtures/batch_indicator.tsv",
  "cell_type_deconvolution.docker_image" = "celltype-deconvolution:test",
  "cell_type_deconvolution.pipeline_version" = "smoke-20260901",
  "cell_type_deconvolution.min_lm22_overlap" = 0.80,
  "cell_type_deconvolution.dtangle_marker_fraction" = 0.10,
  "cell_type_deconvolution.dtangle_quantile_normalize" = FALSE,
  "cell_type_deconvolution.group_mean_threshold" = 0.0001,
  "cell_type_deconvolution.zero_floor" = 0.000001,
  "cell_type_deconvolution.tca_max_iters" = 10L,
  "cell_type_deconvolution.random_seed" = seed,
  "cell_type_deconvolution.prepare_cpu" = 4L,
  "cell_type_deconvolution.prepare_memory" = "64 GB",
  "cell_type_deconvolution.prepare_disk_gb" = 400L,
  "cell_type_deconvolution.prepare_preemptible_attempts" = 2L,
  "cell_type_deconvolution.prepare_max_retries" = 2L,
  "cell_type_deconvolution.dtangle_cpu" = 4L,
  "cell_type_deconvolution.dtangle_memory" = "32 GB",
  "cell_type_deconvolution.dtangle_disk_gb" = 100L,
  "cell_type_deconvolution.dtangle_preemptible_attempts" = 2L,
  "cell_type_deconvolution.dtangle_max_retries" = 2L,
  "cell_type_deconvolution.proportions_cpu" = 2L,
  "cell_type_deconvolution.proportions_memory" = "16 GB",
  "cell_type_deconvolution.proportions_disk_gb" = 50L,
  "cell_type_deconvolution.proportions_preemptible_attempts" = 2L,
  "cell_type_deconvolution.proportions_max_retries" = 2L,
  "cell_type_deconvolution.fit_cpu" = 16L,
  "cell_type_deconvolution.fit_memory" = "192 GB",
  "cell_type_deconvolution.fit_disk_gb" = 750L,
  "cell_type_deconvolution.fit_preemptible_attempts" = 0L,
  "cell_type_deconvolution.fit_max_retries" = 1L,
  "cell_type_deconvolution.export_cpu" = 8L,
  "cell_type_deconvolution.export_memory" = "128 GB",
  "cell_type_deconvolution.export_disk_gb" = 500L,
  "cell_type_deconvolution.export_preemptible_attempts" = 0L,
  "cell_type_deconvolution.export_max_retries" = 1L,
  "cell_type_deconvolution.manifest_cpu" = 4L,
  "cell_type_deconvolution.manifest_memory" = "32 GB",
  "cell_type_deconvolution.manifest_disk_gb" = 100L,
  "cell_type_deconvolution.manifest_preemptible_attempts" = 1L,
  "cell_type_deconvolution.manifest_max_retries" = 2L
)
dtangle_inputs <- append(
  common_inputs,
  list(
    "cell_type_deconvolution.lm22" =
      "tests/fixtures/synthetic_signature.tsv"
  ),
  after = 2L
)
restart_inputs <- append(
  common_inputs,
  list(
    "cell_type_deconvolution.precomputed_proportions" =
      "tests/fixtures/precomputed_proportions.tsv"
  ),
  after = 2L
)
jsonlite::write_json(
  dtangle_inputs,
  file.path(output_directory, "dtangle.inputs.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)
jsonlite::write_json(
  restart_inputs,
  file.path(output_directory, "restart.inputs.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)
