source(testthat::test_path("../..", "R", "expression.R"), local = TRUE)

wdl_text <- function(path) {
  paste(
    readLines(testthat::test_path("../..", path), warn = FALSE),
    collapse = "\n"
  )
}

testthat::test_that("smoke fixture is BED and contains a duplicate gene symbol", {
  bed <- readr::read_tsv(
    testthat::test_path("..", "fixtures", "synthetic_expression.bed"),
    show_col_types = FALSE,
    progress = FALSE
  )
  testthat::expect_identical(
    names(bed)[1:4],
    c("#chr", "start", "end", "gene_id")
  )
  gtf <- read_gtf_gene_annotation(
    testthat::test_path("..", "fixtures", "synthetic.gtf")
  )
  mapped <- gtf$gene_name[match(bed$gene_id, gtf$gene_id)]
  testthat::expect_true(anyDuplicated(mapped) > 0L)
})

testthat::test_that("smoke assertions require unique compressed BED paths", {
  smoke_assertions <- paste(
    readLines(testthat::test_path("..", "smoke", "assert_outputs.R")),
    collapse = "\n"
  )
  purrr::walk(c(
    "cell_type_bed_paths <- normalizePath(cell_type_beds)",
    "inventory_basenames <- inventory$path",
    "all(grepl(\"[.]bed[.]gz$\", cell_type_bed_paths))",
    "anyDuplicated(cell_type_bed_paths) == 0L",
    "all(grepl(\"^[^/]+[.]bed[.]gz$\", inventory_basenames))",
    "anyDuplicated(inventory_basenames) == 0L",
    "setequal(basename(cell_type_bed_paths), inventory_basenames)",
    "public_inventory <- readr::read_tsv(",
    "stopifnot(identical(public_inventory, inventory))"
  ), ~ testthat::expect_match(smoke_assertions, .x, fixed = TRUE))
})

expect_task_contract <- function(path, task, inputs, outputs, defaults) {
  text <- wdl_text(path)
  testthat::expect_match(text, paste0("task ", task))
  purrr::walk(inputs, ~ testthat::expect_match(text, .x, fixed = TRUE))
  purrr::walk(outputs, ~ testthat::expect_match(text, .x, fixed = TRUE))
  purrr::walk(defaults, ~ testthat::expect_match(text, .x, fixed = TRUE))
  testthat::expect_match(text, "set -euo pipefail", fixed = TRUE)
  testthat::expect_match(text, "stage=", fixed = TRUE)
  testthat::expect_match(text, "start_time=", fixed = TRUE)
  testthat::expect_match(text, "completion_time=", fixed = TRUE)
  testthat::expect_match(text, "dimensions=", fixed = TRUE)
  testthat::expect_match(text, "outputs=", fixed = TRUE)
  testthat::expect_match(text, "error_status=", fixed = TRUE)
  testthat::expect_match(text, "preemptible_attempts", fixed = TRUE)
  testthat::expect_match(text, "max_retries", fixed = TRUE)
}

testthat::test_that("the expression WDL requires CPM and GTF", {
  text <- wdl_text("workflows/tasks/expression.wdl")
  testthat::expect_match(text, "File expression", fixed = TRUE)
  testthat::expect_match(text, "File gtf", fixed = TRUE)
  testthat::expect_match(text, "prepared_tca_cpm", fixed = TRUE)
  testthat::expect_match(text, "prepared_tca_log2_cpm", fixed = TRUE)
  testthat::expect_false(grepl("expression_type|gene_length", text))
})

testthat::test_that("top workflow uses one direct BED export", {
  text <- wdl_text("workflows/cell_type_deconvolution.wdl")
  testthat::expect_match(text, "call tca_tasks.ExportTcaBeds", fixed = TRUE)
  testthat::expect_identical(
    stringr::str_count(text, "call tca_tasks.ExportTcaBeds"),
    1L
  )
  purrr::walk(c(
    "Int export_cpu = 8", "String export_memory = \"128 GB\"",
    "Int export_disk_gb = 500", "Int export_preemptible_attempts = 0",
    "Int export_max_retries = 1"
  ), ~ testthat::expect_match(text, .x, fixed = TRUE))
  testthat::expect_match(text, "Array[File] cell_type_beds", fixed = TRUE)
  testthat::expect_match(
    text,
    "prepared_log2_cpm = PrepareExpression.prepared_tca_log2_cpm",
    fixed = TRUE
  )
  testthat::expect_match(
    text,
    "prepared_log2_cpm = PrepareExpression.prepared_dtangle_log2_cpm",
    fixed = TRUE
  )
  testthat::expect_false(grepl(
    "scatter|shard|hdf5|group_tsv|write_tsv",
    text,
    ignore.case = TRUE
  ))
  testthat::expect_false(grepl(
    "prepared_log2_cpm = RunDtangle.shared_bulk",
    text,
    fixed = TRUE
  ))
})

testthat::test_that("top workflow validates exactly one proportion mode before selection", {
  text <- wdl_text("workflows/cell_type_deconvolution.wdl")
  proportions_text <- wdl_text("workflows/tasks/proportions.wdl")

  testthat::expect_match(
    proportions_text,
    "task ValidateProportionMode",
    fixed = TRUE
  )
  purrr::walk(c(
    "File? lm22", "File? precomputed_proportions",
    "String selected_mode", "Boolean estimate_proportions"
  ), ~ testthat::expect_match(proportions_text, .x, fixed = TRUE))
  testthat::expect_match(
    text,
    "call proportion_tasks.ValidateProportionMode",
    fixed = TRUE
  )
  testthat::expect_match(
    text,
    "if (ValidateProportionMode.estimate_proportions)",
    fixed = TRUE
  )
  testthat::expect_match(
    text,
    paste0(
      "File proportions_for_processing = if ",
      "ValidateProportionMode.estimate_proportions"
    ),
    fixed = TRUE
  )
  testthat::expect_false(grepl(
    "Boolean estimate_proportions = !defined(precomputed_proportions)",
    text,
    fixed = TRUE
  ))
})

testthat::test_that("top-level manifest records effective parameters with a pinned TCA version", {
  text <- wdl_text("workflows/cell_type_deconvolution.wdl")
  workflow_inputs <- stringr::str_match(
    text,
    "workflow cell_type_deconvolution \\{\\n  input \\{([\\s\\S]*?)\\n  \\}"
  )[, 2L]
  testthat::expect_false(grepl("parameters_json", workflow_inputs, fixed = TRUE))
  testthat::expect_false(grepl("tca_version", workflow_inputs, fixed = TRUE))
  testthat::expect_match(
    text,
    "String tca_version = \"1.2.1\"",
    fixed = TRUE
  )
  testthat::expect_match(
    text,
    "EffectiveParameters effective_parameters = EffectiveParameters {",
    fixed = TRUE
  )
  testthat::expect_match(
    text,
    "File effective_parameters_json = write_json(effective_parameters)",
    fixed = TRUE
  )
  purrr::walk(c(
    "proportion_mode", "min_lm22_overlap", "dtangle_marker_fraction",
    "dtangle_marker_method", "dtangle_quantile_normalize",
    "group_mean_threshold", "zero_floor", "tca_max_iters",
    "random_seed", "scale"
  ), ~ testthat::expect_match(text, .x, fixed = TRUE))
  purrr::walk(c(
    "tca_shard_size", "write_tsv", "hdf5_gene_chunk_max",
    "hdf5_sample_chunk_max", "hdf5_gzip_level"
  ), ~ testthat::expect_false(grepl(.x, text, fixed = TRUE)))
  testthat::expect_match(
    text,
    "parameters_json = effective_parameters_json",
    fixed = TRUE
  )
  testthat::expect_match(text, "tca_version = tca_version", fixed = TRUE)
})

testthat::test_that("synthetic smoke fixtures are deterministic and restart without LM22", {
  generator <- testthat::test_path(
    "../..", "scripts", "generate_synthetic_fixture.R"
  )
  first <- tempfile("fixture-first-")
  second <- tempfile("fixture-second-")
  dir.create(first)
  dir.create(second)

  first_status <- system2("Rscript", c(generator, first))
  second_status <- system2("Rscript", c(generator, second))
  testthat::expect_identical(first_status, 0L)
  testthat::expect_identical(second_status, 0L)

  expected_files <- c(
    "synthetic_expression.bed", "synthetic.gtf", "synthetic_signature.tsv",
    "batch_indicator.tsv", "precomputed_proportions.tsv",
    "expected_samples.txt", "expected_gene_ids.txt", "expected_coordinates.tsv",
    "expected_groups.txt",
    "dtangle.inputs.json", "restart.inputs.json"
  )
  first_hashes <- tools::md5sum(file.path(first, expected_files))
  second_hashes <- tools::md5sum(file.path(second, expected_files))
  testthat::expect_false(anyNA(first_hashes))
  testthat::expect_identical(unname(first_hashes), unname(second_hashes))

  bed <- readr::read_tsv(
    file.path(first, "synthetic_expression.bed"),
    show_col_types = FALSE,
    progress = FALSE
  )
  signature <- readr::read_tsv(
    file.path(first, "synthetic_signature.tsv"),
    show_col_types = FALSE
  )
  values <- as.matrix(bed[-(1:4)])
  testthat::expect_identical(names(bed)[1:4], c("#chr", "start", "end", "gene_id"))
  testthat::expect_identical(dim(bed), c(73L, 16L))
  testthat::expect_true(all(values > 0))
  testthat::expect_lt(
    max(abs(colSums(values) - 1e6)),
    1e-6
  )
  testthat::expect_identical(dim(signature), c(66L, 23L))
  expected_lm22_types <- c(
    "B cells naive", "B cells memory", "Plasma cells", "T cells CD8",
    "T cells CD4 naive", "T cells CD4 memory resting",
    "T cells CD4 memory activated", "T cells follicular helper",
    "T cells regulatory (Tregs)", "T cells gamma delta",
    "NK cells resting", "NK cells activated", "Monocytes",
    "Macrophages M0", "Macrophages M1", "Macrophages M2",
    "Dendritic cells resting", "Dendritic cells activated",
    "Mast cells resting", "Mast cells activated", "Eosinophils", "Neutrophils"
  )
  testthat::expect_identical(names(signature)[-1L], expected_lm22_types)
  signature_values <- as.matrix(signature[-1L])
  designated_markers <- purrr::map2(
    seq_along(expected_lm22_types),
    expected_lm22_types,
    function(cell_type_index, cell_type) {
      observed <- signature$gene_symbol[
        signature_values[, cell_type_index] > 60
      ]
      expected <- sprintf(
        "SYN_GENE_%03d",
        (3L * cell_type_index - 2L):(3L * cell_type_index)
      )
      testthat::expect_identical(observed, expected, info = cell_type)
      observed
    }
  )
  testthat::expect_true(all(lengths(designated_markers) == 3L))
  expected_gene_ids <- readLines(
    file.path(first, "expected_gene_ids.txt"),
    warn = FALSE
  )
  expected_coordinates <- readr::read_tsv(
    file.path(first, "expected_coordinates.tsv"),
    show_col_types = FALSE,
    progress = FALSE
  ) |>
    tibble::as_tibble()
  testthat::expect_identical(bed$gene_id, expected_gene_ids)
  testthat::expect_identical(bed[1:4], expected_coordinates)

  gtf <- readLines(file.path(first, "synthetic.gtf"), warn = FALSE)
  testthat::expect_true(any(grepl('gene_type "protein_coding"', gtf)))
  testthat::expect_true(any(grepl('gene_type "lncRNA"', gtf)))
  restart <- jsonlite::read_json(file.path(first, "restart.inputs.json"))
  dtangle <- jsonlite::read_json(file.path(first, "dtangle.inputs.json"))
  testthat::expect_false("cell_type_deconvolution.lm22" %in% names(restart))
  testthat::expect_false(
    "cell_type_deconvolution.parameters_json" %in% names(restart)
  )
  testthat::expect_false(
    "cell_type_deconvolution.tca_version" %in% names(restart)
  )
  testthat::expect_false(
    "cell_type_deconvolution.parameters_json" %in% names(dtangle)
  )
  testthat::expect_false(
    "cell_type_deconvolution.tca_version" %in% names(dtangle)
  )
  testthat::expect_identical(
    restart$cell_type_deconvolution.precomputed_proportions,
    "tests/fixtures/precomputed_proportions.tsv"
  )
  expected_expression <- "tests/fixtures/synthetic_expression.bed"
  testthat::expect_identical(
    dtangle$cell_type_deconvolution.expression,
    expected_expression
  )
  testthat::expect_identical(
    restart$cell_type_deconvolution.expression,
    expected_expression
  )
  purrr::walk(c("export_cpu", "export_memory", "export_disk_gb",
    "export_preemptible_attempts", "export_max_retries"), function(name) {
    input_name <- paste0("cell_type_deconvolution.", name)
    testthat::expect_true(input_name %in% names(dtangle))
    testthat::expect_true(input_name %in% names(restart))
  })
  purrr::walk(c("shard", "hdf5", "tsv", "extract", "assemble"), function(term) {
    testthat::expect_false(any(grepl(
      term,
      names(dtangle),
      ignore.case = TRUE
    )))
    testthat::expect_false(any(grepl(
      term,
      names(restart),
      ignore.case = TRUE
    )))
  })
})

testthat::test_that("all modular WDL tasks declare CLI-aligned interfaces", {
  expect_task_contract(
    "workflows/tasks/expression.wdl", "PrepareExpression",
    c("File expression", "File gtf", "String docker_image"),
    c("File prepared_tca_cpm", "File prepared_tca_log2_cpm",
      "File prepared_dtangle_cpm", "File prepared_dtangle_log2_cpm",
      "File prepared_coordinates", "File mapping_report",
      "File excluded_genes", "File log"),
    c("Int cpu = 4", "String memory = \"64 GB\"", "Int disk_gb = 400",
      "Int preemptible_attempts = 2", "Int max_retries = 2")
  )
  expect_task_contract(
    "workflows/tasks/dtangle.wdl", "RunDtangle",
    c("File prepared_log2_cpm", "File lm22", "Float min_overlap",
      "Float marker_fraction", "String marker_method", "Boolean quantile_normalize"),
    c("File proportions", "File markers", "File metadata", "File overlap_report",
      "File transformed_lm22", "File shared_bulk", "File log"),
    c("Int cpu = 4", "String memory = \"32 GB\"", "Int disk_gb = 100",
      "Int preemptible_attempts = 2", "Int max_retries = 2")
  )
  expect_task_contract(
    "workflows/tasks/proportions.wdl", "ProcessProportions",
    c("File proportions", "Float mean_threshold", "Float zero_floor"),
    c("File original", "File combined", "File tca_weights", "File filter_report",
      "File log"),
    c("Int cpu = 2", "String memory = \"16 GB\"", "Int disk_gb = 50",
      "Int preemptible_attempts = 2", "Int max_retries = 2")
  )
  proportions_text <- wdl_text("workflows/tasks/proportions.wdl")
  purrr::walk(c(
    "task ValidateProportionMode", "File? lm22",
    "File? precomputed_proportions", "String selected_mode",
    "Boolean estimate_proportions", "File log"
  ), ~ testthat::expect_match(proportions_text, .x, fixed = TRUE))
  tca_text <- wdl_text("workflows/tasks/tca.wdl")
  testthat::expect_match(tca_text, "task FitTca", fixed = TRUE)
  purrr::walk(c(
    "File prepared_log2_cpm", "File tca_weights", "File? covariates",
    "Int max_iters", "Int random_seed", "File model",
    "File model_log", "File tca_expression", "File excluded_genes",
    "Int cpu = 16", "String memory = \"192 GB\"", "Int disk_gb = 750",
    "Int preemptible_attempts = 0", "Int max_retries = 1"
  ), ~ testthat::expect_match(tca_text, .x, fixed = TRUE))
  testthat::expect_match(tca_text, "task ExportTcaBeds", fixed = TRUE)
  purrr::walk(c(
    "File tca_expression", "File coordinates", "File model",
    "File tca_weights", "Array[File] cell_type_beds",
    "File cell_type_bed_inventory", "File reconstruction_by_sample",
    "File qc_summary", "File qc_plots", "File export_detail_log"
  ), ~ testthat::expect_match(tca_text, .x, fixed = TRUE))
  qc_text <- wdl_text("workflows/tasks/qc.wdl")
  purrr::walk(c(
    "task BuildManifest", "Array[File] cell_type_beds",
    "File cell_type_bed_inventory", "File export_qc_summary",
    "File export_qc_plots", "File output_manifest",
    "File qc_summary", "File qc_plots", "File provenance",
    "Int cpu = 4", "String memory = \"32 GB\"", "Int disk_gb = 100",
    "Int preemptible_attempts = 1", "Int max_retries = 2"
  ), ~ testthat::expect_match(qc_text, .x, fixed = TRUE))
})

testthat::test_that("optional WDL file arguments use quoted shell arrays", {
  tca_text <- wdl_text("workflows/tasks/tca.wdl")
  qc_text <- wdl_text("workflows/tasks/qc.wdl")

  testthat::expect_identical(
    stringr::str_count(
      tca_text,
      stringr::fixed('covariates_arguments=()')
    ),
    2L
  )
  testthat::expect_identical(
    stringr::str_count(
      tca_text,
      stringr::fixed(
        'covariates_arguments=(--covariates "$covariates_path")'
      )
    ),
    2L
  )
  testthat::expect_identical(
    stringr::str_count(
      tca_text,
      stringr::fixed('"${covariates_arguments[@]}"')
    ),
    2L
  )
  testthat::expect_match(qc_text, "parameters_arguments=()", fixed = TRUE)
  testthat::expect_match(
    qc_text,
    'parameters_arguments=(--parameters-json "$parameters_path")',
    fixed = TRUE
  )
  testthat::expect_match(
    qc_text,
    '"${parameters_arguments[@]}"',
    fixed = TRUE
  )
  testthat::expect_match(
    qc_text,
    'dtangle_metadata_arguments=()',
    fixed = TRUE
  )
  testthat::expect_match(
    qc_text,
    'dtangle_metadata_arguments=(--dtangle-metadata "$dtangle_metadata_path")',
    fixed = TRUE
  )
  testthat::expect_match(
    qc_text,
    '"${dtangle_metadata_arguments[@]}"',
    fixed = TRUE
  )
  testthat::expect_false(grepl(
    'then "--covariates " + select_first',
    tca_text,
    fixed = TRUE
  ))
  testthat::expect_false(grepl(
    'then "--parameters-json " + select_first',
    qc_text,
    fixed = TRUE
  ))
})

testthat::test_that("FitTca logs decompressed gene and sample dimensions", {
  text <- wdl_text("workflows/tasks/tca.wdl")

  testthat::expect_match(
    text,
    "gzip -cd outputs/tca_expression.tsv.gz",
    fixed = TRUE
  )
  testthat::expect_match(text, "genes:%s,samples:%s", fixed = TRUE)
  testthat::expect_false(grepl(
    "wc -l < outputs/tca_expression.tsv.gz",
    text,
    fixed = TRUE
  ))
})

testthat::test_that("the logging checker rejects incomplete command logging", {
  fixture <- tempfile(fileext = ".wdl")
  checker_path <- testthat::test_path(
    "../..", "tools", "check_wdl_logging.py"
  )
  writeLines(c(
    "version 1.1", "task MissingLogging {", "  command <<<", "    echo hello", "  >>>", "}"
  ), fixture)
  status <- system2(
    "python3", c(checker_path, fixture), stdout = FALSE,
    stderr = FALSE
  )
  testthat::expect_false(identical(status, 0L))
})

testthat::test_that("the logging checker checks brace-form commands", {
  incomplete_fixture <- tempfile(fileext = ".wdl")
  complete_fixture <- tempfile(fileext = ".wdl")
  checker_path <- testthat::test_path(
    "../..", "tools", "check_wdl_logging.py"
  )
  writeLines(c(
    "version 1.1", "task MissingBraceLogging {", "  command {",
    "    echo 'stage=brace start_time=now'", "  }", "}"
  ), incomplete_fixture)
  writeLines(c(
    "version 1.1", "task CompleteBraceLogging {", "  command {",
    "    echo 'stage=brace start_time=now completion_time=now dimensions=1 outputs=result'",
    "  }", "}"
  ), complete_fixture)

  incomplete_status <- system2(
    "python3", c(checker_path, incomplete_fixture),
    stdout = FALSE, stderr = FALSE
  )
  complete_status <- system2(
    "python3", c(checker_path, complete_fixture),
    stdout = FALSE, stderr = FALSE
  )
  miniwdl_status <- system2(
    "miniwdl", c("check", complete_fixture), stdout = FALSE, stderr = FALSE
  )

  testthat::expect_false(identical(incomplete_status, 0L))
  testthat::expect_identical(complete_status, 0L)
  testthat::expect_identical(miniwdl_status, 0L)
})

testthat::test_that("pipeline CI keeps the runner group for smoke assertions", {
  workflow_path <- testthat::test_path(
    "../..", ".github", "workflows", "pipeline-ci.yml"
  )
  workflow <- paste(readLines(workflow_path, warn = FALSE), collapse = "\n")

  assertions <- stringr::str_extract_all(
    workflow,
    "(?s)- name: Assert .*? smoke outputs\\n.*?(?=\\n      - name:|\\z)"
  )[[1L]]

  testthat::expect_length(assertions, 2L)
  purrr::walk(
    assertions,
    ~ testthat::expect_match(.x, "--group-add \\\"\\$\\(id -g\\)\\\"", fixed = FALSE)
  )
})
