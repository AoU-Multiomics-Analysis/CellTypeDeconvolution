# Direct Cell-Type BED Outputs Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace matrix-only input, HDF5 assembly, and tensor sharding with coordinate-preserving BED input and one direct BED.GZ expression output per retained major cell group.

**Architecture:** Expression preparation creates two views from one BED input. The gene-ID view preserves coordinates for full-matrix TCA, while the gene-symbol view aggregates duplicate symbols only for dtangle. One export task extracts the complete TCA tensor, validates it, calculates reconstruction QC, and writes every retained group as BED.GZ.

**Tech Stack:** R 4.5.3, tidyverse 2.0.0, dtangle 2.0.10, TCA 1.2.1, WDL 1.1, miniwdl 1.15.0 in GitHub Actions, Micromamba, Terra, GitHub Actions.

**Spec:** `docs/superpowers/specs/2026-09-01-bed-output-revision-design.md`

## Global Constraints

- The input expression file is BED or BED.GZ with leading columns `#chr`, `start`, `end`, and `gene_id` in that exact order.
- Sample values are normalized linear CPM. They must be finite and strictly positive.
- Do not normalize CPM and do not add a pseudocount.
- Match GTF `gene_id` values by exact text after trimming surrounding white space.
- Do not filter on `gene_type` or `gene_biotype`.
- Keep distinct mapped gene IDs in the full TCA view, even when they share a gene name.
- Aggregate duplicate gene names only in the dtangle view by summing linear CPM before `log2()`.
- Fit TCA once on the full mapped, nonconstant gene-ID matrix in `log2(CPM)` space.
- Extract the complete tensor once. Do not shard it.
- Name tensor matrices from `colnames(W)`.
- Write one `.bed.gz` output per retained major group with scale `log2_cpm`.
- Preserve BED coordinates, modeled-gene order, and sample order.
- Keep the ten-group mapping and remove groups with cohort mean below `0.0001`.
- Remove `hdf5r`, HDF5 outputs, optional matrix TSV outputs, tensor shards, and shard parameters from the active pipeline.
- Use tidyverse syntax in R scripts.
- Use clean minimal plots without titles or subtitles.
- Every WDL command logs stage, UTC start, dimensions, outputs, completion, and errors.
- Use the pinned Micromamba image. Do not require a local Docker build.
- GitHub Actions is the authoritative container and workflow smoke test.

## File Responsibility Map

- `R/expression_bed.R`: BED/BED.GZ input parsing, coordinate validation, and BED output writing.
- `R/expression.R`: GTF parsing, mapped gene-ID preparation, and dtangle-only gene-symbol aggregation.
- `R/tca_stage.R`: TCA validation, constant-gene removal, and one cohort-level fit. It contains no shard logic.
- `R/bed_outputs.R`: Full tensor extraction, tensor contract validation, cell-group BED writing, and output inventory creation.
- `R/qc.R`: Reconstruction and manifest functions that operate on the complete in-memory tensor and BED inventory.
- `scripts/prepare_expression.R`: Produce coordinate, gene-ID TCA, and gene-symbol dtangle preparation artifacts.
- `scripts/fit_tca.R`: Fit TCA and write model, filtered gene-ID expression, and constant-gene report.
- `scripts/export_tca_beds.R`: Extract the full tensor, write BED outputs, and write reconstruction and QC artifacts.
- `scripts/build_manifest.R`: Build provenance from the localized BED inventory.
- `workflows/tasks/*.wdl`: Expose preparation, fit, direct export, and manifest tasks with logging and Terra resources.
- `workflows/cell_type_deconvolution.wdl`: Wire both proportion modes to one full-matrix TCA and direct BED export path.
- `tests/fixtures` and `tests/smoke`: Prove coordinate, gene-ID, sample, group, scale, and duplicate-symbol behavior end to end.
- `README.md`, `docs`, `examples`, and `.dockstore.yml`: Document the final Terra interface.

---

### Task 1: Parse BED Expression and Build Separate TCA and dtangle Views

**Files:**
- Create: `R/expression_bed.R`
- Modify: `R/expression.R`
- Modify: `scripts/prepare_expression.R`
- Modify: `tests/testthat/test-expression.R`

**Interfaces:**
- Consumes: BED/BED.GZ with `#chr`, `start`, `end`, `gene_id`, then sample columns; GTF/GTF.GZ.
- Produces: `read_expression_bed(path) -> list(coordinates, cpm)`, `prepare_expression_bed(expression, annotation) -> list(tca_cpm, tca_log_expression, dtangle_cpm, dtangle_log_expression, coordinates, mapping_report, excluded_genes)`, and seven CLI artifacts.

- [ ] **Step 1: Write failing BED parsing and two-view tests**

Add tests that use two gene IDs with one shared gene name:

```r
testthat::test_that("BED preparation preserves gene IDs and collapses only dtangle symbols", {
  bed_path <- tempfile(fileext = ".bed.gz")
  input <- tibble::tibble(
    `#chr` = c("chr1", "chr2", "chr3"),
    start = c(0L, 100L, 200L),
    end = c(50L, 150L, 250L),
    gene_id = c("g1", "g2", "g3"),
    S1 = c(2, 3, 5),
    S2 = c(7, 11, 13)
  )
  readr::write_tsv(input, bed_path)
  annotation <- tibble::tibble(
    gene_id = c("g1", "g2", "g3"),
    gene_name = c("A", "A", "B"),
    gene_type = c("protein_coding", "lncRNA", NA_character_)
  )

  expression <- read_expression_bed(bed_path)
  result <- prepare_expression_bed(expression, annotation)

  testthat::expect_identical(rownames(result$tca_cpm), c("g1", "g2", "g3"))
  testthat::expect_identical(result$coordinates, input[1:4])
  testthat::expect_identical(rownames(result$dtangle_cpm), c("A", "B"))
  testthat::expect_equal(unname(result$dtangle_cpm["A", ]), c(5, 18))
  testthat::expect_equal(result$tca_log_expression, log2(result$tca_cpm))
  testthat::expect_equal(result$dtangle_log_expression, log2(result$dtangle_cpm))
})

testthat::test_that("BED validation rejects invalid coordinates and duplicate IDs", {
  invalid <- tibble::tibble(
    `#chr` = c("chr1", "chr1"), start = c(10L, 20L),
    end = c(10L, 30L), gene_id = c("g1", "g1"), S1 = c(1, 2)
  )
  path <- tempfile(fileext = ".bed")
  readr::write_tsv(invalid, path)
  testthat::expect_error(read_expression_bed(path), "start.*less than.*end")
})
```

Extend the CLI test to require these outputs:

```r
c(
  "prepared_tca_cpm.tsv.gz",
  "prepared_tca_log2_cpm.tsv.gz",
  "prepared_dtangle_cpm.tsv.gz",
  "prepared_dtangle_log2_cpm.tsv.gz",
  "prepared_coordinates.tsv",
  "gene_mapping_report.tsv",
  "excluded_genes.tsv"
)
```

- [ ] **Step 2: Run the expression tests and verify the intended failure**

Run: `Rscript tests/testthat.R tests/testthat/test-expression.R`

Expected: FAIL because `read_expression_bed()` and `prepare_expression_bed()` do not exist and the CLI still expects a gene-ID-first matrix.

- [ ] **Step 3: Implement BED parsing and validation**

Create `R/expression_bed.R` with this contract:

```r
expression_bed_columns <- function() c("#chr", "start", "end", "gene_id")

read_expression_bed <- function(path) {
  table <- readr::read_tsv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    name_repair = "minimal",
    show_col_types = FALSE,
    progress = FALSE
  )
  required <- expression_bed_columns()
  if (ncol(table) <= length(required) ||
      !identical(names(table)[seq_along(required)], required)) {
    stop(
      "Expression BED must start with #chr, start, end, and gene_id",
      call. = FALSE
    )
  }
  coordinates <- table |>
    dplyr::transmute(
      `#chr` = trimws(.data[["#chr"]]),
      start = readr::parse_integer(.data$start),
      end = readr::parse_integer(.data$end),
      gene_id = trimws(.data$gene_id)
    )
  if (anyNA(coordinates$start) || anyNA(coordinates$end) ||
      any(coordinates$start >= coordinates$end)) {
    stop("BED start must be an integer less than end", call. = FALSE)
  }
  if (any(!nzchar(coordinates[["#chr"]])) ||
      any(!nzchar(coordinates$gene_id)) ||
      anyDuplicated(coordinates$gene_id) > 0L) {
    stop("BED chromosome and gene_id values must be non-empty and unique", call. = FALSE)
  }
  sample_ids <- names(table)[-(seq_along(required))]
  if (any(!nzchar(sample_ids)) || anyDuplicated(sample_ids) > 0L) {
    stop("BED sample identifiers must be non-empty and unique", call. = FALSE)
  }
  cpm <- table |>
    dplyr::select(dplyr::all_of(sample_ids)) |>
    dplyr::mutate(dplyr::across(
      dplyr::everything(),
      ~ readr::parse_double(.x, na = character())
    )) |>
    as.matrix()
  rownames(cpm) <- coordinates$gene_id
  cpm <- validate_cpm_matrix(cpm)
  if (any(cpm <= 0)) {
    stop("CPM values must be strictly positive", call. = FALSE)
  }
  list(coordinates = coordinates, cpm = cpm)
}
```

Add `write_expression_bed(path, coordinates, matrix)` for Task 3:

```r
write_expression_bed <- function(path, coordinates, matrix) {
  required <- expression_bed_columns()
  if (!inherits(coordinates, "data.frame") ||
      !identical(names(coordinates), required)) {
    stop("coordinates must have the exact BED columns", call. = FALSE)
  }
  if (!is.matrix(matrix) || !is.numeric(matrix) ||
      !identical(rownames(matrix), coordinates$gene_id) ||
      any(!is.finite(matrix))) {
    stop("BED matrix genes and finite values must match coordinates", call. = FALSE)
  }
  output <- dplyr::bind_cols(
    tibble::as_tibble(coordinates),
    tibble::as_tibble(matrix, .name_repair = "minimal")
  )
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(output, path, na = "")
  invisible(path)
}
```

- [ ] **Step 4: Implement separate preparation views**

Refactor `R/expression.R` so `prepare_expression_bed()` keeps mapped gene IDs for TCA and uses `collapse_cpm_to_gene_names()` only for dtangle:

```r
prepare_expression_bed <- function(expression, annotation) {
  if (!is.list(expression) ||
      !all(c("coordinates", "cpm") %in% names(expression))) {
    stop("expression must contain coordinates and cpm", call. = FALSE)
  }
  cpm <- validate_cpm_matrix(expression$cpm)
  annotation <- validate_gtf_gene_annotation(annotation)
  mapping_report <- make_cpm_mapping_report(cpm, annotation)
  mapped_ids <- mapping_report |>
    dplyr::filter(!is.na(.data$gene_name), nzchar(.data$gene_name)) |>
    dplyr::pull("gene_id")
  if (length(mapped_ids) == 0L) {
    stop("No expression genes have a usable GTF gene_name", call. = FALSE)
  }
  tca_cpm <- cpm[mapped_ids, , drop = FALSE]
  coordinate_index <- match(mapped_ids, expression$coordinates$gene_id)
  coordinates <- expression$coordinates[coordinate_index, , drop = FALSE]
  dtangle_cpm <- collapse_cpm_to_gene_names(cpm, annotation)
  list(
    tca_cpm = tca_cpm,
    tca_log_expression = log2(tca_cpm),
    dtangle_cpm = dtangle_cpm,
    dtangle_log_expression = log2(dtangle_cpm),
    coordinates = coordinates,
    mapping_report = mapping_report,
    excluded_genes = make_excluded_genes(mapping_report)
  )
}
```

Change duplicate-name mapping text from `duplicate_gene_name_collapsed` to `duplicate_gene_name_aggregated_for_dtangle`. Do not classify those rows as excluded.

- [ ] **Step 5: Update the preparation CLI**

Read the BED through `read_expression_bed()`. Write:

```r
output_paths <- list(
  tca_cpm = file.path(options$output_dir, "prepared_tca_cpm.tsv.gz"),
  tca_log = file.path(options$output_dir, "prepared_tca_log2_cpm.tsv.gz"),
  dtangle_cpm = file.path(options$output_dir, "prepared_dtangle_cpm.tsv.gz"),
  dtangle_log = file.path(options$output_dir, "prepared_dtangle_log2_cpm.tsv.gz"),
  coordinates = file.path(options$output_dir, "prepared_coordinates.tsv"),
  mapping_report = file.path(options$output_dir, "gene_mapping_report.tsv"),
  excluded_genes = file.path(options$output_dir, "excluded_genes.tsv")
)
write_numeric_matrix(result$tca_cpm, output_paths$tca_cpm, "gene_id")
write_numeric_matrix(result$tca_log_expression, output_paths$tca_log, "gene_id")
write_numeric_matrix(result$dtangle_cpm, output_paths$dtangle_cpm, "gene_name")
write_numeric_matrix(result$dtangle_log_expression, output_paths$dtangle_log, "gene_name")
readr::write_tsv(result$coordinates, output_paths$coordinates, na = "")
```

Log input rows, mapped TCA rows, dtangle symbols, samples, and every output path.

- [ ] **Step 6: Run focused and lint checks**

Run: `Rscript tests/testthat.R tests/testthat/test-expression.R`

Expected: PASS.

Run: `Rscript tools/lint_r.R`

Expected: `R lint ok`.

- [ ] **Step 7: Commit Task 1**

```bash
git add R/expression_bed.R R/expression.R scripts/prepare_expression.R tests/testthat/test-expression.R
git commit -m "feat: preserve BED coordinates for TCA"
```

---

### Task 2: Remove Gene Shards from TCA Fitting

**Files:**
- Modify: `R/constants.R`
- Modify: `R/tca_stage.R`
- Modify: `scripts/fit_tca.R`
- Modify: `tests/testthat/test-tca-stage.R`

**Interfaces:**
- Consumes: gene-ID-by-sample `prepared_tca_log2_cpm.tsv.gz`, TCA weights, and optional covariates.
- Produces: `fit_tca_stage(X, W, C2, num_cores, max_iters, random_seed, log_file) -> list(model, X, excluded_genes)` and four files: model RDS, model log, filtered gene-ID expression, and constant-gene report.

- [ ] **Step 1: Replace shard tests with gene-ID fit-output tests**

Delete tests for `make_gene_shard_manifest()` and `write_gene_shards()`. Add:

```r
testthat::test_that("constant-gene removal reports gene_id and preserves order", {
  X <- matrix(
    c(1, 1, 2, 3, 4, 6),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("g1", "g2", "g3"), c("S1", "S2"))
  )
  result <- remove_constant_features(X)
  testthat::expect_identical(rownames(result$matrix), c("g2", "g3"))
  testthat::expect_identical(result$report$gene_id, "g1")
  testthat::expect_identical(result$report$reason, "constant_expression")
})

testthat::test_that("the TCA CLI does not create shard artifacts", {
  text <- paste(readLines(
    testthat::test_path("../..", "scripts", "fit_tca.R"),
    warn = FALSE
  ), collapse = "\n")
  testthat::expect_false(grepl("--shard-size", text, fixed = TRUE))
})
```

- [ ] **Step 2: Run the focused TCA test and verify failure**

Run: `Rscript tests/testthat.R tests/testthat/test-tca-stage.R`

Expected: FAIL because the report still uses `gene_name` and shard functions and CLI options remain.

- [ ] **Step 3: Remove shard code from constants and TCA stage**

Remove `tensor_shard_size` from `pipeline_defaults()`. Delete `make_gene_shard_manifest()` and `write_gene_shards()` from `R/tca_stage.R`. Change the report to:

```r
report = tibble::tibble(
  gene_id = rownames(X)[constant],
  reason = rep("constant_expression", sum(constant))
)
```

Keep all TCA 1.2.1 validation and model-fit behavior unchanged.

- [ ] **Step 4: Simplify the fit CLI**

Remove `--shard-size`, shard logging, the shard manifest, and shard text files. Read the expression matrix with `gene_id`:

```r
X <- read_numeric_matrix(options$expression_log, "gene_id")
output_paths <- list(
  model = file.path(options$output_dir, "tca_model.rds"),
  model_log = file.path(options$output_dir, "tca_model.log"),
  expression = file.path(options$output_dir, "tca_expression.tsv.gz"),
  excluded_genes = file.path(options$output_dir, "tca_excluded_genes.tsv")
)
write_numeric_matrix(result$X, output_paths$expression, "gene_id")
```

The completion log must report genes, samples, retained groups, and excluded constant genes. It must not report shards.

- [ ] **Step 5: Run focused and lint checks**

Run: `Rscript tests/testthat.R tests/testthat/test-tca-stage.R`

Expected: PASS, with dependency-based TCA skips only when TCA 1.2.1 is absent from the host.

Run: `Rscript tools/lint_r.R`

Expected: `R lint ok`.

- [ ] **Step 6: Commit Task 2**

```bash
git add R/constants.R R/tca_stage.R scripts/fit_tca.R tests/testthat/test-tca-stage.R
git commit -m "refactor: fit TCA without gene shards"
```

---

### Task 3: Extract the Full Tensor and Write Cell-Group BED Files

**Files:**
- Create: `R/bed_outputs.R`
- Create: `scripts/export_tca_beds.R`
- Create: `tests/testthat/test-bed-outputs.R`
- Modify: `R/qc.R`
- Delete: `R/tensor_outputs.R`
- Delete: `scripts/extract_tca_shard.R`
- Delete: `scripts/assemble_tca_outputs.R`
- Delete: `tests/testthat/test-tensor-outputs.R`

**Interfaces:**
- Consumes: filtered gene-ID `log2(CPM)` matrix, fitted TCA model, retained weights, prepared BED coordinates, and optional covariates.
- Produces: `extract_full_tensor(X, model, num_cores, log_file) -> named list`, `validate_tensor_contract(tensor, gene_ids, sample_ids, cell_groups)`, `write_cell_type_beds(tensor, coordinates, output_dir) -> list(paths, inventory)`, `reconstruct_tensor(tensor, weights, C2, deltas_hat) -> matrix`, reconstruction/QC files, and one BED.GZ per group.

- [ ] **Step 1: Write failing direct tensor and BED writer tests**

Create `tests/testthat/test-bed-outputs.R`:

```r
testthat::test_that("cell-type BED files preserve coordinates and sample order", {
  coordinates <- tibble::tibble(
    `#chr` = c("chr2", "chr1"), start = c(20L, 10L),
    end = c(30L, 15L), gene_id = c("g2", "g1")
  )
  tensor <- list(
    "CD4 T cells" = matrix(
      c(1, 2, 3, 4), nrow = 2, byrow = TRUE,
      dimnames = list(c("g2", "g1"), c("S2", "S1"))
    ),
    "CD8 T cells" = matrix(
      c(5, 6, 7, 8), nrow = 2, byrow = TRUE,
      dimnames = list(c("g2", "g1"), c("S2", "S1"))
    )
  )
  result <- write_cell_type_beds(tensor, coordinates, tempfile())
  observed <- readr::read_tsv(
    result$paths[["CD4 T cells"]],
    show_col_types = FALSE,
    progress = FALSE
  )
  testthat::expect_identical(
    names(observed),
    c("#chr", "start", "end", "gene_id", "S2", "S1")
  )
  testthat::expect_identical(observed$gene_id, c("g2", "g1"))
  testthat::expect_equal(as.matrix(observed[c("S2", "S1")]), tensor[[1L]])
  testthat::expect_identical(
    result$inventory$cell_group,
    c("CD4 T cells", "CD8 T cells")
  )
  testthat::expect_true(all(result$inventory$scale == "log2_cpm"))
})

testthat::test_that("tensor validation rejects group and identifier mismatches", {
  tensor <- list(A = matrix(
    1:4, nrow = 2,
    dimnames = list(c("g1", "g2"), c("S1", "S2"))
  ))
  testthat::expect_error(
    validate_tensor_contract(tensor, c("g1", "g2"), c("S1", "S2"), c("A", "B")),
    "cell groups"
  )
})
```

Add a TCA 1.2.1 test that calls `extract_full_tensor()` on `TCA::test_data()` and checks names equal `colnames(model$W)` and all matrices use full gene and sample order.

- [ ] **Step 2: Run the BED output test and verify failure**

Run: `Rscript tests/testthat.R tests/testthat/test-bed-outputs.R`

Expected: FAIL because the direct tensor and BED functions do not exist.

- [ ] **Step 3: Implement full tensor extraction and validation**

Create `R/bed_outputs.R`. Adapt the existing TCA version and logging checks, but call `TCA::tensor()` once on all rows:

```r
extract_full_tensor <- function(X, model, num_cores = 1L, log_file = NULL) {
  validate_tensor_matrix(X, "TCA expression")
  num_cores <- validate_positive_integer(num_cores, "num_cores")
  if (!is.list(model) || !is.matrix(model$W) ||
      !identical(colnames(X), rownames(model$W))) {
    stop("Model sample order must match expression sample order exactly", call. = FALSE)
  }
  validate_tensor_tca_version()
  tensor <- TCA::tensor(
    X = X,
    tca.mdl = model,
    scale = FALSE,
    parallel = num_cores > 1L,
    num_cores = num_cores,
    log_file = log_file,
    verbose = TRUE
  )
  if (!is.list(tensor) || length(tensor) != ncol(model$W)) {
    stop("TCA tensor source dimension does not match the model", call. = FALSE)
  }
  names(tensor) <- colnames(model$W)
  tensor <- purrr::map(tensor, function(source_matrix) {
    if (!identical(dim(source_matrix), dim(X))) {
      stop("TCA tensor returned an invalid gene-by-sample dimension", call. = FALSE)
    }
    dimnames(source_matrix) <- dimnames(X)
    source_matrix
  })
  validate_tensor_contract(
    tensor, rownames(X), colnames(X), colnames(model$W)
  )
  tensor
}
```

Implement exact group, gene, sample, and value validation:

```r
validate_tensor_contract <- function(
    tensor,
    gene_ids,
    sample_ids,
    cell_groups) {
  if (!is.list(tensor) || !identical(names(tensor), cell_groups)) {
    stop("Tensor cell groups must match the TCA weights exactly", call. = FALSE)
  }
  slugs <- slugify_cell_group(cell_groups)
  if (any(!nzchar(slugs)) || anyDuplicated(slugs) > 0L) {
    stop("Tensor cell groups must have unique non-empty slugs", call. = FALSE)
  }
  valid <- purrr::map_lgl(tensor, function(source_matrix) {
    is.matrix(source_matrix) && is.numeric(source_matrix) &&
      identical(rownames(source_matrix), gene_ids) &&
      identical(colnames(source_matrix), sample_ids) &&
      identical(dim(source_matrix), c(length(gene_ids), length(sample_ids))) &&
      all(is.finite(source_matrix))
  })
  if (!all(valid)) {
    stop("Tensor genes, samples, dimensions, and values must match", call. = FALSE)
  }
  invisible(TRUE)
}
```

- [ ] **Step 4: Implement direct BED writing and inventory**

Use `write_expression_bed()` from Task 1:

```r
write_cell_type_beds <- function(tensor, coordinates, output_dir) {
  validate_tensor_contract(
    tensor,
    coordinates$gene_id,
    colnames(tensor[[1L]]),
    names(tensor)
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- stats::setNames(
    file.path(output_dir, paste0(slugify_cell_group(names(tensor)), ".bed.gz")),
    names(tensor)
  )
  purrr::iwalk(paths, function(path, cell_group) {
    write_expression_bed(path, coordinates, tensor[[cell_group]])
  })
  inventory <- tibble::tibble(
    logical_name = paste0(slugify_cell_group(names(paths)), "_expression"),
    path = unname(normalizePath(paths)),
    n_genes = nrow(tensor[[1L]]),
    n_samples = ncol(tensor[[1L]]),
    scale = "log2_cpm",
    cell_group = names(paths),
    slug = slugify_cell_group(names(paths))
  )
  list(paths = paths, inventory = inventory)
}
```

- [ ] **Step 5: Implement the direct export CLI and whole-matrix QC**

Create `scripts/export_tca_beds.R` with options:

```text
--expression-log
--coordinates
--model
--weights
--covariates
--num-cores
--output-dir
--log-file
```

Read expression with `gene_id`, coordinates with exact BED columns, weights with `sample_id`, and optional covariates with `sample_id`. When `--covariates` is absent, use `model$C2`. Require coordinate gene order, weight sample order, model group order, and optional covariate sample order to match exactly.

The preparation coordinate table can still contain genes that the fit removed as constant. Select modeled rows before tensor validation:

```r
coordinate_index <- match(rownames(X), coordinates$gene_id)
if (anyNA(coordinate_index)) {
  stop("Every modeled gene_id must have BED coordinates", call. = FALSE)
}
coordinates <- coordinates[coordinate_index, , drop = FALSE]
if (!identical(coordinates$gene_id, rownames(X))) {
  stop("Coordinate gene order must match TCA expression exactly", call. = FALSE)
}
```

Extract once:

```r
tensor <- extract_full_tensor(X, model, options$num_cores, export_log_path)
bed_result <- write_cell_type_beds(tensor, coordinates, options$output_dir)
reconstructed <- reconstruct_tensor(
  tensor,
  weights,
  C2,
  if (is.null(C2)) NULL else model$deltas_hat[
    rownames(X), colnames(C2), drop = FALSE
  ]
)
statistics <- initialize_reconstruction_stats(colnames(X)) |>
  update_reconstruction_stats(X, reconstructed)
reconstruction_by_sample <- finalize_reconstruction_stats(statistics)
qc_points <- sample_qc_values(X, reconstructed)
qc_paths <- write_qc_reports(
  weights = weights,
  reconstruction_by_sample = reconstruction_by_sample,
  observed = qc_points$observed,
  reconstructed = qc_points$reconstructed,
  gene_count = nrow(X),
  output_dir = options$output_dir
)
```

Move `sample_qc_values()` from the old assembly script into `R/qc.R`. Write `cell_type_bed_inventory.tsv`, `reconstruction_by_sample.tsv`, `qc_summary.tsv`, and `qc_plots.pdf`. Log complete input and output dimensions and paths. Use a nonzero exit for every validation or write failure.

- [ ] **Step 6: Remove HDF5 and shard implementation files**

Delete `R/tensor_outputs.R`, `scripts/extract_tca_shard.R`, `scripts/assemble_tca_outputs.R`, and `tests/testthat/test-tensor-outputs.R`. Preserve reusable reconstruction, QC, hashing, and manifest functions in `R/qc.R` before deletion. Rename `reconstruct_tensor_shard()` to `reconstruct_tensor()` and update its tests so no shard-named API remains.

- [ ] **Step 7: Run direct output, QC, full R, and lint tests**

Run: `Rscript tests/testthat.R tests/testthat/test-bed-outputs.R`

Expected: PASS.

Run: `Rscript tests/testthat.R`

Expected: PASS, with dependency-based skips only when dtangle or TCA 1.2.1 is absent from the host.

Run: `Rscript tools/lint_r.R`

Expected: `R lint ok`.

- [ ] **Step 8: Commit Task 3**

```bash
git add R/bed_outputs.R R/qc.R scripts/export_tca_beds.R tests/testthat/test-bed-outputs.R
git add -u R/tensor_outputs.R scripts/extract_tca_shard.R scripts/assemble_tca_outputs.R tests/testthat/test-tensor-outputs.R
git commit -m "feat: write direct cell-type BED outputs"
```

---

### Task 4: Rewire WDL for One Direct BED Export Task

**Files:**
- Modify: `workflows/tasks/expression.wdl`
- Modify: `workflows/tasks/tca.wdl`
- Modify: `workflows/tasks/qc.wdl`
- Modify: `workflows/cell_type_deconvolution.wdl`
- Modify: `scripts/build_manifest.R`
- Modify: `tests/testthat/test-wdl-contract.R`
- Modify: `tests/testthat/test-wdl-manifest-boundary.R`

**Interfaces:**
- Consumes: Task 1–3 CLI contracts.
- Produces: one `ExportTcaBeds` call, `Array[File] cell_type_beds`, BED inventory, reconstruction/QC outputs, and a localized BED manifest. It produces no scatter, shards, HDF5, or optional matrix TSV array.

- [ ] **Step 1: Write failing WDL contract tests**

Replace active shard and HDF5 expectations with:

```r
testthat::test_that("top workflow uses one direct BED export", {
  text <- wdl_text("workflows/cell_type_deconvolution.wdl")
  testthat::expect_match(text, "call tca_tasks.ExportTcaBeds", fixed = TRUE)
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
})

testthat::test_that("manifest localizes direct BED outputs", {
  text <- wdl_text("workflows/tasks/qc.wdl")
  testthat::expect_match(text, "Array[File] cell_type_beds", fixed = TRUE)
  testthat::expect_match(text, "cell_type_bed_inventory", fixed = TRUE)
  testthat::expect_match(text, "assembled_outputs.localized.tsv", fixed = TRUE)
  testthat::expect_false(grepl("group_hdf5|group_tsv", text))
})
```

Keep the common WDL logging assertions.

- [ ] **Step 2: Run WDL contract tests and verify failure**

Run: `Rscript tests/testthat.R tests/testthat/test-wdl-contract.R tests/testthat/test-wdl-manifest-boundary.R`

Expected: FAIL because the active workflow still scatters shard extraction and exposes HDF5.

- [ ] **Step 3: Update expression and fit task contracts**

`PrepareExpression` outputs:

```wdl
File prepared_tca_cpm = "outputs/prepared_tca_cpm.tsv.gz"
File prepared_tca_log2_cpm = "outputs/prepared_tca_log2_cpm.tsv.gz"
File prepared_dtangle_cpm = "outputs/prepared_dtangle_cpm.tsv.gz"
File prepared_dtangle_log2_cpm = "outputs/prepared_dtangle_log2_cpm.tsv.gz"
File prepared_coordinates = "outputs/prepared_coordinates.tsv"
File mapping_report = "outputs/gene_mapping_report.tsv"
File excluded_genes = "outputs/excluded_genes.tsv"
File log = "prepare_expression.log"
```

`FitTca` consumes `prepared_tca_log2_cpm`. Remove `shard_size`, `shard_manifest`, and `shards`. Delete `ExtractTcaShard`.

- [ ] **Step 4: Add the direct export WDL task**

Add `ExportTcaBeds` to `workflows/tasks/tca.wdl`:

```wdl
task ExportTcaBeds {
  input {
    File tca_expression
    File coordinates
    File model
    File tca_weights
    File? covariates
    String docker_image
    Int cpu = 8
    String memory = "128 GB"
    Int disk_gb = 500
    Int preemptible_attempts = 0
    Int max_retries = 1
  }

  String covariates_argument = if defined(covariates) then "--covariates " + select_first([covariates]) else ""

  command <<<
    set -euo pipefail
    stage="export_tca_beds"
    log="$stage.log"
    status=0
    printf 'stage=%s start_time=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
    trap 'status=$?; printf "stage=%s error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    Rscript /opt/celltype/scripts/export_tca_beds.R \
      --expression-log '~{tca_expression}' \
      --coordinates '~{coordinates}' \
      --model '~{model}' \
      --weights '~{tca_weights}' \
      ~{covariates_argument} \
      --num-cores '~{cpu}' \
      --output-dir outputs \
      --log-file outputs/export_tca_beds.log 2>&1 | tee -a "$log"
    printf 'stage=%s dimensions=%s outputs=%s completion_time=%s\n' \
      "$stage" "$(wc -l < outputs/cell_type_bed_inventory.tsv)" \
      "cell_type_beds,inventory,reconstruction,qc_summary,qc_plots" \
      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>

  output {
    Array[File] cell_type_beds = glob("outputs/*.bed.gz")
    File cell_type_bed_inventory = "outputs/cell_type_bed_inventory.tsv"
    File reconstruction_by_sample = "outputs/reconstruction_by_sample.tsv"
    File qc_summary = "outputs/qc_summary.tsv"
    File qc_plots = "outputs/qc_plots.pdf"
    File export_detail_log = "outputs/export_tca_beds.log"
    File log = "export_tca_beds.log"
  }

  runtime {
    docker: docker_image
    cpu: cpu
    memory: memory
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: preemptible_attempts
    maxRetries: max_retries
  }
}
```

- [ ] **Step 5: Update manifest localization and top-level wiring**

Change `BuildManifest` to accept `cell_type_beds`, `cell_type_bed_inventory`, `export_qc_summary`, and `export_qc_plots`. Stage each BED by safe basename. Rewrite only the inventory `path` column to the localized files before calling `scripts/build_manifest.R`. Replace the script's fallback parameter object with `list(scale = "log2_cpm")`; the top workflow will continue to supply the complete effective-parameter JSON.

In the top workflow:

- Send `PrepareExpression.prepared_dtangle_log2_cpm` to `RunDtangle`.
- Send `PrepareExpression.prepared_tca_log2_cpm` to `FitTca`.
- Call `ExportTcaBeds` once with `FitTca.tca_expression`, `PrepareExpression.prepared_coordinates`, the model, weights, and covariates.
- Send the export inventory and BED array to `BuildManifest`.
- Expose `cell_type_beds`, `cell_type_bed_inventory`, reconstruction, QC, model, proportions, reports, logs, manifest, and effective parameters.
- Use `export_log` and `export_detail_log` output names. Remove assembly-named top-level outputs.

Remove these workflow inputs and effective-parameter fields:

```text
tca_shard_size
write_tsv
extract_cpu
extract_memory
extract_disk_gb
extract_preemptible_attempts
extract_max_retries
assemble_cpu
assemble_memory
assemble_disk_gb
assemble_preemptible_attempts
assemble_max_retries
hdf5_gene_chunk_max
hdf5_sample_chunk_max
hdf5_gzip_level
```

Add `export_cpu`, `export_memory`, `export_disk_gb`, `export_preemptible_attempts`, and `export_max_retries` with the values in `ExportTcaBeds` above.

- [ ] **Step 6: Run WDL tests and static validation**

Run: `Rscript tests/testthat.R tests/testthat/test-wdl-contract.R tests/testthat/test-wdl-manifest-boundary.R`

Expected: PASS.

Run: `miniwdl check workflows/cell_type_deconvolution.wdl`

Expected: PASS.

Run: `python3 tools/check_wdl_logging.py workflows`

Expected: all task command blocks pass.

- [ ] **Step 7: Commit Task 4**

```bash
git add workflows scripts/build_manifest.R tests/testthat/test-wdl-contract.R tests/testthat/test-wdl-manifest-boundary.R
git commit -m "feat: wire direct BED export in WDL"
```

---

### Task 5: Replace Smoke Fixtures and Assertions with BED Contracts

**Files:**
- Modify: `scripts/generate_synthetic_fixture.R`
- Modify: `tests/fixtures/dtangle.inputs.json`
- Modify: `tests/fixtures/restart.inputs.json`
- Create: `tests/fixtures/synthetic_expression.bed`
- Create: `tests/fixtures/expected_gene_ids.txt`
- Create: `tests/fixtures/expected_coordinates.tsv`
- Modify: `tests/fixtures/synthetic.gtf`
- Modify: `tests/smoke/assert_outputs.R`
- Modify: `tests/testthat/test-wdl-contract.R`
- Delete: `tests/fixtures/synthetic_cpm.tsv`
- Delete: `tests/fixtures/expected_genes.txt`

**Interfaces:**
- Consumes: final WDL contract from Task 4.
- Produces: deterministic dtangle and restart fixtures that prove coordinate preservation, duplicate-symbol separation, direct BED outputs, and full-matrix TCA.

- [ ] **Step 1: Write a failing fixture contract test**

Add:

```r
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
```

- [ ] **Step 2: Run the fixture contract and verify failure**

Run: `Rscript tests/testthat.R tests/testthat/test-wdl-contract.R`

Expected: FAIL because the BED fixture does not exist.

- [ ] **Step 3: Generate a deterministic BED fixture**

Change the generator so it writes:

```r
coordinates <- tibble::tibble(
  `#chr` = rep(c("chr1", "chr2"), length.out = nrow(expression_cpm)),
  start = as.integer((seq_len(nrow(expression_cpm)) - 1L) * 1000L),
  end = as.integer((seq_len(nrow(expression_cpm)) - 1L) * 1000L + 500L),
  gene_id = rownames(expression_cpm)
)
expression_bed <- dplyr::bind_cols(
  coordinates,
  tibble::as_tibble(expression_cpm, .name_repair = "minimal")
)
readr::write_tsv(
  expression_bed,
  file.path(output_directory, "synthetic_expression.bed")
)
```

Add one extra gene ID that maps to the same gene name as one LM22 signature gene before CPM scaling:

```r
duplicate_gene_id <- "ENSGDUP000001"
duplicate_gene_name <- signature_gene_names[[1L]]
duplicate_linear <- matrix(
  stats::runif(length(sample_ids), min = 0.75, max = 3.25),
  nrow = 1L,
  dimnames = list(duplicate_gene_name, sample_ids)
)
bulk_linear <- rbind(
  signature %*% t(weights),
  extra_gene_expression,
  duplicate_linear
)
all_gene_ids <- c(signature_gene_ids, extra_gene_ids, duplicate_gene_id)
all_gene_names <- c(signature_gene_names, extra_gene_names, duplicate_gene_name)
bulk_cpm <- sweep(bulk_linear, 2L, colSums(bulk_linear), "/") * 1e6
expression_cpm <- bulk_cpm
rownames(expression_cpm) <- all_gene_ids
```

Give the two rows different coordinates through their different row indices. Expected TCA gene IDs must include both rows. The expected dtangle shared-symbol matrix must contain the shared LM22 symbol once.

Write expected gene IDs in input order and exact expected coordinate rows. Update both JSON files to use `synthetic_expression.bed`, remove shard, TSV, extract, and assemble inputs, and add the export resource inputs.

- [ ] **Step 4: Replace HDF5 smoke assertions with BED assertions**

Remove the `hdf5r` requirement. Read `cell_type_bed_inventory` and `cell_type_beds`. For every file:

```r
bed <- readr::read_tsv(
  path,
  col_types = readr::cols(.default = readr::col_character()),
  name_repair = "minimal",
  show_col_types = FALSE,
  progress = FALSE
)
stopifnot(
  identical(names(bed)[1:4], c("#chr", "start", "end", "gene_id")),
  identical(bed$gene_id, expected_gene_ids),
  identical(names(bed)[-(1:4)], expected_samples),
  identical(bed[["#chr"]], expected_coordinates[["#chr"]]),
  identical(readr::parse_integer(bed$start), expected_coordinates$start),
  identical(readr::parse_integer(bed$end), expected_coordinates$end)
)
values <- bed[-(1:4)] |>
  dplyr::mutate(dplyr::across(dplyr::everything(), readr::parse_double)) |>
  as.matrix()
stopifnot(all(is.finite(values)))
```

Assert one BED per expected group, unique group slugs, `log2_cpm` inventory scale, full gene count, and exact sample count. Assert the TCA expression uses `gene_id`, includes both duplicate-symbol gene IDs, and includes genes outside LM22. Assert the dtangle shared bulk remains gene-symbol keyed and contains the shared symbol once.

The manifest parameters must contain only the final scientific fields and `scale`. It must not contain shard, HDF5, or optional-TSV parameters.

- [ ] **Step 5: Regenerate fixtures and run smoke static checks**

Run: `Rscript scripts/generate_synthetic_fixture.R tests/fixtures`

Expected: fixture files and input JSON are deterministic.

Run the generator a second time, then run: `git diff --exit-code -- tests/fixtures`

Expected: no new diff after the second generation.

Run: `Rscript tests/testthat.R tests/testthat/test-wdl-contract.R`

Expected: PASS.

Run: `miniwdl check workflows/cell_type_deconvolution.wdl`

Expected: PASS.

- [ ] **Step 6: Commit Task 5**

```bash
git add scripts/generate_synthetic_fixture.R tests/fixtures tests/smoke/assert_outputs.R tests/testthat/test-wdl-contract.R
git add -u tests/fixtures/synthetic_cpm.tsv tests/fixtures/expected_genes.txt
git commit -m "test: smoke test coordinate-preserving BED outputs"
```

---

### Task 6: Remove HDF5 Dependency and Revise User Documentation

**Files:**
- Modify: `envs/environment.yml`
- Modify: `envs/Dockerfile`
- Modify: `README.md`
- Modify: `docs/terra.md`
- Modify: `docs/data-dictionary.md`
- Create: `examples/bed.inputs.json`
- Delete: `examples/cpm.inputs.json`
- Modify: `examples/precomputed-proportions.inputs.json`
- Modify: `.dockstore.yml`
- Modify: `tests/testthat/test-documentation.R`
- Modify: `tests/testthat/test-wdl-contract.R`

**Interfaces:**
- Consumes: final WDL names, resources, fixtures, and BED output schema.
- Produces: pinned image without `hdf5r`, Terra-ready BED examples, and exact BED input/output documentation.

- [ ] **Step 1: Write failing dependency and documentation tests**

Add:

```r
testthat::test_that("active dependencies and guides use direct BED outputs", {
  environment <- paste(readLines(
    testthat::test_path("../..", "envs", "environment.yml"),
    warn = FALSE
  ), collapse = "\n")
  guides <- paste(vapply(
    c("README.md", "docs/terra.md", "docs/data-dictionary.md"),
    function(path) paste(readLines(
      testthat::test_path("../..", path), warn = FALSE
    ), collapse = "\n"),
    character(1)
  ), collapse = "\n")
  testthat::expect_false(grepl("r-hdf5r", environment, fixed = TRUE))
  testthat::expect_match(guides, "#chr.*start.*end.*gene_id")
  testthat::expect_match(guides, "(?i)one.*BED[.]GZ.*retained")
  testthat::expect_match(guides, "log2_cpm", fixed = TRUE)
  testthat::expect_match(guides, "does not produce HDF5", fixed = TRUE)
  testthat::expect_false(grepl(
    "group_hdf5|tensor_shards|gene_shard_manifest|write_tsv",
    guides
  ))
})
```

Require both example JSON files to use the exact `expression` and `gtf` WDL inputs and no removed parameters.

- [ ] **Step 2: Run documentation tests and verify failure**

Run: `Rscript tests/testthat.R tests/testthat/test-documentation.R tests/testthat/test-wdl-contract.R`

Expected: FAIL because `r-hdf5r`, HDF5 documentation, and removed parameters remain.

- [ ] **Step 3: Remove the pinned HDF5 dependency**

Delete this exact dependency from `envs/environment.yml`:

```yaml
  - r-hdf5r=1.3.12
```

Keep all other package versions and Micromamba image pins unchanged. Remove any Dockerfile version check or label that refers to `hdf5r`. Keep explicit checks for dtangle 2.0.10 and TCA 1.2.1.

- [ ] **Step 4: Revise README and Terra guidance**

Use short, direct ASD-STE100-style sentences. Document:

- Exact BED/BED.GZ input columns and positive linear CPM values.
- Required GTF and no gene-type filter.
- The separate gene-ID TCA and gene-symbol dtangle views.
- Exact trimmed ID matching and duplicate-symbol behavior.
- Both proportion modes, ratio markers, ten groups, and `0.0001` filter.
- One full-matrix TCA fit and one direct tensor extraction.
- One BED.GZ per retained group on the `log2_cpm` model scale.
- Coordinate, modeled-gene, and sample-order preservation.
- Exclusion of unmapped and constant genes.
- Direct export resource defaults for about 9,000 samples.
- Manifest, QC, logs, GitHub Actions smoke coverage, and LM22 licensing.
- A migration note: the workflow does not produce HDF5 or tensor shards.
- Existing scientific limitations, including LM22 platform limits, no erythrocyte or platelet model, cohort-specific groups, small-proportion instability, and statistical TCA estimates.

Do not add plot titles or subtitles.

- [ ] **Step 5: Revise the exact data dictionary and examples**

For each final input and output, list exact WDL name, WDL type, mode, default, scale, validation, and description. Remove all shard, HDF5, write-TSV, extract, and assemble fields. Add the five export resource inputs and these outputs:

```text
prepared_tca_cpm
prepared_tca_log2_cpm
prepared_dtangle_cpm
prepared_dtangle_log2_cpm
prepared_coordinates
cell_type_beds
cell_type_bed_inventory
reconstruction_by_sample
qc_summary
qc_plots
export_detail_log
export_log
```

Replace `examples/cpm.inputs.json` with `examples/bed.inputs.json`. Update both examples to use cloud BED paths and the final resource fields. Keep exact namespace `cell_type_deconvolution`.

Update `.dockstore.yml` only if the fixture test parameter path changed.

- [ ] **Step 6: Run documentation, environment, WDL, and lint checks**

Run: `Rscript tests/testthat.R tests/testthat/test-documentation.R tests/testthat/test-wdl-contract.R`

Expected: PASS.

Run: `miniwdl check workflows/cell_type_deconvolution.wdl`

Expected: PASS.

Run: `Rscript tools/lint_r.R`

Expected: `R lint ok`.

Run: `git diff --check`

Expected: no output.

- [ ] **Step 7: Commit Task 6**

```bash
git add envs README.md docs examples .dockstore.yml tests/testthat/test-documentation.R tests/testthat/test-wdl-contract.R
git add -u examples/cpm.inputs.json
git commit -m "docs: publish direct BED output contract"
```

---

### Task 7: Final Local Verification and GitHub Actions Handoff

**Files:**
- Modify only files required by verified failures.

**Interfaces:**
- Consumes: complete BED-output implementation.
- Produces: clean local evidence and, after explicit user approval to push, passing GitHub Actions evidence for both workflow modes.

- [ ] **Step 1: Run the complete local R suite**

Run: `Rscript tests/testthat.R`

Expected: PASS. On a host without dtangle 2.0.10 or TCA 1.2.1, only tests with explicit dependency guards may skip. No HDF5 dependency skip is permitted.

- [ ] **Step 2: Run all static checks**

Run: `Rscript tools/lint_r.R`

Expected: `R lint ok`.

Run: `miniwdl check workflows/cell_type_deconvolution.wdl`

Expected: PASS.

Run: `python3 tools/check_wdl_logging.py workflows`

Expected: every WDL command block passes.

Run: `git diff --check`

Expected: no output.

- [ ] **Step 3: Confirm removal of active HDF5 and shard implementation**

Run:

```bash
rg -n -i 'hdf5|hdf5r|tensor_shard|gene_shard|shard_size|write_tsv|group_hdf5' \
  R scripts workflows envs tests README.md docs/terra.md docs/data-dictionary.md examples
```

Expected: no active implementation contract. Negative tests and a README or guide migration sentence that says HDF5 and sharding were removed are allowed. Historical files under `docs/superpowers` are outside this scan.

- [ ] **Step 4: Confirm no LM22 matrix or secrets are tracked**

Run: `git ls-files | rg -i 'lm22[.]txt|token|credential|secret'`

Expected: no LM22 data matrix, token, credential, or secret file. Review path-name matches without printing file contents.

- [ ] **Step 5: Inspect the final branch diff and status**

Run: `git status --short`

Expected: no output.

Run: `git diff --check main...HEAD`

Expected: no output.

Run: `git diff --stat main...HEAD`

Expected: only pipeline, tests, fixtures, CI, environment, documentation, specifications, and plans in scope.

- [ ] **Step 6: Stop and request approval before the GitHub push**

Do not run `git push` without explicit user approval. Report the branch name, commit range, local checks, dependency-based skips, and the exact GitHub workflow that will run.

- [ ] **Step 7: After approval, push and wait for pipeline CI**

Run: `git push -u origin feat/dtangle-tca-pipeline`

Run: `gh run list --workflow pipeline-ci.yml --branch feat/dtangle-tca-pipeline --limit 1`

Expected: one completed successful run.

- [ ] **Step 8: Record authoritative runtime evidence**

Run:

```bash
gh run view "$(gh run list --workflow pipeline-ci.yml --branch feat/dtangle-tca-pipeline --limit 1 --json databaseId --jq '.[0].databaseId')" --log-failed
```

Expected: no failed logs. The successful job must show the Micromamba image build, full R tests with active TCA 1.2.1 tensor tests, R lint, WDL validation, logging validation, dtangle smoke, restart smoke, and BED output assertions.

- [ ] **Step 9: Commit only corrections proven by verification**

If a check exposes a repository defect, add a focused regression test, implement the smallest correction, rerun the focused test and all Steps 1–5, then commit:

```bash
git add R scripts workflows tests envs .github README.md docs examples .dockstore.yml tools
git commit -m "fix: resolve BED pipeline verification failures"
```

If verification changes no file, do not create a verification commit.
