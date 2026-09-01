# Whole-Blood dtangle and TCA Pipeline Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a tested Terra WDL workflow that estimates LM22 cell proportions with dtangle, combines and filters the proportions, and uses TCA to produce sample-specific expression matrices for each retained major cell group.

**Architecture:** Focused R modules implement matrix input, expression preparation, dtangle, proportion processing, TCA, HDF5 assembly, and quality control. Small command-line R scripts expose those modules to separate WDL tasks. One top-level WDL runs dtangle or accepts precomputed proportions, fits one cohort-wide TCA model, scatters tensor extraction by gene, and assembles one HDF5 matrix per retained group.

**Tech Stack:** WDL 1.1, miniwdl 1.15.0, R 4.5.3, tidyverse 2.0.0, dtangle 2.0.10, TCA 1.2.1, hdf5r 1.3.12, testthat 3.3.2, micromamba 2.9.0, GitHub Actions, GHCR, Terra, and Dockstore.

**Spec:** `docs/superpowers/specs/2026-09-01-whole-blood-dtangle-tca-design.md`

## Global Constraints

- The target data are human whole-blood RNA-seq from approximately 9,000 samples.
- Accept a gene-by-sample matrix as raw counts or TPM.
- Require `gene_id`, `gene_symbol`, and `gene_length_bp` annotation for count input.
- Accept the standard linear LM22 matrix at run time and do not redistribute LM22.
- Use `log2(TPM + 1)` for the bulk mixtures and `log2(LM22)` for the reference.
- Stop when any LM22 value is zero, negative, missing, or nonfinite.
- Keep joint quantile normalization off by default. When enabled, use `limma::normalizeBetweenArrays()` once on the joined log-scale matrix.
- Run `dtangle()` with all 22 LM22 types, `data_type = "rna-seq"`, `marker_method = "ratio"`, and a marker fraction of `0.10`.
- Combine the 22 proportions into the nine approved groups after dtangle.
- Remove every group whose cohort mean is less than `0.0001`; no group is mandatory.
- Replace exact zeros in retained groups with `1e-6`, then normalize each sample row to one.
- Fit one cohort-wide TCA model with `refit_W = FALSE` on `log2(TPM + 1)`.
- Label TCA outputs as `log2_tpm_plus_1`, not TPM.
- Use a default tensor shard size of 500 genes and preserve deterministic gene and sample order.
- Use random seed `20260901` for every R stage that can call a stochastic package operation.
- Use tidyverse syntax in R where it improves clarity.
- Use clean, minimal plots without titles or subtitles.
- Add stage, start-time, completion-time, dimensions, outputs, and error context to every WDL command log.
- Use a micromamba base image and exact package pins from conda-forge and bioconda.
- Build and run the end-to-end smoke test in GitHub Actions. Do not require a local Docker build.
- Execute the plan on branch `feat/dtangle-tca-pipeline` in the isolated worktree created by the execution skill.

## Planned File Structure

```text
R/
  constants.R                LM22 names, group map, and pipeline defaults
  io.R                       TSV, JSON, sample-ID, and matrix I/O
  expression.R               count-to-TPM and gene-symbol processing
  dtangle_stage.R            LM22 validation, alignment, and dtangle fitting
  proportions.R              grouping, filtering, zero floor, and normalization
  tca_stage.R                TCA validation, fitting, and gene-shard creation
  tensor_outputs.R           tensor extraction and HDF5 assembly
  qc.R                       reconstruction statistics, plots, and manifest
scripts/
  bootstrap.R                resolve the project root and source R modules
  prepare_expression.R       expression CLI
  run_dtangle.R              dtangle CLI
  process_proportions.R      proportion CLI
  fit_tca.R                  TCA fit and shard-list CLI
  extract_tca_shard.R        tensor extraction CLI
  assemble_tca_outputs.R     HDF5 assembly and reconstruction CLI
  build_manifest.R           final QC and provenance CLI
  generate_synthetic_fixture.R  deterministic non-LM22 test-data generator
tests/
  testthat.R                 testthat entry point
  testthat/helper-load.R     source modules and shared helpers
  testthat/test-io.R
  testthat/test-expression.R
  testthat/test-dtangle-stage.R
  testthat/test-proportions.R
  testthat/test-tca-stage.R
  testthat/test-tensor-outputs.R
  testthat/test-wdl-contract.R
  fixtures/                  generated synthetic expression, signature, and inputs
  smoke/assert_outputs.R     end-to-end output assertions
workflows/
  cell_type_deconvolution.wdl
  tasks/expression.wdl
  tasks/dtangle.wdl
  tasks/proportions.wdl
  tasks/tca.wdl
  tasks/qc.wdl
envs/
  environment.yml
  Dockerfile
examples/
  counts.inputs.json
  tpm.inputs.json
  precomputed-proportions.inputs.json
docs/
  terra.md
  data-dictionary.md
.dockstore.yml
.github/workflows/
  pipeline-ci.yml
  docker-image.yml
```

---

### Task 1: Core Constants, Matrix I/O, and Test Harness

**Files:**
- Create: `R/constants.R`
- Create: `R/io.R`
- Create: `scripts/bootstrap.R`
- Create: `tests/testthat.R`
- Create: `tests/testthat/helper-load.R`
- Create: `tests/testthat/test-io.R`

**Interfaces:**
- Consumes: TSV or TSV.GZ files whose first column contains row identifiers.
- Produces: `lm22_cell_types() -> character(22)`, `lm22_group_map() -> named list(9)`, `slugify_cell_group(x) -> character`, `read_numeric_matrix(path, id_column) -> numeric matrix`, `write_numeric_matrix(x, path, id_column) -> invisible(path)`, and `assert_identical_ids(expected, observed, label) -> invisible(TRUE)`.

- [ ] **Step 1: Write the failing I/O and constant tests**

```r
testthat::test_that("the group map partitions all LM22 types once", {
  members <- unlist(lm22_group_map(), use.names = FALSE)
  testthat::expect_setequal(members, lm22_cell_types())
  testthat::expect_length(members, 22L)
  testthat::expect_false(anyDuplicated(members) > 0L)
})

testthat::test_that("numeric matrix round-trip preserves IDs and values", {
  x <- matrix(c(1, 2, 3, 4), nrow = 2,
              dimnames = list(c("G1", "G2"), c("S1", "S2")))
  path <- tempfile(fileext = ".tsv.gz")
  write_numeric_matrix(x, path, "gene_symbol")
  observed <- read_numeric_matrix(path, "gene_symbol")
  testthat::expect_identical(dimnames(observed), dimnames(x))
  testthat::expect_equal(observed, x)
})

testthat::test_that("sample mismatch reports the first difference", {
  testthat::expect_error(
    assert_identical_ids(c("S1", "S2"), c("S2", "S1"), "proportions"),
    "proportions.*S1.*S2"
  )
})
```

- [ ] **Step 2: Run the tests and verify the expected failure**

Run: `Rscript tests/testthat.R tests/testthat/test-io.R`

Expected: FAIL because `lm22_group_map()` and `read_numeric_matrix()` do not exist.

- [ ] **Step 3: Implement the constants and strict matrix I/O**

Use these exact group names and members in `R/constants.R`:

```r
lm22_group_map <- function() {
  list(
    "B cells" = c("B cells naive", "B cells memory", "Plasma cells"),
    "CD4 T cells" = c(
      "T cells CD4 naive", "T cells CD4 memory resting",
      "T cells CD4 memory activated", "T cells follicular helper",
      "T cells regulatory (Tregs)"
    ),
    "CD8 T cells" = "T cells CD8",
    "NK cells" = c("NK cells resting", "NK cells activated"),
    "Monocyte/myeloid" = c(
      "Monocytes", "Macrophages M0", "Macrophages M1", "Macrophages M2"
    ),
    "Neutrophils" = "Neutrophils",
    "Eosinophils" = "Eosinophils",
    "Dendritic cells" = c("Dendritic cells resting", "Dendritic cells activated"),
    "Other LM22" = c(
      "T cells gamma delta", "Mast cells resting", "Mast cells activated"
    )
  )
}

lm22_cell_types <- function() {
  c(
    "B cells naive", "B cells memory", "Plasma cells", "T cells CD8",
    "T cells CD4 naive", "T cells CD4 memory resting",
    "T cells CD4 memory activated", "T cells follicular helper",
    "T cells regulatory (Tregs)", "T cells gamma delta",
    "NK cells resting", "NK cells activated", "Monocytes",
    "Macrophages M0", "Macrophages M1", "Macrophages M2",
    "Dendritic cells resting", "Dendritic cells activated",
    "Mast cells resting", "Mast cells activated", "Eosinophils", "Neutrophils"
  )
}

pipeline_defaults <- function() {
  list(
    min_lm22_overlap = 0.80,
    marker_fraction = 0.10,
    group_mean_threshold = 0.0001,
    zero_floor = 1e-6,
    tensor_shard_size = 500L
  )
}

slugify_cell_group <- function(x) {
  x |>
    stringr::str_to_lower() |>
    stringr::str_replace_all("[^a-z0-9]+", "_") |>
    stringr::str_remove("^_") |>
    stringr::str_remove("_$")
}
```

Implement `read_numeric_matrix()` with `readr::read_tsv()`, explicit first-column lookup, duplicate-ID checks, numeric-column checks, finite-value checks, and stable input order. Implement `write_numeric_matrix()` with `tibble::rownames_to_column()` and `readr::write_tsv()`.

Implement `scripts/bootstrap.R` so every CLI resolves and sources the same modules:

```r
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]))
pipeline_root <- normalizePath(Sys.getenv(
  "CELLTYPE_ROOT",
  unset = file.path(dirname(script_path), "..")
))
r_files <- list.files(file.path(pipeline_root, "R"), pattern = "[.]R$", full.names = TRUE)
invisible(lapply(sort(r_files), source))
```

Implement `tests/testthat.R` so the focused commands in this plan work:

```r
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L) {
  testthat::test_dir("tests/testthat", reporter = "summary")
} else {
  purrr::walk(args, ~ testthat::test_file(.x, reporter = "summary"))
}
```

- [ ] **Step 4: Run the focused tests and lint**

Run: `Rscript tests/testthat.R tests/testthat/test-io.R`

Expected: PASS.

Run: `Rscript tools/lint_r.R`

Expected: `R lint ok`.

- [ ] **Step 5: Commit the core interfaces**

```bash
git add R/constants.R R/io.R scripts/bootstrap.R tests/testthat.R tests/testthat/helper-load.R tests/testthat/test-io.R
git commit -m "feat: add pipeline constants and matrix IO"
```

---

### Task 2: Expression Preparation

**Files:**
- Create: `R/expression.R`
- Create: `scripts/prepare_expression.R`
- Create: `tests/testthat/test-expression.R`

**Interfaces:**
- Consumes: A nonnegative gene-by-sample matrix, `expression_type` equal to `counts` or `tpm`, and optional annotation with `gene_id`, `gene_symbol`, and `gene_length_bp`.
- Produces: `read_optional_annotation(path) -> tibble or NULL`, `counts_to_tpm(counts, gene_length_bp) -> numeric matrix`, `collapse_to_symbols(values, annotation) -> numeric matrix`, `prepare_expression(expression, expression_type, annotation = NULL) -> list(tpm, log_expression, mapping_report, excluded_genes)`, and CLI files `prepared_tpm.tsv.gz`, `prepared_log2_tpm_plus_1.tsv.gz`, `gene_mapping_report.tsv`, and `excluded_genes.tsv`.

- [ ] **Step 1: Write failing count, TPM, mapping, and validation tests**

```r
testthat::test_that("counts convert to TPM by gene length", {
  counts <- matrix(c(100, 100, 200, 100), nrow = 2,
                   dimnames = list(c("g1", "g2"), c("s1", "s2")))
  observed <- counts_to_tpm(counts, c(g1 = 1000, g2 = 2000))
  testthat::expect_equal(colSums(observed), c(s1 = 1e6, s2 = 1e6))
  testthat::expect_equal(observed[, "s1"], c(g1 = 2 / 3 * 1e6, g2 = 1 / 3 * 1e6))
})

testthat::test_that("duplicate symbols are summed and renormalized", {
  values <- matrix(c(1, 2, 3, 4), nrow = 2,
                   dimnames = list(c("ENSG1", "ENSG2"), c("s1", "s2")))
  annotation <- tibble::tibble(
    gene_id = c("ENSG1", "ENSG2"), gene_symbol = c("A", "A"),
    gene_length_bp = c(1000, 1000)
  )
  observed <- collapse_to_symbols(values, annotation)
  testthat::expect_identical(rownames(observed), "A")
  testthat::expect_equal(as.numeric(observed), c(3, 7))
})

testthat::test_that("log expression is log2 TPM plus one", {
  x <- matrix(c(0, 1e6), nrow = 2,
              dimnames = list(c("A", "B"), "s1"))
  result <- prepare_expression(x, "tpm")
  testthat::expect_equal(result$log_expression, log2(result$tpm + 1))
})

testthat::test_that("the expression CLI writes its four declared outputs", {
  input <- tempfile(fileext = ".tsv")
  output_dir <- tempfile()
  dir.create(output_dir)
  x <- matrix(c(0, 1e6), nrow = 2,
              dimnames = list(c("A", "B"), "s1"))
  write_numeric_matrix(x, input, "gene_id")
  status <- system2(
    "Rscript",
    c(
      "scripts/prepare_expression.R", "--expression", input,
      "--expression-type", "tpm", "--output-dir", output_dir
    )
  )
  testthat::expect_equal(status, 0L)
  testthat::expect_true(all(file.exists(file.path(output_dir, c(
    "prepared_tpm.tsv.gz", "prepared_log2_tpm_plus_1.tsv.gz",
    "gene_mapping_report.tsv", "excluded_genes.tsv"
  ))))
})
```

Add tests that reject negative values, missing gene lengths, nonpositive gene lengths, an invalid expression type, and duplicated sample identifiers.

- [ ] **Step 2: Run the focused tests and verify failure**

Run: `Rscript tests/testthat.R tests/testthat/test-expression.R`

Expected: FAIL because `counts_to_tpm()` is not defined.

- [ ] **Step 3: Implement expression preparation**

Use sample-wise rate normalization:

```r
counts_to_tpm <- function(counts, gene_length_bp) {
  stopifnot(identical(rownames(counts), names(gene_length_bp)))
  rates <- counts / (gene_length_bp / 1000)
  totals <- colSums(rates)
  if (any(totals <= 0)) stop("Each sample must have a positive total expression rate.")
  sweep(rates, 2, totals, "/") * 1e6
}
```

Implement symbol collapse with a long tibble, `dplyr::group_by(gene_symbol, sample_id)`, `dplyr::summarise(value = sum(value))`, and `tidyr::pivot_wider()`. After symbol collapse, normalize each TPM column to one million. Record unmapped genes, duplicate-symbol actions, constant genes, and the reason for each exclusion.

Implement `read_optional_annotation()` so an empty string or `NULL` returns `NULL`; otherwise, read a TSV and require `gene_id` and `gene_symbol`. Make `prepare_expression()` require `gene_length_bp` when `expression_type` is `counts`.

- [ ] **Step 4: Add the command-line wrapper**

Use `optparse` options `--expression`, `--expression-type`, `--annotation`, and `--output-dir`. Log the stage name, UTC start time, input dimensions, output paths, and UTC completion time with `message()`.

```r
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]))
source(file.path(dirname(script_path), "bootstrap.R"))

result <- prepare_expression(
  expression = read_numeric_matrix(options$expression, "gene_id"),
  expression_type = options$expression_type,
  annotation = read_optional_annotation(options$annotation)
)
```

- [ ] **Step 5: Run focused tests and a CLI round trip**

Run: `Rscript tests/testthat.R tests/testthat/test-expression.R`

Expected: PASS, including the CLI test that creates temporary input and checks all four outputs.

- [ ] **Step 6: Commit expression preparation**

```bash
git add R/expression.R scripts/prepare_expression.R tests/testthat/test-expression.R
git commit -m "feat: prepare count and TPM expression inputs"
```

---

### Task 3: LM22 Validation and dtangle Estimation

**Files:**
- Create: `R/dtangle_stage.R`
- Create: `scripts/run_dtangle.R`
- Create: `tests/testthat/test-dtangle-stage.R`

**Interfaces:**
- Consumes: `log2(TPM + 1)` bulk expression, the standard positive linear LM22 matrix, minimum overlap, marker fraction, and a quantile-normalization flag.
- Produces: `validate_lm22(lm22_linear) -> invisible(TRUE)`, `transform_lm22(lm22_linear) -> numeric matrix`, `prepare_dtangle_inputs(bulk_log, lm22_linear, min_overlap, quantile_normalize) -> list(Y, references, transformed_lm22, shared_bulk, overlap_report)`, `estimate_dtangle(inputs, marker_fraction) -> list(proportions, markers, metadata)`, and CLI outputs for those objects.

- [ ] **Step 1: Write failing LM22 scale and structure tests**

```r
testthat::test_that("standard LM22 is logged without a pseudocount", {
  lm22 <- matrix(
    rep(c(0.25, 4, 16), each = 22), nrow = 3, byrow = TRUE,
    dimnames = list(c("G1", "G2", "G3"), lm22_cell_types())
  )
  observed <- transform_lm22(lm22)
  testthat::expect_equal(observed, log2(lm22))
})

testthat::test_that("LM22 rejects zero and missing cell types", {
  lm22 <- matrix(1, nrow = 3, ncol = 21,
                 dimnames = list(c("G1", "G2", "G3"), lm22_cell_types()[-1]))
  testthat::expect_error(validate_lm22(lm22), "22 standard LM22 columns")
})
```

Add tests for duplicate genes, negative values, zero values, nonfinite values, and overlap below 0.80.

- [ ] **Step 2: Write the failing dtangle integration test**

Build a deterministic synthetic reference with three unique markers for every LM22 type:

```r
make_synthetic_lm22 <- function(markers_per_type = 3L, baseline = 4, marker = 1024) {
  genes <- paste0("G", seq_len(22L * markers_per_type))
  reference <- matrix(
    baseline,
    nrow = length(genes),
    ncol = 22L,
    dimnames = list(genes, lm22_cell_types())
  )
  for (cell_index in seq_len(22L)) {
    first <- (cell_index - 1L) * markers_per_type + 1L
    last <- cell_index * markers_per_type
    reference[first:last, cell_index] <- marker
  }
  reference
}

reference <- make_synthetic_lm22(markers_per_type = 3L, baseline = 4, marker = 1024)
weights <- matrix(rexp(8 * 22, rate = 1), nrow = 8,
                  dimnames = list(paste0("S", 1:8), lm22_cell_types()))
weights <- weights / rowSums(weights)
bulk_linear <- reference %*% t(weights)
bulk_log <- log2(bulk_linear + 1)

inputs <- prepare_dtangle_inputs(bulk_log, reference, 0.80, FALSE)
fit <- estimate_dtangle(inputs, marker_fraction = 0.10)
testthat::expect_equal(dim(fit$proportions), c(8L, 22L))
testthat::expect_true(all(fit$proportions >= 0))
testthat::expect_equal(rowSums(fit$proportions), rep(1, 8), tolerance = 1e-8)
```

- [ ] **Step 3: Run the tests and verify failure**

Run: `Rscript tests/testthat.R tests/testthat/test-dtangle-stage.R`

Expected: FAIL because `transform_lm22()` and `estimate_dtangle()` do not exist.

- [ ] **Step 4: Implement LM22 validation and alignment**

Implement these rules:

```r
transform_lm22 <- function(lm22_linear) {
  validate_lm22(lm22_linear)
  log2(lm22_linear)
}

common_genes <- intersect(rownames(lm22_linear), rownames(bulk_log))
overlap_fraction <- length(common_genes) / nrow(lm22_linear)
if (overlap_fraction < min_overlap) {
  stop(sprintf("LM22 overlap %.3f is below %.3f.", overlap_fraction, min_overlap))
}
```

Keep the LM22 gene order. When `quantile_normalize` is true, combine `transformed_lm22` and `bulk_log` as a genes-by-profiles matrix and call `limma::normalizeBetweenArrays()` once. Return `Y` as samples by genes and `references` as cell types by genes.

- [ ] **Step 5: Implement dtangle and its reports**

```r
fit <- dtangle::dtangle(
  Y = inputs$Y,
  references = inputs$references,
  n_markers = marker_fraction,
  data_type = "rna-seq",
  marker_method = "ratio"
)

proportions <- fit$estimates
colnames(proportions) <- rownames(inputs$references)
markers <- purrr::map2_dfr(fit$markers, rownames(inputs$references), function(indices, cell_type) {
  tibble::tibble(
    cell_type = cell_type,
    marker_rank = seq_along(indices),
    gene_symbol = colnames(inputs$Y)[indices]
  )
})
```

Set metadata fields for dtangle version, gamma, marker method, marker fraction, overlap count, overlap fraction, quantile-normalization status, sample count, and reference dimensions.

Stop if any cell type has zero selected markers. Validate that every estimate is finite and nonnegative and that each row sum differs from one by no more than `1e-8`.

- [ ] **Step 6: Add the dtangle command-line wrapper**

Use options `--bulk-log`, `--lm22`, `--min-overlap`, `--marker-fraction`, `--quantile-normalize`, and `--output-dir`. Write `dtangle_proportions.tsv`, `dtangle_markers.tsv`, `dtangle_metadata.json`, `dtangle_overlap.tsv`, `dtangle_lm22_log.tsv.gz`, and `dtangle_shared_bulk.tsv.gz`.

- [ ] **Step 7: Run focused tests and lint**

Run: `Rscript tests/testthat.R tests/testthat/test-dtangle-stage.R`

Expected: PASS, with all eight synthetic sample rows summing to one.

Run: `Rscript tools/lint_r.R`

Expected: `R lint ok`.

- [ ] **Step 8: Commit dtangle support**

```bash
git add R/dtangle_stage.R scripts/run_dtangle.R tests/testthat/test-dtangle-stage.R
git commit -m "feat: estimate LM22 proportions with dtangle"
```

---

### Task 4: Major-Group Proportion Processing

**Files:**
- Create: `R/proportions.R`
- Create: `scripts/process_proportions.R`
- Create: `tests/testthat/test-proportions.R`

**Interfaces:**
- Consumes: A sample-by-22-type proportion matrix from dtangle or restart input.
- Produces: `combine_lm22_proportions(proportions) -> numeric matrix`, `filter_and_adjust_groups(combined, mean_threshold, zero_floor) -> list(weights, report)`, `process_proportions(proportions, mean_threshold, zero_floor) -> list(original, combined, tca_weights, report)`, and four TSV outputs.

- [ ] **Step 1: Write failing aggregation and filter tests**

```r
testthat::test_that("all 22 LM22 columns map to the approved nine groups", {
  x <- matrix(seq_len(22), nrow = 1,
              dimnames = list("S1", lm22_cell_types()))
  observed <- combine_lm22_proportions(x)
  expected <- vapply(lm22_group_map(), function(members) sum(x[1, members]), numeric(1))
  testthat::expect_equal(as.numeric(observed[1, ]), unname(expected))
  testthat::expect_identical(colnames(observed), names(lm22_group_map()))
})

testthat::test_that("filter, zero floor, and normalization follow the specification", {
  combined <- matrix(
    c(0.7, 0.29995, 0.00005, 0.6, 0.4, 0), nrow = 2, byrow = TRUE,
    dimnames = list(c("S1", "S2"), c("B cells", "CD4 T cells", "Eosinophils"))
  )
  result <- filter_and_adjust_groups(combined, 0.0001, 1e-6)
  testthat::expect_false("Eosinophils" %in% colnames(result$weights))
  testthat::expect_true(all(result$weights > 0))
  testthat::expect_equal(rowSums(result$weights), c(S1 = 1, S2 = 1))
})
```

Add tests for missing LM22 columns, negative values, nonfinite values, invalid input row sums, and fewer than two retained groups.

Use an input proportion row-sum tolerance of `1e-6` so restart files can contain ordinary decimal rounding.

- [ ] **Step 2: Run the tests and verify failure**

Run: `Rscript tests/testthat.R tests/testthat/test-proportions.R`

Expected: FAIL because `combine_lm22_proportions()` does not exist.

- [ ] **Step 3: Implement grouping, filtering, and adjustment**

Use `vapply()` over `lm22_group_map()` for stable group order. Retain groups with `colMeans(combined) >= mean_threshold`. Replace only exact zeros with `zero_floor`, then normalize with `sweep(adjusted, 1, rowSums(adjusted), "/")`.

The report must contain `cell_group`, `cohort_mean`, `threshold`, `retained`, `filter_reason`, `zero_count_before`, and `zero_floor`.

- [ ] **Step 4: Add the command-line wrapper**

Use options `--proportions`, `--mean-threshold`, `--zero-floor`, and `--output-dir`. Write `proportions_lm22.tsv`, `proportions_combined.tsv`, `proportions_tca_weights.tsv`, and `cell_group_filter_report.tsv`.

- [ ] **Step 5: Run tests and commit**

Run: `Rscript tests/testthat.R tests/testthat/test-proportions.R`

Expected: PASS.

```bash
git add R/proportions.R scripts/process_proportions.R tests/testthat/test-proportions.R
git commit -m "feat: combine and filter LM22 proportions"
```

---

### Task 5: Cohort-Wide TCA Fit and Deterministic Gene Shards

**Files:**
- Create: `R/tca_stage.R`
- Create: `scripts/fit_tca.R`
- Create: `tests/testthat/test-tca-stage.R`

**Interfaces:**
- Consumes: Gene-by-sample `log2(TPM + 1)`, sample-by-group positive weights, optional sample-by-covariate `C2`, CPU count, maximum iterations, shard size, and random seed.
- Produces: `validate_tca_inputs(X, W, C2 = NULL) -> invisible(TRUE)`, `remove_constant_features(X) -> list(matrix, report)`, `make_gene_shard_manifest(genes, shard_size) -> tibble`, `fit_tca_stage(X, W, C2, num_cores, max_iters, random_seed, log_file) -> list(model, X, excluded_genes)`, `write_gene_shards(genes, shard_size, output_dir) -> tibble`, `tca_model.rds`, `tca_expression.tsv.gz`, `tca_excluded_genes.tsv`, `gene_shard_manifest.tsv`, and zero-padded shard files.

- [ ] **Step 1: Write failing TCA input and shard tests**

```r
testthat::test_that("TCA rejects reordered samples", {
  X <- matrix(rnorm(40), nrow = 4,
              dimnames = list(paste0("G", 1:4), paste0("S", 1:10)))
  W <- matrix(0.5, nrow = 10, ncol = 2,
              dimnames = list(rev(colnames(X)), c("A", "B")))
  testthat::expect_error(validate_tca_inputs(X, W), "sample order")
})

testthat::test_that("constant genes are removed and reported", {
  X <- rbind(variable = 1:10, constant = rep(2, 10))
  result <- remove_constant_features(X)
  testthat::expect_identical(rownames(result$matrix), "variable")
  testthat::expect_equal(result$report$gene_symbol, "constant")
})

testthat::test_that("shards preserve deterministic gene order", {
  manifest <- make_gene_shard_manifest(paste0("G", 1:7), 3L)
  testthat::expect_equal(manifest$shard_id, c(1, 1, 1, 2, 2, 2, 3))
  testthat::expect_identical(manifest$gene_symbol, paste0("G", 1:7))
})
```

- [ ] **Step 2: Write the failing small TCA fit test**

```r
testthat::test_that("one model fits all eligible genes without refitting W", {
  set.seed(20260901)
  data <- TCA::test_data(20, 24, 3, 0, 0, 0.01)
  result <- fit_tca_stage(
    X = data$X, W = data$W, C2 = NULL,
    num_cores = 1L, max_iters = 2L, random_seed = 20260901L,
    log_file = tempfile()
  )
  testthat::expect_identical(result$model$W, data$W)
  testthat::expect_equal(dim(result$model$mus_hat), c(20L, 3L))
  testthat::expect_true(is.finite(result$model$tau_hat))
})
```

- [ ] **Step 3: Run the tests and verify failure**

Run: `Rscript tests/testthat.R tests/testthat/test-tca-stage.R`

Expected: FAIL because `validate_tca_inputs()` and `fit_tca_stage()` do not exist.

- [ ] **Step 4: Implement strict TCA validation and optional covariates**

Require identical sample IDs and order across `X`, `W`, and `C2`. Require finite values, positive weights, row sums within `1e-8` of one, at least two retained groups, and no intercept column in `C2`. Reject missing covariate values.

- [ ] **Step 5: Implement the cohort-wide model fit**

```r
set.seed(random_seed)
model <- TCA::tca(
  X = X,
  W = W,
  C2 = C2,
  refit_W = FALSE,
  vars.mle = FALSE,
  constrain_mu = FALSE,
  parallel = num_cores > 1L,
  num_cores = num_cores,
  max_iters = max_iters,
  log_file = log_file,
  verbose = TRUE
)
```

Do not split samples or genes before this fit. Save the model with `saveRDS(model, "tca_model.rds", compress = "xz")`.

- [ ] **Step 6: Implement zero-padded shard lists and the CLI**

Write `shard_00001.txt`, `shard_00002.txt`, and later files in input gene order. The manifest columns must be `gene_index`, `gene_symbol`, `shard_id`, `shard_name`, and `index_within_shard`.

Use CLI options `--expression-log`, `--weights`, `--covariates`, `--num-cores`, `--max-iters`, `--shard-size`, `--random-seed`, and `--output-dir`.

- [ ] **Step 7: Run tests and commit**

Run: `Rscript tests/testthat.R tests/testthat/test-tca-stage.R`

Expected: PASS.

```bash
git add R/tca_stage.R scripts/fit_tca.R tests/testthat/test-tca-stage.R
git commit -m "feat: fit cohort-wide TCA model"
```

---

### Task 6: Tensor Extraction, HDF5 Assembly, and Quality Control

**Files:**
- Create: `R/tensor_outputs.R`
- Create: `R/qc.R`
- Create: `scripts/extract_tca_shard.R`
- Create: `scripts/assemble_tca_outputs.R`
- Create: `scripts/build_manifest.R`
- Create: `tests/testthat/test-tensor-outputs.R`

**Interfaces:**
- Consumes: The complete log-expression matrix, TCA model RDS, one gene list, gene-shard manifest, filtered weights, optional C2, and all shard HDF5 files.
- Produces: `extract_tensor_shard(X, model, genes, num_cores, log_file) -> named list of gene-by-sample matrices`, `write_tensor_shard(path, tensor, shard_id) -> invisible(path)`, `assemble_hdf5_shards(shard_paths, manifest, output_dir, pipeline_version, tca_version) -> named character`, one HDF5 file per shard, one final HDF5 file per cell group, optional TSV.GZ files, per-sample reconstruction metrics, QC plots, and `output_manifest.json`.

- [ ] **Step 1: Write failing tensor and HDF5 contract tests**

```r
testthat::test_that("tensor extraction keeps feature, sample, and source order", {
  set.seed(20260901)
  data <- TCA::test_data(12, 16, 3, 0, 0, 0.01)
  model <- TCA::tca(data$X, data$W, refit_W = FALSE, max_iters = 2, verbose = FALSE)
  genes <- rownames(data$X)[1:5]
  tensor <- extract_tensor_shard(data$X, model, genes, 1L, tempfile())
  testthat::expect_identical(names(tensor), colnames(data$W))
  testthat::expect_true(all(vapply(tensor, identical, logical(1),
                                   data$X[genes, , drop = FALSE]) == FALSE))
  testthat::expect_identical(rownames(tensor[[1]]), genes)
  testthat::expect_identical(colnames(tensor[[1]]), colnames(data$X))
})

testthat::test_that("assembled HDF5 follows the documented schema", {
  samples <- c("S1", "S2")
  shard_1 <- list("B cells" = matrix(
    1:4, nrow = 2, dimnames = list(c("G1", "G2"), samples)
  ))
  shard_2 <- list("B cells" = matrix(
    5:8, nrow = 2, dimnames = list(c("G3", "G4"), samples)
  ))
  paths <- c(tempfile(fileext = ".h5"), tempfile(fileext = ".h5"))
  write_tensor_shard(paths[[1]], shard_1, 1L)
  write_tensor_shard(paths[[2]], shard_2, 2L)
  manifest <- tibble::tibble(
    gene_index = 1:4,
    gene_symbol = paste0("G", 1:4),
    shard_id = c(1L, 1L, 2L, 2L),
    shard_name = c("shard_00001", "shard_00001", "shard_00002", "shard_00002"),
    index_within_shard = c(1L, 2L, 1L, 2L)
  )
  output <- assemble_hdf5_shards(
    paths, manifest, tempfile(), "test", as.character(utils::packageVersion("TCA"))
  )
  h5 <- hdf5r::H5File$new(output[["B cells"]], mode = "r")
  on.exit(h5$close_all())
  testthat::expect_true(all(c("expression", "gene_id", "sample_id") %in% names(h5)))
  testthat::expect_identical(h5$attr_open("scale")$read(), "log2_tpm_plus_1")
})
```

Add tests for a missing shard, duplicate gene, wrong sample order, nonfinite tensor value, and an output dimension mismatch.

- [ ] **Step 2: Run the tests and verify failure**

Run: `Rscript tests/testthat.R tests/testthat/test-tensor-outputs.R`

Expected: FAIL because `extract_tensor_shard()` does not exist.

- [ ] **Step 3: Implement shard extraction**

```r
subset_model <- TCA::tcasub(model, genes, log_file = log_file, verbose = TRUE)
tensor <- TCA::tensor(
  X = X[genes, , drop = FALSE],
  tca.mdl = subset_model,
  scale = FALSE,
  parallel = num_cores > 1L,
  num_cores = num_cores,
  log_file = log_file,
  verbose = TRUE
)
names(tensor) <- colnames(model$W)
```

Create a group named with the exact value from `slugify_cell_group(source)` under `/sources`, then create its `expression` dataset. For example, B cells use `/sources/b_cells/expression` and CD4 T cells use `/sources/cd4_t_cells/expression`. Also create `/gene_id` and `/sample_id`, and add `shard_id`, `source_names`, and `scale` attributes.

- [ ] **Step 4: Implement streaming final HDF5 assembly**

Create one file per retained group with fixed dimensions from the manifest. Use a chunk size of `min(500, n_genes)` by `min(256, n_samples)` and gzip level 6. Fill the final dataset shard by shard at each gene's manifest position. Store:

```text
/expression     numeric, genes x samples
/gene_id        UTF-8 strings
/sample_id      UTF-8 strings
attributes:
  cell_group
  scale = log2_tpm_plus_1
  pipeline_version
  tca_version
```

When `--write-tsv true`, stream rows through a `gzfile()` connection instead of holding another full copy in memory.

- [ ] **Step 5: Implement reconstruction statistics**

For every shard, compute:

```r
source_reconstruction <- Reduce(`+`, purrr::map2(
  tensor,
  seq_len(ncol(weights)),
  function(source_matrix, source_index) {
    sweep(source_matrix, 2, weights[, source_index], "*")
  }
))

covariate_term <- if (is.null(C2) || ncol(C2) == 0L) {
  matrix(0, nrow(source_reconstruction), ncol(source_reconstruction))
} else {
  t(C2 %*% t(subset_model$deltas_hat))
}

reconstructed <- source_reconstruction + covariate_term
```

Accumulate `n`, `sum_x`, `sum_y`, `sum_x2`, `sum_y2`, `sum_xy`, and `sum_squared_error` for each sample. Derive Pearson correlation and RMSE after all shards without rereading the complete assembled matrices.

- [ ] **Step 6: Implement QC plots and the output manifest**

Write `qc_summary.tsv`, `reconstruction_by_sample.tsv`, `qc_plots.pdf`, and `output_manifest.json`. Use `ggplot2::theme_minimal()`, no title, and no subtitle. The manifest must include logical output name, file name, SHA-256, dimensions, scale, cell group, software versions, parameters, and container image input.

- [ ] **Step 7: Add the three command-line wrappers**

Each wrapper must log stage name, UTC start and completion times, dimensions, files written, and errors. Use nonzero exit status for every validation failure.

- [ ] **Step 8: Run tests and commit**

Run: `Rscript tests/testthat.R tests/testthat/test-tensor-outputs.R`

Expected: PASS with finite reconstruction metrics and exact HDF5 order.

```bash
git add R/tensor_outputs.R R/qc.R scripts/extract_tca_shard.R scripts/assemble_tca_outputs.R scripts/build_manifest.R tests/testthat/test-tensor-outputs.R
git commit -m "feat: extract and assemble TCA expression matrices"
```

---

### Task 7: Modular WDL Tasks with Required Logging

**Files:**
- Create: `workflows/tasks/expression.wdl`
- Create: `workflows/tasks/dtangle.wdl`
- Create: `workflows/tasks/proportions.wdl`
- Create: `workflows/tasks/tca.wdl`
- Create: `workflows/tasks/qc.wdl`
- Create: `tests/testthat/test-wdl-contract.R`
- Create: `tools/check_wdl_logging.py`

**Interfaces:**
- Consumes: The R command-line interfaces from Tasks 2 through 6 and a task-level `docker_image` input.
- Produces: WDL 1.1 tasks `PrepareExpression`, `RunDtangle`, `ProcessProportions`, `FitTca`, `ExtractTcaShard`, `AssembleTca`, and `BuildManifest`.

- [ ] **Step 1: Write failing WDL contract and logging tests**

```r
testthat::test_that("all pipeline WDL tasks exist", {
  expected <- c(
    "PrepareExpression", "RunDtangle", "ProcessProportions",
    "FitTca", "ExtractTcaShard", "AssembleTca", "BuildManifest"
  )
  files <- list.files("workflows/tasks", pattern = "[.]wdl$", full.names = TRUE)
  text <- files |>
    purrr::map_chr(~ paste(readLines(.x), collapse = "\n")) |>
    paste(collapse = "\n")
  purrr::walk(expected, ~ testthat::expect_match(text, paste0("task ", .x)))
})
```

Implement `tools/check_wdl_logging.py` so it extracts every `command <<< ... >>>` block and fails unless that block contains `stage=`, `start_time=`, and `completion_time=`.

- [ ] **Step 2: Run the tests and verify failure**

Run: `Rscript tests/testthat.R tests/testthat/test-wdl-contract.R`

Expected: FAIL because the WDL task files do not exist.

- [ ] **Step 3: Implement the WDL task interfaces**

Use these task inputs and outputs:

```text
PrepareExpression:
  in  expression, expression_type, annotation?, docker_image, cpu, memory, disk_gb
  out prepared_tpm, prepared_log, mapping_report, excluded_genes, log
RunDtangle:
  in  prepared_log, lm22, min_overlap, marker_fraction, quantile_normalize,
      docker_image, cpu, memory, disk_gb
  out proportions, markers, metadata, overlap_report, transformed_lm22,
      shared_bulk, log
ProcessProportions:
  in  proportions, mean_threshold, zero_floor, docker_image, cpu, memory, disk_gb
  out original, combined, tca_weights, filter_report, log
FitTca:
  in  prepared_log, tca_weights, covariates?, shard_size, max_iters, random_seed,
      docker_image, cpu, memory, disk_gb
  out model, tca_expression, excluded_genes, shard_manifest, Array[File] shards, log
ExtractTcaShard:
  in  tca_expression, model, shard, docker_image, cpu, memory, disk_gb
  out shard_hdf5, log
AssembleTca:
  in  Array[File] shard_hdf5, shard_manifest, tca_expression, model,
      tca_weights, covariates?, write_tsv, pipeline_version,
      docker_image, cpu, memory, disk_gb
  out Array[File] group_hdf5, Array[File] group_tsv,
      reconstruction_by_sample, assembly_qc, log
BuildManifest:
  in  all primary/supporting files and parameter values
  out output_manifest, qc_summary, qc_plots, provenance, log
```

- [ ] **Step 4: Add the logging template to every command block**

Use this pattern with a task-specific `stage` value:

```wdl
command <<<
  set -euo pipefail
  stage="prepare_expression"
  log="$stage.log"
  printf 'stage=%s start_time=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  trap 'status=$?; printf "stage=%s error_status=%s time=%s\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
  Rscript /opt/celltype/scripts/prepare_expression.R \
    --expression '~{expression}' \
    --expression-type '~{expression_type}' \
    --annotation '~{if defined(annotation) then annotation else ""}' \
    --output-dir outputs 2>&1 | tee -a "$log"
  printf 'stage=%s completion_time=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
>>>
```

Use `runtime { docker: docker_image }` and explicit `cpu`, `memory`, and `disks` settings in every task.

Every task also accepts `preemptible_attempts` and `max_retries` integer inputs and passes them to the Terra/Cromwell `preemptible` and `maxRetries` runtime keys. This makes every value in the following table overridable without a WDL source edit.

Use these overridable defaults for the 9,000-sample starting profile:

| Task | CPU | Memory | Disk | Preemptible attempts | Retries |
| --- | ---: | ---: | ---: | ---: | ---: |
| PrepareExpression | 4 | 64 GB | 400 GB | 2 | 2 |
| RunDtangle | 4 | 32 GB | 100 GB | 2 | 2 |
| ProcessProportions | 2 | 16 GB | 50 GB | 2 | 2 |
| FitTca | 16 | 192 GB | 750 GB | 0 | 1 |
| ExtractTcaShard | 8 | 64 GB | 200 GB | 2 | 2 |
| AssembleTca | 8 | 128 GB | 500 GB | 0 | 1 |
| BuildManifest | 4 | 32 GB | 100 GB | 1 | 2 |

- [ ] **Step 5: Validate task syntax and logging**

Run: `miniwdl check workflows/tasks/expression.wdl workflows/tasks/dtangle.wdl workflows/tasks/proportions.wdl workflows/tasks/tca.wdl workflows/tasks/qc.wdl`

Expected: No syntax or type errors.

Run: `python3 tools/check_wdl_logging.py workflows`

Expected: `All WDL command blocks contain required logging fields.`

- [ ] **Step 6: Run tests and commit**

Run: `Rscript tests/testthat.R tests/testthat/test-wdl-contract.R`

Expected: PASS.

```bash
git add workflows/tasks tests/testthat/test-wdl-contract.R tools/check_wdl_logging.py
git commit -m "feat: add logged WDL pipeline tasks"
```

---

### Task 8: Top-Level Workflow and Deterministic Smoke Fixtures

**Files:**
- Create: `workflows/cell_type_deconvolution.wdl`
- Create: `scripts/generate_synthetic_fixture.R`
- Create: `tests/smoke/assert_outputs.R`
- Create: generated files under `tests/fixtures/`
- Create: `tests/fixtures/dtangle.inputs.json`
- Create: `tests/fixtures/restart.inputs.json`

**Interfaces:**
- Consumes: Counts or TPM, optional annotation, optional LM22, optional precomputed 22-type proportions, optional C2 covariates, scientific thresholds, and run-time settings.
- Produces: All primary and supporting outputs in the specification through one workflow namespace.

- [ ] **Step 1: Add a failing top-level workflow contract test**

```r
testthat::test_that("top-level workflow exposes both proportion modes", {
  text <- paste(readLines("workflows/cell_type_deconvolution.wdl"), collapse = "\n")
  testthat::expect_match(text, "File[?] lm22")
  testthat::expect_match(text, "File[?] precomputed_proportions")
  testthat::expect_match(text, "if .*precomputed_proportions")
  testthat::expect_match(text, "scatter .*shard")
})
```

- [ ] **Step 2: Run the test and verify failure**

Run: `Rscript tests/testthat.R tests/testthat/test-wdl-contract.R`

Expected: FAIL because `workflows/cell_type_deconvolution.wdl` does not exist.

- [ ] **Step 3: Implement the top-level workflow**

Import all task WDLs. Always call `PrepareExpression`. Use this branch pattern:

```wdl
Boolean estimate_proportions = !defined(precomputed_proportions)

if (estimate_proportions) {
  call dtangle_tasks.RunDtangle {
    input:
      prepared_log = PrepareExpression.prepared_log,
      lm22 = select_first([lm22]),
      min_overlap = min_lm22_overlap,
      marker_fraction = dtangle_marker_fraction,
      quantile_normalize = dtangle_quantile_normalize,
      docker_image = docker_image
  }
}

File proportions_for_processing = select_first([
  precomputed_proportions,
  RunDtangle.proportions
])
```

Then call `ProcessProportions`, `FitTca`, scatter `ExtractTcaShard` over `FitTca.shards`, call `AssembleTca`, and call `BuildManifest`. Expose the three proportion tables, HDF5 array, optional TSV array, QC, model, logs, and manifest.

- [ ] **Step 4: Create a deterministic synthetic fixture generator**

Use seed `20260901`. Generate:

- A 66-gene positive linear signature with all 22 LM22 column names and three strong markers per type.
- Twelve bulk TPM samples generated by linear mixtures of the synthetic references plus small positive noise.
- A count matrix and annotation that reproduce the same sample and gene identifiers.
- One technical covariate named `batch_indicator`.
- A precomputed 22-type proportion table.
- Expected sample, gene, and group files for assertions.

The fixture must not copy, approximate, or derive values from LM22.

Run: `Rscript scripts/generate_synthetic_fixture.R tests/fixtures`

Expected: Deterministic TSV files and both JSON input files.

- [ ] **Step 5: Implement smoke-output assertions**

`tests/smoke/assert_outputs.R` must check:

```r
stopifnot(all(is.finite(proportions)))
stopifnot(all(proportions >= 0))
stopifnot(max(abs(rowSums(proportions) - 1)) < 1e-8)
stopifnot(all(tca_weights > 0))
stopifnot(max(abs(rowSums(tca_weights) - 1)) < 1e-8)
stopifnot(all(expected_samples == observed_samples))
stopifnot(all(expected_genes == observed_genes))
stopifnot(all(is.finite(reconstruction$correlation)))
stopifnot(all(is.finite(reconstruction$rmse)))
```

Also open every final HDF5 and check dimensions, identifiers, `cell_group`, and `scale = log2_tpm_plus_1`.

- [ ] **Step 6: Validate the complete WDL**

Run: `miniwdl check workflows/cell_type_deconvolution.wdl`

Expected: No syntax, import, or type errors.

- [ ] **Step 7: Commit the workflow and fixtures**

```bash
git add workflows/cell_type_deconvolution.wdl scripts/generate_synthetic_fixture.R tests/fixtures tests/smoke tests/testthat/test-wdl-contract.R
git commit -m "feat: add complete deconvolution workflow"
```

---

### Task 9: Pinned Micromamba Image and GitHub Actions Smoke Test

**Files:**
- Create: `envs/environment.yml`
- Modify: `envs/Dockerfile`
- Create: `.github/workflows/pipeline-ci.yml`
- Modify: `.github/workflows/docker-image.yml`
- Delete: `.github/workflows/r-lint.yml`
- Delete: `.github/workflows/wdl-validation.yml`
- Modify: `README.md` badge block by running `tools/update_workflow_badges.py`

**Interfaces:**
- Consumes: Repository source and the synthetic input JSON files.
- Produces: A pinned test image, passing R/WDL/smoke checks on pull requests, and versioned GHCR images on `main`.

- [ ] **Step 1: Write the exact environment lock input**

```yaml
name: celltype-deconvolution
channels:
  - conda-forge
  - bioconda
  - nodefaults
dependencies:
  - r-base=4.5.3
  - r-dtangle=2.0.10
  - r-tca=1.2.1
  - r-tidyverse=2.0.0
  - r-optparse=1.8.2
  - r-jsonlite=2.0.0
  - r-digest=0.6.39
  - r-hdf5r=1.3.12
  - r-testthat=3.3.2
  - r-lintr=3.4.0
  - bioconductor-limma=3.66.0
```

- [ ] **Step 2: Implement the micromamba Dockerfile**

```dockerfile
FROM mambaorg/micromamba:2.9.0-ubuntu22.04

COPY --chown=$MAMBA_USER:$MAMBA_USER envs/environment.yml /tmp/environment.yml
RUN micromamba install --yes --name base --file /tmp/environment.yml && \
    micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER R /opt/celltype/R
COPY --chown=$MAMBA_USER:$MAMBA_USER scripts /opt/celltype/scripts
ENV R_LIBS_USER=/opt/conda/lib/R/library

RUN Rscript -e 'stopifnot(as.character(packageVersion("dtangle")) == "2.0.10")' && \
    Rscript -e 'stopifnot(as.character(packageVersion("TCA")) == "1.2.1")'

WORKDIR /workspace
```

- [ ] **Step 3: Add the pull-request CI workflow**

`pipeline-ci.yml` must:

1. Check out the repository.
2. Install `miniwdl==1.15.0` with pipx.
3. Run `miniwdl check workflows/cell_type_deconvolution.wdl`.
4. Run `python3 tools/check_wdl_logging.py workflows`.
5. Build `celltype-deconvolution:test` from `envs/Dockerfile`.
6. Run `Rscript tools/lint_r.R` and `Rscript tests/testthat.R` inside that image.
7. Run the complete dtangle-mode smoke workflow with miniwdl and the local image tag.
8. Run `tests/smoke/assert_outputs.R` on its outputs.
9. Run the restart-mode smoke workflow and its assertions.
10. Upload miniwdl logs and QC files when a step fails.

Use `docker_image = "celltype-deconvolution:test"` in both smoke JSON files.

- [ ] **Step 4: Update the main-branch image workflow**

Keep `docker/login-action` and `docker/build-push-action`. Set `push: true` on `main`. Publish a `sha-` tag that contains the commit SHA and publish the `latest` tag. Add OCI source, revision, and created labels. Do not run a separate local smoke test in this workflow; `pipeline-ci.yml` is the required gate.

- [ ] **Step 5: Consolidate redundant validation workflows and update badges**

Delete the old separate R-lint and WDL-validation workflows after their checks exist in `pipeline-ci.yml`.

Run: `python3 tools/update_workflow_badges.py`

Expected: README badges list Docker Image CI, Pipeline CI, and Update README workflow badges.

- [ ] **Step 6: Perform static verification without a local image build**

Run: `miniwdl check workflows/cell_type_deconvolution.wdl`

Run: `python3 tools/check_wdl_logging.py workflows`

Run: `git diff --check`

Expected: All commands pass. Do not run `docker build` locally.

- [ ] **Step 7: Commit the environment and CI**

```bash
git add envs .github/workflows README.md
git commit -m "ci: build and smoke test the deconvolution image"
```

---

### Task 10: Terra, Dockstore, Data Dictionary, and User Documentation

**Files:**
- Modify: `README.md`
- Create: `.dockstore.yml`
- Create: `docs/terra.md`
- Create: `docs/data-dictionary.md`
- Create: `examples/counts.inputs.json`
- Create: `examples/tpm.inputs.json`
- Create: `examples/precomputed-proportions.inputs.json`
- Create: `tests/testthat/test-documentation.R`

**Interfaces:**
- Consumes: Final WDL input and output names.
- Produces: Importable Dockstore metadata, Terra-ready examples, and exact input/output documentation.

- [ ] **Step 1: Write failing documentation contract tests**

```r
testthat::test_that("all example JSON files parse and use workflow input names", {
  files <- list.files("examples", pattern = "[.]inputs[.]json$", full.names = TRUE)
  testthat::expect_length(files, 3L)
  parsed <- purrr::map(files, jsonlite::read_json)
  testthat::expect_true(all(purrr::map_lgl(parsed, ~ "cell_type_deconvolution.expression" %in% names(.x))))
})

testthat::test_that("LM22 is documented as user supplied and linear", {
  text <- paste(readLines("docs/terra.md"), collapse = "\n")
  testthat::expect_match(text, "user-supplied")
  testthat::expect_match(text, "linear LM22")
  testthat::expect_match(text, "log2\\(LM22\\)")
})
```

- [ ] **Step 2: Run the documentation tests and verify failure**

Run: `Rscript tests/testthat.R tests/testthat/test-documentation.R`

Expected: FAIL because the three examples and Terra guide do not exist.

- [ ] **Step 3: Write Dockstore metadata**

```yaml
version: 1.2
workflows:
  - subclass: WDL
    primaryDescriptorPath: /workflows/cell_type_deconvolution.wdl
    testParameterFiles:
      - /tests/fixtures/dtangle.inputs.json
    name: cell-type-deconvolution
```

- [ ] **Step 4: Replace the template README**

Document the method, accepted matrices, required LM22 file, two proportion modes, major-group mapping, output scale, quick Terra import, local static checks, GitHub Actions smoke test, and license limitations. State that the repository does not contain LM22.

- [ ] **Step 5: Write the Terra run guide**

Include:

- How to upload or reference counts, TPM, annotation, LM22, covariates, and precomputed proportions.
- Which inputs are required in each mode.
- A 9,000-sample starting profile for CPU, memory, disk, retry, and preemptible settings.
- How to resume with precomputed proportions.
- How to read the manifest and HDF5 outputs.
- How to interpret `log2_tpm_plus_1` and the limitations of cross-platform LM22 use.

- [ ] **Step 6: Write the data dictionary and examples**

List every input and output field with type, required mode, default, units or scale, validation rule, and description. Ensure JSON keys exactly match the top-level WDL namespace.

- [ ] **Step 7: Run documentation and full static checks**

Run: `Rscript tests/testthat.R tests/testthat/test-documentation.R`

Expected: PASS.

Run: `miniwdl check workflows/cell_type_deconvolution.wdl`

Expected: PASS.

Run: `git diff --check`

Expected: No output.

- [ ] **Step 8: Commit documentation**

```bash
git add README.md .dockstore.yml docs/terra.md docs/data-dictionary.md examples tests/testthat/test-documentation.R
git commit -m "docs: add Terra and Dockstore usage guides"
```

---

### Task 11: Final Verification and GitHub Actions Evidence

**Files:**
- Modify only files required by verified failures.

**Interfaces:**
- Consumes: The complete implementation branch.
- Produces: A clean branch and passing GitHub Actions evidence for R tests, WDL validation, both smoke modes, and image construction.

- [ ] **Step 1: Run all non-container local checks**

Run: `Rscript tests/testthat.R`

Expected: All tests pass when the pinned R environment is active.

Run: `Rscript tools/lint_r.R`

Expected: `R lint ok`.

Run: `miniwdl check workflows/cell_type_deconvolution.wdl`

Expected: No errors.

Run: `python3 tools/check_wdl_logging.py workflows`

Expected: All command blocks pass.

Run: `git diff --check`

Expected: No output.

- [ ] **Step 2: Confirm no licensed or secret material is tracked**

Run: `git ls-files | rg -i 'lm22[.]txt|token|credential|secret'`

Expected: No LM22 matrix, token, credential, or secret file. Documentation references to LM22 are allowed only when their path does not match `LM22.txt`.

- [ ] **Step 3: Push the implementation branch and wait for CI**

Run: `git push -u origin feat/dtangle-tca-pipeline`

Run: `gh run list --workflow pipeline-ci.yml --branch feat/dtangle-tca-pipeline --limit 1`

Expected: One completed successful run.

- [ ] **Step 4: Inspect the successful workflow checks**

Run: `gh run view "$(gh run list --workflow pipeline-ci.yml --branch feat/dtangle-tca-pipeline --limit 1 --json databaseId --jq '.[0].databaseId')" --log-failed`

Expected: No failed logs. The run summary must show the image build, R tests, R lint, WDL validation, dtangle-mode smoke test, restart-mode smoke test, and output assertions.

- [ ] **Step 5: Record verification evidence**

Add the successful Actions run URL and tested image identifier to the pull-request description or final handoff. Do not claim the smoke workflow passed until the GitHub Actions run is complete and successful.

- [ ] **Step 6: Commit any verification-only corrections**

If verification required a correction, rerun the failed and full checks, then commit only the corrected files:

```bash
git add R scripts workflows tests envs .github README.md docs examples .dockstore.yml tools
git commit -m "fix: resolve pipeline verification failures"
```

If no correction was needed, do not create an empty commit.
