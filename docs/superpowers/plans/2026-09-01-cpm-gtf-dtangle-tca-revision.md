# CPM/GTF dtangle and TCA Pipeline Revision Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a tested Terra WDL workflow that accepts positive normalized human whole-blood CPM plus a GTF, estimates LM22 proportions with dtangle, combines and filters major lineages, and uses TCA to write one cell-type-specific gene-by-sample matrix for each retained lineage.

**Architecture:** Focused R modules prepare CPM, parse the GTF, run dtangle, process proportions, fit one cohort-wide TCA model, extract tensors by gene shard, and assemble final HDF5 outputs. Separate WDL tasks call these modules. The dtangle path uses only the LM22 gene intersection in LM22 order, while TCA uses the complete mapped expression matrix after constant-gene removal.

**Tech Stack:** WDL 1.1, miniwdl 1.15.0, R 4.5.3, tidyverse 2.0.0, dtangle 2.0.10, TCA 1.2.1, limma 3.66.0, hdf5r 1.3.12, testthat 3.3.1, micromamba 2.9.0, GitHub Actions, GHCR, Terra, and Dockstore.

**Spec:** `docs/superpowers/specs/2026-09-01-whole-blood-dtangle-tca-design.md`

## Global Constraints

- The target data are human whole-blood RNA-seq from approximately 9,000 samples.
- Require a gene-by-sample linear CPM matrix whose first column is `gene_id`.
- Require a GTF or GTF.GZ with `gene` records and `gene_id` and `gene_name` attributes.
- Do not filter by `gene_type` or `gene_biotype`.
- Do not calculate CPM, renormalize CPM, or derive gene lengths.
- Sum CPM rows when multiple input gene IDs map to the same gene name.
- Require every prepared CPM value to be finite and strictly positive. Do not add a pseudocount.
- Use `log2(CPM)` for dtangle mixtures and TCA.
- Accept the standard positive linear LM22 matrix at run time and do not redistribute LM22.
- Use `log2(LM22)` without a pseudocount.
- Keep joint quantile normalization off by default. When enabled, apply `limma::normalizeBetweenArrays()` once to the joined dtangle-only reference and mixture profiles.
- Preserve LM22 gene order in both aligned dtangle matrices and verify identical gene columns.
- Run dtangle 2.0.10 with all 22 LM22 types, `data_type = "rna-seq"`, `marker_method = "ratio"`, and marker fraction `0.10`.
- Combine the 22 proportions into ten groups. Keep gamma-delta T cells and mast cells separate.
- Remove every group whose cohort mean is less than `0.0001`. No group is mandatory.
- Replace exact zeros in retained groups with `1e-6`, then normalize each sample row to one.
- Fit one cohort-wide TCA model with `refit_W = FALSE` on the complete mapped `log2(CPM)` matrix after constant-gene removal.
- Label TCA outputs as `log2_cpm`, not linear CPM.
- Use a default tensor shard size of 500 genes and stable gene and sample order.
- Use random seed `20260901` for stochastic R operations.
- Use tidyverse syntax in R where it improves clarity.
- Use clean minimal plots without titles or subtitles.
- Log stage, start time, completion time, dimensions, output paths, and error context in every WDL command.
- Use a micromamba base image and exact package pins from conda-forge and bioconda.
- Build and smoke-test the image in GitHub Actions. Do not require a local Docker build.
- Preserve unrelated user changes and commit only the files for the active task.

## Starting State

- Worktree: `.worktrees/dtangle-tca-pipeline`
- Branch: `feat/dtangle-tca-pipeline`
- Approved design commit: `86143a5`
- Core I/O, the old counts/TPM expression stage, and the initial dtangle stage already exist.
- `R/dtangle_stage.R` and `tests/testthat/test-dtangle-stage.R` contain uncommitted fixes for dtangle version validation and quantile-normalization column slicing. Keep those fixes.

## Planned File Structure

```text
R/
  constants.R                 LM22 names, ten-group map, and defaults
  io.R                        strict matrix and table I/O
  expression.R                GTF parsing, CPM mapping, collapse, and log2 conversion
  dtangle_stage.R             LM22 validation, alignment, optional QN, and dtangle
  proportions.R               grouping, cohort filter, zero floor, and normalization
  tca_stage.R                 TCA validation, fit, and deterministic gene shards
  tensor_outputs.R            tensor extraction and HDF5 assembly
  qc.R                        reconstruction metrics, plots, and provenance manifest
scripts/
  bootstrap.R
  prepare_expression.R
  run_dtangle.R
  process_proportions.R
  fit_tca.R
  extract_tca_shard.R
  assemble_tca_outputs.R
  build_manifest.R
  generate_synthetic_fixture.R
tests/
  testthat/test-expression.R
  testthat/test-dtangle-stage.R
  testthat/test-proportions.R
  testthat/test-tca-stage.R
  testthat/test-tensor-outputs.R
  testthat/test-wdl-contract.R
  testthat/test-documentation.R
  fixtures/
  smoke/assert_outputs.R
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
  cpm.inputs.json
  precomputed-proportions.inputs.json
docs/
  terra.md
  data-dictionary.md
.github/workflows/
  pipeline-ci.yml
  docker-image.yml
```

---

### Task 1: Replace Counts/TPM Preparation with CPM and GTF Preparation

**Files:**
- Modify: `R/expression.R`
- Modify: `scripts/prepare_expression.R`
- Modify: `tests/testthat/test-expression.R`

**Interfaces:**
- Consumes: `cpm` as a numeric gene-by-sample matrix and `gtf_path` as a GTF or GTF.GZ path.
- Produces: `extract_gtf_attribute(attributes, key) -> character`, `validate_gtf_gene_annotation(annotation) -> tibble`, `read_gtf_gene_annotation(path, chunk_size = 100000L) -> tibble(gene_id, gene_name, gene_type)`, `collapse_cpm_to_gene_names(cpm, annotation) -> numeric matrix`, `prepare_expression(cpm, annotation) -> list(cpm, log_expression, mapping_report, excluded_genes)`, and four CLI files.

- [ ] **Step 1: Replace the old tests with failing CPM/GTF tests**

```r
testthat::test_that("GTF parsing retains all gene types and ignores non-gene records", {
  gtf <- tempfile(fileext = ".gtf")
  writeLines(c(
    "1\tsrc\tgene\t1\t10\t.\t+\t.\tgene_id \"g1\"; gene_name \"A\"; gene_type \"protein_coding\";",
    "1\tsrc\tgene\t20\t30\t.\t+\t.\tgene_id \"g2\"; gene_name \"B\"; gene_type \"lncRNA\";",
    "1\tsrc\ttranscript\t1\t10\t.\t+\t.\tgene_id \"g1\"; gene_name \"A\";"
  ), gtf)
  observed <- read_gtf_gene_annotation(gtf, chunk_size = 1L)
  testthat::expect_identical(observed$gene_id, c("g1", "g2"))
  testthat::expect_identical(observed$gene_name, c("A", "B"))
  testthat::expect_identical(observed$gene_type, c("protein_coding", "lncRNA"))
})

testthat::test_that("duplicate gene names are summed without CPM renormalization", {
  cpm <- matrix(c(2, 3, 5, 7), nrow = 2,
                dimnames = list(c("g1", "g2"), c("s1", "s2")))
  annotation <- tibble::tibble(
    gene_id = c("g1", "g2"), gene_name = c("A", "A"), gene_type = c("x", "y")
  )
  result <- prepare_expression(cpm, annotation)
  testthat::expect_identical(rownames(result$cpm), "A")
  testthat::expect_equal(unname(result$cpm), matrix(c(5, 12), nrow = 1))
  testthat::expect_equal(result$log_expression, log2(result$cpm))
})

testthat::test_that("prepared CPM must be strictly positive", {
  cpm <- matrix(c(1, 0), nrow = 1,
                dimnames = list("g1", c("s1", "s2")))
  annotation <- tibble::tibble(gene_id = "g1", gene_name = "A", gene_type = NA_character_)
  testthat::expect_error(prepare_expression(cpm, annotation), "strictly positive")
})
```

The same test file must verify compressed input, trimmed matching, missing names, duplicate GTF IDs, and the CLI contract:

```r
testthat::test_that("compressed GTF and trimmed IDs are supported", {
  gtf <- tempfile(fileext = ".gtf.gz")
  connection <- gzfile(gtf, "wt")
  writeLines(
    "1\tsrc\tgene\t1\t10\t.\t+\t.\tgene_id \"g1\"; gene_name \"A\"; gene_type \"lncRNA\";",
    connection
  )
  close(connection)
  cpm <- matrix(c(2, 4), nrow = 1,
                dimnames = list(" g1 ", c("s1", "s2")))
  result <- prepare_expression(cpm, read_gtf_gene_annotation(gtf))
  testthat::expect_identical(rownames(result$cpm), "A")
  testthat::expect_equal(as.numeric(result$cpm), c(2, 4))
})

testthat::test_that("missing gene names are reported and duplicate GTF IDs stop", {
  cpm <- matrix(1:4, nrow = 2,
                dimnames = list(c("g1", "g2"), c("s1", "s2")))
  annotation <- tibble::tibble(
    gene_id = c("g1", "g2"), gene_name = c("A", NA_character_),
    gene_type = c("protein_coding", "lncRNA")
  )
  result <- prepare_expression(cpm, annotation)
  testthat::expect_true(any(
    result$mapping_report$gene_id == "g2" &
      result$mapping_report$mapping_action == "missing_gene_name"
  ))
  testthat::expect_error(
    validate_gtf_gene_annotation(dplyr::bind_rows(annotation, annotation[1, ])),
    "gene_id.*unique"
  )
})

testthat::test_that("malformed GTF records stop parsing", {
  gtf <- tempfile(fileext = ".gtf")
  writeLines("1\tsrc\tgene\t1\t10", gtf)
  testthat::expect_error(read_gtf_gene_annotation(gtf), "nine.*fields")
})

testthat::test_that("the CPM CLI writes its four declared outputs", {
  input <- tempfile(fileext = ".tsv")
  gtf <- tempfile(fileext = ".gtf")
  output_dir <- tempfile()
  dir.create(output_dir)
  write_numeric_matrix(
    matrix(c(2, 4), nrow = 1, dimnames = list("g1", c("s1", "s2"))),
    input,
    "gene_id"
  )
  writeLines(
    "1\tsrc\tgene\t1\t10\t.\t+\t.\tgene_id \"g1\"; gene_name \"A\";",
    gtf
  )
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c("scripts/prepare_expression.R", "--expression", shQuote(input),
      "--gtf", shQuote(gtf), "--output-dir", shQuote(output_dir))
  )
  expected <- c(
    "prepared_cpm.tsv.gz", "prepared_log2_cpm.tsv.gz",
    "gene_mapping_report.tsv", "excluded_genes.tsv"
  )
  testthat::expect_equal(status, 0L)
  testthat::expect_true(all(file.exists(file.path(output_dir, expected))))
})
```

- [ ] **Step 2: Run the focused test and confirm the old interface fails**

Run: `Rscript tests/testthat.R tests/testthat/test-expression.R`

Expected: FAIL because `read_gtf_gene_annotation()` does not exist and the current code still requires `expression_type`.

- [ ] **Step 3: Implement chunked GTF gene parsing**

Use a text connection for GTF and `gzfile(path, "rt")` for GTF.GZ. Read at most `chunk_size` lines per iteration. Ignore comments and non-gene records. Parse only the required attributes plus an optional type attribute.

```r
extract_gtf_attribute <- function(attributes, key) {
  pattern <- paste0("(?:^|;)[[:space:]]*", key, "[[:space:]]+\\\"([^\\\"]+)\\\"")
  stringr::str_match(attributes, pattern)[, 2L]
}

gene_rows <- tibble::tibble(
  gene_id = extract_gtf_attribute(fields[, 9L], "gene_id"),
  gene_name = extract_gtf_attribute(fields[, 9L], "gene_name"),
  gene_type = dplyr::coalesce(
    extract_gtf_attribute(fields[, 9L], "gene_type"),
    extract_gtf_attribute(fields[, 9L], "gene_biotype")
  )
)
```

Require nine tab-separated GTF fields for each non-comment record. Trim `gene_id` and `gene_name`. Stop on a malformed record, missing `gene_id`, duplicate `gene_id`, zero gene records, or zero usable `gene_name` values. Do not filter on `gene_type`.

- [ ] **Step 4: Implement CPM mapping and log conversion**

Trim expression `gene_id` row names and reject duplicates created by trimming. Match them to trimmed GTF IDs by exact text. Use the expression row order to form the mapping report. Collapse duplicated gene names with tidyverse operations and preserve sample order.

```r
collapsed <- tibble::as_tibble(cpm, rownames = "gene_id") |>
  tidyr::pivot_longer(-"gene_id", names_to = "sample_id", values_to = "cpm") |>
  dplyr::inner_join(mapped_annotation, by = "gene_id") |>
  dplyr::group_by(.data$gene_name, .data$sample_id) |>
  dplyr::summarise(cpm = sum(.data$cpm), .groups = "drop") |>
  tidyr::pivot_wider(names_from = "sample_id", values_from = "cpm") |>
  dplyr::select("gene_name", dplyr::all_of(colnames(cpm)))
```

Do not divide by column sums. Stop if the collapsed matrix has a missing, nonfinite, zero, or negative value. Return `log_expression = log2(prepared_cpm)`. Record `mapped`, `duplicate_gene_name_collapsed`, `missing_gtf_gene_id`, and `missing_gene_name` actions.

- [ ] **Step 5: Replace the expression CLI contract**

Use only `--expression`, `--gtf`, and `--output-dir` as required options.

```r
result <- prepare_expression(
  cpm = read_numeric_matrix(options$expression, "gene_id"),
  annotation = read_gtf_gene_annotation(options$gtf)
)
write_numeric_matrix(result$cpm, file.path(options$output_dir, "prepared_cpm.tsv.gz"), "gene_name")
write_numeric_matrix(
  result$log_expression,
  file.path(options$output_dir, "prepared_log2_cpm.tsv.gz"),
  "gene_name"
)
```

Log the input dimensions, GTF gene count, mapped gene-name count, collapsed count, output paths, and UTC completion time. Wrap the command body with the same nonzero error handler used by `run_dtangle.R`.

- [ ] **Step 6: Run tests, lint, and commit**

Run: `Rscript tests/testthat.R tests/testthat/test-expression.R`

Expected: PASS.

Run: `Rscript tools/lint_r.R`

Expected: `R lint ok`.

```bash
git add R/expression.R scripts/prepare_expression.R tests/testthat/test-expression.R
git commit -m "feat: prepare positive CPM with GTF annotation"
```

---

### Task 2: Complete the dtangle Stage for log2(CPM)

**Files:**
- Modify: `R/dtangle_stage.R`
- Modify: `scripts/run_dtangle.R`
- Modify: `tests/testthat/test-dtangle-stage.R`

**Interfaces:**
- Consumes: complete mapped gene-by-sample `log2(CPM)`, positive linear LM22, minimum overlap, marker fraction, and the optional QN flag.
- Produces: `validate_lm22(lm22_linear) -> invisible(TRUE)`, `transform_lm22(lm22_linear) -> numeric matrix`, `prepare_dtangle_inputs(bulk_log, lm22_linear, min_overlap, quantile_normalize) -> list(Y, references, transformed_lm22, shared_bulk, overlap_report)`, and `estimate_dtangle(inputs, marker_fraction, marker_method = "ratio") -> list(proportions, markers, metadata)`. `Y` and `references` have identical LM22-ordered genes. `overlap_report` has one row per LM22 gene with `gene_symbol`, `reference_index`, and `matched`.

- [ ] **Step 1: Add failing scale and order tests while preserving current fixes**

```r
testthat::test_that("finite negative log2 CPM values are accepted", {
  lm22 <- make_synthetic_lm22()
  bulk_log <- matrix(
    seq(-3, 3, length.out = nrow(lm22)), ncol = 1,
    dimnames = list(rev(rownames(lm22)), "S1")
  )
  inputs <- prepare_dtangle_inputs(bulk_log, lm22, 0.80, FALSE)
  testthat::expect_identical(colnames(inputs$Y), rownames(lm22))
  testthat::expect_identical(colnames(inputs$references), rownames(lm22))
})

testthat::test_that("synthetic CPM mixtures use log2 without a pseudocount", {
  reference <- make_synthetic_lm22()
  weights <- matrix(rexp(8 * 22), nrow = 8,
                    dimnames = list(paste0("S", 1:8), lm22_cell_types()))
  weights <- weights / rowSums(weights)
  bulk_cpm <- reference %*% t(weights)
  inputs <- prepare_dtangle_inputs(log2(bulk_cpm), reference, 0.80, FALSE)
  testthat::expect_equal(inputs$shared_bulk, log2(bulk_cpm))
})

testthat::test_that("overlap report lists every LM22 gene in reference order", {
  lm22 <- make_synthetic_lm22()
  kept <- rownames(lm22)[-c(2L, 5L)]
  bulk_log <- matrix(1, nrow = length(kept), ncol = 1,
                     dimnames = list(rev(kept), "S1"))
  inputs <- prepare_dtangle_inputs(bulk_log, lm22, 0.80, FALSE)
  testthat::expect_identical(inputs$overlap_report$gene_symbol, rownames(lm22))
  testthat::expect_identical(
    inputs$overlap_report$matched,
    rownames(lm22) %in% kept
  )
})
```

Keep the uncommitted tests for dtangle version `2.0.10` and QN profile slicing by column position.

- [ ] **Step 2: Run the focused test and confirm the negative-value failure**

Run: `Rscript tests/testthat.R tests/testthat/test-dtangle-stage.R`

Expected: FAIL because `validate_bulk_log()` currently rejects negative log values.

- [ ] **Step 3: Update validation and user-facing scale text**

Require finite log values, but do not require nonnegative log values.

```r
if (any(!is.finite(bulk_log))) {
  stop("Bulk log2(CPM) values must be finite", call. = FALSE)
}
```

Keep `validate_dtangle_version()` and the positional QN split. After alignment, assert:

```r
if (!identical(colnames(t(shared_bulk)), colnames(t(shared_lm22)))) {
  stop("Aligned dtangle gene columns are not identical", call. = FALSE)
}
```

Make `transformed_lm22` the final aligned reference matrix that enters dtangle, after optional QN. Build `overlap_report` in complete LM22 row order and derive overlap counts from its `matched` column. Validate `marker_method` as the supported value `ratio`. Keep marker counts and fitted `gamma` in metadata.

Change CLI help from `log2(TPM + 1)` to `log2(CPM)`. Read the bulk file with first column `gene_name` and LM22 with first column `gene_symbol`. Add `--marker-method` with default `ratio`. Write the final aligned matrices with first column `gene_symbol`. Keep `dtangle_shared_bulk.tsv.gz` as the dtangle-only intersection output. Log the shared-gene count, marker count per cell type, dtangle version, marker settings, QN setting, output paths, and UTC completion time.

- [ ] **Step 4: Run tests, lint, and commit the complete dtangle change set**

Run: `Rscript tests/testthat.R tests/testthat/test-dtangle-stage.R`

Expected: PASS with dtangle 2.0.10.

Run: `Rscript tools/lint_r.R`

Expected: `R lint ok`.

```bash
git add R/dtangle_stage.R scripts/run_dtangle.R tests/testthat/test-dtangle-stage.R
git commit -m "fix: align dtangle inputs on log2 CPM"
```

---

### Task 3: Ten-Group Proportion Processing

**Files:**
- Modify: `R/constants.R`
- Create: `R/proportions.R`
- Create: `scripts/process_proportions.R`
- Create: `tests/testthat/test-proportions.R`

**Interfaces:**
- Consumes: sample-by-22-type dtangle or restart proportions.
- Produces: `lm22_group_map() -> named list(10)`, `combine_lm22_proportions(proportions) -> numeric matrix`, `filter_and_adjust_groups(combined, mean_threshold, zero_floor) -> list(weights, report)`, `process_proportions(proportions, mean_threshold, zero_floor) -> list(original, combined, tca_weights, report)`, and four TSV files.

- [ ] **Step 1: Write failing map, aggregation, and threshold tests**

```r
testthat::test_that("the ten groups partition LM22 exactly once", {
  members <- unlist(lm22_group_map(), use.names = FALSE)
  testthat::expect_length(lm22_group_map(), 10L)
  testthat::expect_setequal(members, lm22_cell_types())
  testthat::expect_false(anyDuplicated(members) > 0L)
  testthat::expect_identical(lm22_group_map()[["Gamma-delta T cells"]], "T cells gamma delta")
  testthat::expect_identical(
    lm22_group_map()[["Mast cells"]],
    c("Mast cells resting", "Mast cells activated")
  )
})

testthat::test_that("filter and zero adjustment follow the cohort rule", {
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

The same file must contain explicit validation tests:

```r
testthat::test_that("22-type inputs require exact columns and valid proportions", {
  valid <- matrix(1 / 22, nrow = 2, ncol = 22,
                  dimnames = list(c("S1", "S2"), lm22_cell_types()))
  testthat::expect_error(
    process_proportions(valid[, -1, drop = FALSE], 0.0001, 1e-6),
    "LM22 columns"
  )
  extra <- cbind(valid, unknown = 0)
  testthat::expect_error(process_proportions(extra, 0.0001, 1e-6), "LM22 columns")
  negative <- valid
  negative[1, 1] <- -1
  testthat::expect_error(process_proportions(negative, 0.0001, 1e-6), "nonnegative")
  nonfinite <- valid
  nonfinite[1, 1] <- Inf
  testthat::expect_error(process_proportions(nonfinite, 0.0001, 1e-6), "finite")
  bad_sum <- valid
  bad_sum[1, 1] <- bad_sum[1, 1] + 1e-4
  testthat::expect_error(process_proportions(bad_sum, 0.0001, 1e-6), "sum to one")
})

testthat::test_that("no lineage is mandatory and TCA still needs two", {
  combined <- matrix(c(0.99995, 0.00005), nrow = 1,
                     dimnames = list("S1", c("CD4 T cells", "Mast cells")))
  testthat::expect_error(filter_and_adjust_groups(combined, 0.01, 1e-6), "two")
})
```

- [ ] **Step 2: Run the focused test and verify failure**

Run: `Rscript tests/testthat.R tests/testthat/test-proportions.R`

Expected: FAIL because the map still has `Other LM22` and proportion functions do not exist.

- [ ] **Step 3: Implement the ten-group map and processing**

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
    "Gamma-delta T cells" = "T cells gamma delta",
    "NK cells" = c("NK cells resting", "NK cells activated"),
    "Monocyte/myeloid" = c(
      "Monocytes", "Macrophages M0", "Macrophages M1", "Macrophages M2"
    ),
    "Neutrophils" = "Neutrophils",
    "Eosinophils" = "Eosinophils",
    "Dendritic cells" = c("Dendritic cells resting", "Dendritic cells activated"),
    "Mast cells" = c("Mast cells resting", "Mast cells activated")
  )
}
```

Aggregate in map order with `vapply()`. Retain groups with `colMeans(combined) >= mean_threshold`. Replace exact zeros only, then renormalize rows.

```r
retained <- colMeans(combined) >= mean_threshold
adjusted <- combined[, retained, drop = FALSE]
adjusted[adjusted == 0] <- zero_floor
weights <- sweep(adjusted, 1L, rowSums(adjusted), "/")
```

The filter report columns must be `cell_group`, `cohort_mean`, `threshold`, `retained`, `filter_reason`, `zero_count_before`, and `zero_floor`. Write `proportions_lm22.tsv`, `proportions_combined.tsv`, `proportions_tca_weights.tsv`, and `cell_group_filter_report.tsv`.

- [ ] **Step 4: Run tests, lint, and commit**

Run: `Rscript tests/testthat.R tests/testthat/test-proportions.R`

Expected: PASS.

Run: `Rscript tools/lint_r.R`

Expected: `R lint ok`.

```bash
git add R/constants.R R/proportions.R scripts/process_proportions.R tests/testthat/test-proportions.R
git commit -m "feat: combine LM22 into major blood lineages"
```

---

### Task 4: Cohort-Wide TCA Fit on the Full Mapped Matrix

**Files:**
- Create: `R/tca_stage.R`
- Create: `scripts/fit_tca.R`
- Create: `tests/testthat/test-tca-stage.R`

**Interfaces:**
- Consumes: complete mapped gene-by-sample `log2(CPM)`, positive sample-by-group weights, optional sample-by-covariate `C2`, CPU count, maximum iterations, shard size, and random seed.
- Produces: `validate_tca_inputs(X, W, C2 = NULL) -> invisible(TRUE)`, `remove_constant_features(X) -> list(matrix, report)`, `fit_tca_stage(X, W, C2, num_cores, max_iters, random_seed, log_file) -> list(model, X, excluded_genes)`, `make_gene_shard_manifest(genes, shard_size) -> tibble`, `write_gene_shards(genes, shard_size, output_dir) -> tibble`, `tca_model.rds`, `tca_model.log`, `tca_expression.tsv.gz`, `tca_excluded_genes.tsv`, `gene_shard_manifest.tsv`, and shard text files.

- [ ] **Step 1: Write failing validation, full-feature, and shard tests**

```r
testthat::test_that("TCA input keeps non-LM22 genes", {
  X <- matrix(rnorm(50), nrow = 5,
              dimnames = list(c("LM1", "LM2", "EXTRA1", "EXTRA2", "EXTRA3"), paste0("S", 1:10)))
  W <- matrix(0.5, nrow = 10, ncol = 2,
              dimnames = list(colnames(X), c("B cells", "CD4 T cells")))
  result <- remove_constant_features(X)
  testthat::expect_identical(rownames(result$matrix), rownames(X))
})

testthat::test_that("TCA rejects reordered samples", {
  X <- matrix(rnorm(40), nrow = 4,
              dimnames = list(paste0("G", 1:4), paste0("S", 1:10)))
  W <- matrix(0.5, nrow = 10, ncol = 2,
              dimnames = list(rev(colnames(X)), c("A", "B")))
  testthat::expect_error(validate_tca_inputs(X, W), "sample order")
})

testthat::test_that("constant genes are removed and shard order is stable", {
  X <- rbind(variable = 1:7, constant = rep(2, 7))
  filtered <- remove_constant_features(X)
  manifest <- make_gene_shard_manifest(rownames(filtered$matrix), 3L)
  testthat::expect_identical(filtered$report$gene_name, "constant")
  testthat::expect_identical(manifest$gene_name, "variable")
})

testthat::test_that("one model fits all genes without refitting weights", {
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

- [ ] **Step 2: Run the focused test and verify failure**

Run: `Rscript tests/testthat.R tests/testthat/test-tca-stage.R`

Expected: FAIL because TCA functions do not exist.

- [ ] **Step 3: Implement validation and constant-gene removal**

Require identical sample IDs and order across `X`, `W`, and `C2`; finite values; strictly positive weights; weight row sums within `1e-8`; at least two groups; no missing covariates; and no intercept column in `C2`. Do not intersect `X` with LM22.

```r
constant <- apply(X, 1L, function(values) length(unique(values)) == 1L)
list(
  matrix = X[!constant, , drop = FALSE],
  report = tibble::tibble(gene_name = rownames(X)[constant], reason = "constant_expression")
)
```

- [ ] **Step 4: Implement one cohort-wide fit and deterministic shards**

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

Fit after constant-gene removal. Write `shard_00001.txt` and later files in full prepared gene order. Do not split samples or fit separate TCA models per shard.

- [ ] **Step 5: Add the CLI and run the small package fit test**

Use `--expression-log`, `--weights`, `--covariates`, `--num-cores`, `--max-iters`, `--shard-size`, `--random-seed`, and `--output-dir`. Read the expression file with first column `gene_name`. Log `scale=log2_cpm`.

Run: `Rscript tests/testthat.R tests/testthat/test-tca-stage.R`

Expected: PASS, including a `TCA::test_data()` fit with `refit_W = FALSE`.

- [ ] **Step 6: Commit the TCA fit stage**

```bash
git add R/tca_stage.R scripts/fit_tca.R tests/testthat/test-tca-stage.R
git commit -m "feat: fit TCA on the full log2 CPM matrix"
```

---

### Task 5: Tensor Extraction, HDF5 Assembly, and QC

**Files:**
- Create: `R/tensor_outputs.R`
- Create: `R/qc.R`
- Create: `scripts/extract_tca_shard.R`
- Create: `scripts/assemble_tca_outputs.R`
- Create: `scripts/build_manifest.R`
- Create: `tests/testthat/test-tensor-outputs.R`

**Interfaces:**
- Consumes: TCA expression matrix, model RDS, one gene list, shard manifest, filtered weights, optional `C2`, and all shard HDF5 files.
- Produces: `extract_tensor_shard(X, model, genes, num_cores, log_file) -> named list of gene-by-sample matrices`, `write_tensor_shard(path, tensor, shard_id) -> invisible(path)`, `assemble_hdf5_shards(shard_paths, manifest, output_dir, pipeline_version, tca_version) -> named character`, per-sample reconstruction metrics, minimal QC plots, provenance, and one final HDF5 per retained group.

- [ ] **Step 1: Write failing order, scale, and error tests**

```r
testthat::test_that("assembled HDF5 uses the log2 CPM contract", {
  samples <- c("S1", "S2")
  shard <- list("B cells" = matrix(1:4, nrow = 2,
                                    dimnames = list(c("G1", "G2"), samples)))
  shard_path <- tempfile(fileext = ".h5")
  write_tensor_shard(shard_path, shard, 1L)
  manifest <- tibble::tibble(
    gene_index = 1:2, gene_name = c("G1", "G2"), shard_id = 1L,
    shard_name = "shard_00001", index_within_shard = 1:2
  )
  output <- assemble_hdf5_shards(
    shard_path, manifest, tempfile(), "test", as.character(packageVersion("TCA"))
  )
  h5 <- hdf5r::H5File$new(output[["B cells"]], mode = "r")
  on.exit(h5$close_all())
  testthat::expect_identical(h5$attr_open("scale")$read(), "log2_cpm")
  testthat::expect_identical(h5[["gene_name"]][], c("G1", "G2"))
})
```

The same file must contain explicit failure tests for malformed shard sets:

```r
testthat::test_that("assembly rejects missing shards and duplicated genes", {
  manifest <- tibble::tibble(
    gene_index = 1:3, gene_name = c("G1", "G1", "G3"),
    shard_id = c(1L, 1L, 2L),
    shard_name = c("shard_00001", "shard_00001", "shard_00002"),
    index_within_shard = c(1L, 2L, 1L)
  )
  testthat::expect_error(
    assemble_hdf5_shards(character(), manifest, tempfile(), "test", "1.2.1"),
    "duplicate.*gene"
  )
  manifest$gene_name[[2L]] <- "G2"
  testthat::expect_error(
    assemble_hdf5_shards(character(), manifest, tempfile(), "test", "1.2.1"),
    "missing.*shard"
  )
})

testthat::test_that("shard writing rejects nonfinite values", {
  bad <- list("B cells" = matrix(
    c(1, Inf), nrow = 1, dimnames = list("G1", c("S1", "S2"))
  ))
  testthat::expect_error(write_tensor_shard(tempfile(fileext = ".h5"), bad, 1L), "finite")
})

testthat::test_that("tensor sources and samples follow the model", {
  set.seed(20260901)
  data <- TCA::test_data(12, 16, 3, 0, 0, 0.01)
  model <- TCA::tca(data$X, data$W, refit_W = FALSE, max_iters = 2, verbose = FALSE)
  tensor <- extract_tensor_shard(data$X, model, rownames(data$X)[1:4], 1L, tempfile())
  testthat::expect_identical(names(tensor), colnames(data$W))
  testthat::expect_true(all(purrr::map_lgl(tensor, ~ identical(colnames(.x), colnames(data$X)))))
})

testthat::test_that("assembly rejects shard sample-order and dimension mismatches", {
  first <- list("B cells" = matrix(
    1:4, nrow = 2, dimnames = list(c("G1", "G2"), c("S1", "S2"))
  ))
  second <- list("B cells" = matrix(
    5:8, nrow = 2, dimnames = list(c("G3", "G4"), c("S2", "S1"))
  ))
  paths <- c(tempfile(fileext = ".h5"), tempfile(fileext = ".h5"))
  write_tensor_shard(paths[[1L]], first, 1L)
  write_tensor_shard(paths[[2L]], second, 2L)
  manifest <- tibble::tibble(
    gene_index = 1:4, gene_name = paste0("G", 1:4),
    shard_id = c(1L, 1L, 2L, 2L),
    shard_name = rep(c("shard_00001", "shard_00002"), each = 2L),
    index_within_shard = rep(1:2, 2L)
  )
  testthat::expect_error(
    assemble_hdf5_shards(paths, manifest, tempfile(), "test", "1.2.1"),
    "sample order"
  )
  second[[1L]] <- second[[1L]][1, , drop = FALSE]
  write_tensor_shard(paths[[2L]], second, 2L)
  testthat::expect_error(
    assemble_hdf5_shards(paths, manifest, tempfile(), "test", "1.2.1"),
    "dimension"
  )
})
```

- [ ] **Step 2: Run the focused test and verify failure**

Run: `Rscript tests/testthat.R tests/testthat/test-tensor-outputs.R`

Expected: FAIL because tensor functions do not exist.

- [ ] **Step 3: Implement tensor extraction and streaming assembly**

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

Write shard data under `/sources/<slug>/expression`. Assemble one final file per group with `/expression`, `/gene_name`, `/sample_id`, and attributes `cell_group`, `scale = log2_cpm`, `pipeline_version`, and `tca_version`. Use chunks of at most 500 genes by 256 samples and gzip level 6. Stream optional TSV rows through `gzfile()`.

- [ ] **Step 4: Implement reconstruction metrics and reports**

For each shard, calculate the proportion-weighted sum of source tensors and add the fitted `C2` term when present. Accumulate sample-wise sums needed for Pearson correlation and RMSE. Write `reconstruction_by_sample.tsv`, `qc_summary.tsv`, `qc_plots.pdf`, and `output_manifest.json`.

```r
source_reconstruction <- Reduce(`+`, purrr::map2(
  tensor,
  seq_len(ncol(weights)),
  ~ sweep(.x, 2L, weights[, .y], "*")
))
```

Use `ggplot2::theme_minimal()` without a title or subtitle. Include a cell-group proportion distribution, observed-versus-reconstructed comparison, and residual distribution. The manifest must contain file SHA-256, dimensions, scale, cell group, software versions, parameters, and the container image.

- [ ] **Step 5: Add the three CLIs, run tests, and commit**

Run: `Rscript tests/testthat.R tests/testthat/test-tensor-outputs.R`

Expected: PASS with exact gene and sample order and finite reconstruction metrics.

```bash
git add R/tensor_outputs.R R/qc.R scripts/extract_tca_shard.R scripts/assemble_tca_outputs.R scripts/build_manifest.R tests/testthat/test-tensor-outputs.R
git commit -m "feat: extract cell-type-specific log2 CPM matrices"
```

---

### Task 6: Logged Modular WDL Tasks

**Files:**
- Create: `workflows/tasks/expression.wdl`
- Create: `workflows/tasks/dtangle.wdl`
- Create: `workflows/tasks/proportions.wdl`
- Create: `workflows/tasks/tca.wdl`
- Create: `workflows/tasks/qc.wdl`
- Create: `tests/testthat/test-wdl-contract.R`
- Create: `tools/check_wdl_logging.py`

**Interfaces:**
- Consumes: R CLIs from Tasks 1 through 5 and a task-level `docker_image`.
- Produces: WDL 1.1 tasks `PrepareExpression`, `RunDtangle`, `ProcessProportions`, `FitTca`, `ExtractTcaShard`, `AssembleTca`, and `BuildManifest`.

- [ ] **Step 1: Write failing task and logging contract tests**

```r
testthat::test_that("the expression WDL requires CPM and GTF", {
  text <- paste(readLines("workflows/tasks/expression.wdl"), collapse = "\n")
  testthat::expect_match(text, "File expression")
  testthat::expect_match(text, "File gtf")
  testthat::expect_match(text, "prepared_cpm")
  testthat::expect_match(text, "prepared_log2_cpm")
  testthat::expect_false(grepl("expression_type|gene_length", text))
})
```

Implement `tools/check_wdl_logging.py` to fail any `command <<< ... >>>` block that lacks `stage=`, `start_time=`, `completion_time=`, `dimensions=`, and `outputs=`.

- [ ] **Step 2: Run the contract test and verify failure**

Run: `Rscript tests/testthat.R tests/testthat/test-wdl-contract.R`

Expected: FAIL because WDL task files do not exist.

- [ ] **Step 3: Implement exact WDL task interfaces**

```text
PrepareExpression:
  in  expression, gtf, docker_image, cpu, memory, disk_gb
  out prepared_cpm, prepared_log2_cpm, mapping_report, excluded_genes, log
RunDtangle:
  in  prepared_log2_cpm, lm22, min_overlap, marker_fraction, marker_method,
      quantile_normalize, docker_image, cpu, memory, disk_gb
  out proportions, markers, metadata, overlap_report, transformed_lm22,
      shared_bulk, log
ProcessProportions:
  in  proportions, mean_threshold, zero_floor, docker_image, cpu, memory, disk_gb
  out original, combined, tca_weights, filter_report, log
FitTca:
  in  prepared_log2_cpm, tca_weights, covariates?, shard_size, max_iters,
      random_seed, docker_image, cpu, memory, disk_gb
  out model, model_log, tca_expression, excluded_genes, shard_manifest,
      Array[File] shards, log
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
  in  primary and supporting files, scientific parameters, container image
  out output_manifest, qc_summary, qc_plots, provenance, log
```

- [ ] **Step 4: Add logging and Terra run-time controls**

Use this pattern in each command with a task-specific stage:

```wdl
command <<<
  set -euo pipefail
  stage="prepare_expression"
  log="$stage.log"
  printf 'stage=%s start_time=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  trap 'status=$?; printf "stage=%s error_status=%s time=%s\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
  Rscript /opt/celltype/scripts/prepare_expression.R \
    --expression '~{expression}' --gtf '~{gtf}' --output-dir outputs 2>&1 | tee -a "$log"
  printf 'stage=%s dimensions=%s outputs=%s completion_time=%s\n' \
    "$stage" "$(wc -l < outputs/gene_mapping_report.tsv)" "outputs" \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
>>>
```

Every task accepts `preemptible_attempts` and `max_retries`. Use the approved starting profile: PrepareExpression 4 CPU/64 GB/400 GB; RunDtangle 4/32 GB/100 GB; ProcessProportions 2/16 GB/50 GB; FitTca 16/192 GB/750 GB; ExtractTcaShard 8/64 GB/200 GB; AssembleTca 8/128 GB/500 GB; BuildManifest 4/32 GB/100 GB.

- [ ] **Step 5: Validate WDL, logging, tests, and commit**

Run: `miniwdl check workflows/tasks/expression.wdl workflows/tasks/dtangle.wdl workflows/tasks/proportions.wdl workflows/tasks/tca.wdl workflows/tasks/qc.wdl`

Expected: No syntax or type errors.

Run: `python3 tools/check_wdl_logging.py workflows`

Expected: `All WDL command blocks contain required logging fields.`

Run: `Rscript tests/testthat.R tests/testthat/test-wdl-contract.R`

Expected: PASS.

```bash
git add workflows/tasks tests/testthat/test-wdl-contract.R tools/check_wdl_logging.py
git commit -m "feat: add logged CPM and TCA WDL tasks"
```

---

### Task 7: Top-Level Workflow and Deterministic Smoke Fixtures

**Files:**
- Create: `workflows/cell_type_deconvolution.wdl`
- Create: `scripts/generate_synthetic_fixture.R`
- Create: generated files under `tests/fixtures/`
- Create: `tests/fixtures/dtangle.inputs.json`
- Create: `tests/fixtures/restart.inputs.json`
- Create: `tests/smoke/assert_outputs.R`
- Modify: `tests/testthat/test-wdl-contract.R`

**Interfaces:**
- Consumes: CPM, GTF, optional LM22, optional precomputed 22-type proportions, optional `C2`, scientific thresholds, and run-time settings.
- Produces: all primary and supporting outputs through workflow `cell_type_deconvolution`.

- [ ] **Step 1: Add a failing top-level workflow contract test**

```r
testthat::test_that("top-level workflow requires CPM and GTF and exposes both proportion modes", {
  text <- paste(readLines("workflows/cell_type_deconvolution.wdl"), collapse = "\n")
  testthat::expect_match(text, "File expression")
  testthat::expect_match(text, "File gtf")
  testthat::expect_match(text, "File[?] lm22")
  testthat::expect_match(text, "File[?] precomputed_proportions")
  testthat::expect_match(text, "scatter .*shard")
  testthat::expect_false(grepl("expression_type|gene_length", text))
})
```

- [ ] **Step 2: Run the contract test and verify failure**

Run: `Rscript tests/testthat.R tests/testthat/test-wdl-contract.R`

Expected: FAIL because the top-level workflow does not exist.

- [ ] **Step 3: Implement the top-level branch and full-matrix flow**

Always call `PrepareExpression`. Use `PrepareExpression.prepared_log2_cpm` for both `RunDtangle` and `FitTca`; do not pass `RunDtangle.shared_bulk` to TCA.

```wdl
Boolean estimate_proportions = !defined(precomputed_proportions)
String dtangle_marker_method = "ratio"

if (estimate_proportions) {
  call dtangle_tasks.RunDtangle {
    input:
      prepared_log2_cpm = PrepareExpression.prepared_log2_cpm,
      lm22 = select_first([lm22]),
      min_overlap = min_lm22_overlap,
      marker_fraction = dtangle_marker_fraction,
      marker_method = dtangle_marker_method,
      quantile_normalize = dtangle_quantile_normalize,
      docker_image = docker_image
  }
}

File proportions_for_processing = select_first([
  precomputed_proportions,
  RunDtangle.proportions
])
```

Then call proportion processing, the one TCA fit, the gene-shard scatter, assembly, and manifest creation. Expose the three proportion tables, group HDF5 files, optional TSV files, model, TCA model log, QC, reports, command logs, and manifest.

- [ ] **Step 4: Generate synthetic CPM, GTF, and both mode inputs**

Use seed `20260901`. Generate a positive 66-gene, 22-type synthetic signature with three marker genes per type. Generate twelve positive bulk CPM samples plus six non-signature genes. Scale each sample to one million without later pipeline renormalization. Generate a GTF that includes all expression gene IDs, both protein-coding and non-protein-coding gene types, and no licensed LM22 values. Generate one `batch_indicator` covariate table and a precomputed 22-type proportion table.

```r
bulk_linear <- rbind(signature %*% t(weights), extra_gene_expression)
bulk_cpm <- sweep(bulk_linear, 2L, colSums(bulk_linear), "/") * 1e6
stopifnot(all(bulk_cpm > 0), max(abs(colSums(bulk_cpm) - 1e6)) < 1e-6)
```

Run: `Rscript scripts/generate_synthetic_fixture.R tests/fixtures`

Expected: deterministic CPM, GTF, signature, covariate, proportion, expected-ID, and two JSON files. `restart.inputs.json` must omit LM22 and use the precomputed 22-type table.

- [ ] **Step 5: Implement smoke assertions**

```r
stopifnot(all(is.finite(proportions)), all(proportions >= 0))
stopifnot(max(abs(rowSums(proportions) - 1)) < 1e-8)
stopifnot(all(tca_weights > 0))
stopifnot(max(abs(rowSums(tca_weights) - 1)) < 1e-8)
stopifnot(identical(expected_samples, observed_samples))
stopifnot(identical(expected_full_mapped_genes, observed_genes))
stopifnot(all(is.finite(reconstruction$correlation)))
stopifnot(all(is.finite(reconstruction$rmse)))
```

Open every final HDF5 and check dimensions, identifiers, `cell_group`, and `scale = log2_cpm`. Confirm that at least one non-LM22 gene is present in every TCA output.

- [ ] **Step 6: Validate and commit**

Run: `miniwdl check workflows/cell_type_deconvolution.wdl`

Expected: No syntax, import, or type errors.

```bash
git add workflows/cell_type_deconvolution.wdl scripts/generate_synthetic_fixture.R tests/fixtures tests/smoke tests/testthat/test-wdl-contract.R
git commit -m "feat: add complete CPM deconvolution workflow"
```

---

### Task 8: Pinned Micromamba Image and GitHub Actions Smoke Test

**Files:**
- Create: `envs/environment.yml`
- Modify: `envs/Dockerfile`
- Create: `.github/workflows/pipeline-ci.yml`
- Modify: `.github/workflows/docker-image.yml`
- Delete: `.github/workflows/r-lint.yml`
- Delete: `.github/workflows/wdl-validation.yml`
- Modify: `README.md` badge block by running `tools/update_workflow_badges.py`

**Interfaces:**
- Consumes: repository source and synthetic input JSON files.
- Produces: a pinned image, pull-request R/WDL/smoke evidence, and GHCR images on `main`.

- [ ] **Step 1: Create the exact environment input**

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
  - r-testthat=3.3.1
  - r-lintr=3.4.0
  - bioconductor-limma=3.66.0
```

- [ ] **Step 2: Implement the pinned micromamba Dockerfile**

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

- [ ] **Step 3: Add the GitHub Actions build and two-mode smoke job**

`pipeline-ci.yml` must check out the repository, install `miniwdl==1.15.0`, validate WDL and logging, build `celltype-deconvolution:test`, run R lint and tests in that image, run dtangle mode, assert its outputs, run restart mode, assert its outputs, and upload logs and QC files on failure. Both fixture JSON files must use `docker_image = "celltype-deconvolution:test"`.

Use GitHub Actions YAML commands with exact checks:

```yaml
- name: Static checks
  run: |
    miniwdl check workflows/cell_type_deconvolution.wdl
    python3 tools/check_wdl_logging.py workflows
- name: R checks
  run: |
    docker run --rm -v "$PWD:/workspace" celltype-deconvolution:test Rscript tools/lint_r.R
    docker run --rm -v "$PWD:/workspace" celltype-deconvolution:test Rscript tests/testthat.R
```

- [ ] **Step 4: Update image publication and consolidate old checks**

Set the main-branch image job to `push: true`. Publish `sha-<commit>` and `latest` tags with OCI source, revision, and created labels. Delete the old R-lint and WDL-validation workflows only after their checks exist in `pipeline-ci.yml`. Run `python3 tools/update_workflow_badges.py`.

- [ ] **Step 5: Run static checks without a local Docker build and commit**

Run: `miniwdl check workflows/cell_type_deconvolution.wdl`

Run: `python3 tools/check_wdl_logging.py workflows`

Run: `git diff --check`

Expected: all commands pass. Do not run `docker build` locally.

```bash
git add envs .github/workflows README.md
git commit -m "ci: build and smoke test the CPM pipeline"
```

---

### Task 9: Terra, Dockstore, and Data Documentation

**Files:**
- Modify: `README.md`
- Create: `.dockstore.yml`
- Create: `docs/terra.md`
- Create: `docs/data-dictionary.md`
- Create: `examples/cpm.inputs.json`
- Create: `examples/precomputed-proportions.inputs.json`
- Create: `tests/testthat/test-documentation.R`

**Interfaces:**
- Consumes: final WDL input and output names.
- Produces: Dockstore metadata, two Terra-ready examples, and exact input/output documentation.

- [ ] **Step 1: Write failing documentation contract tests**

```r
testthat::test_that("both examples use CPM and GTF workflow inputs", {
  files <- list.files("examples", pattern = "[.]inputs[.]json$", full.names = TRUE)
  testthat::expect_length(files, 2L)
  parsed <- purrr::map(files, jsonlite::read_json)
  testthat::expect_true(all(purrr::map_lgl(parsed, ~ all(c(
    "cell_type_deconvolution.expression",
    "cell_type_deconvolution.gtf"
  ) %in% names(.x)))))
})

testthat::test_that("the scale and LM22 rules are explicit", {
  text <- paste(readLines("docs/terra.md"), collapse = "\n")
  testthat::expect_match(text, "user-supplied")
  testthat::expect_match(text, "linear LM22")
  testthat::expect_match(text, "log2\\(CPM\\)")
  testthat::expect_match(text, "strictly positive")
  testthat::expect_false(grepl("expression_type|log2_tpm", text))
})
```

- [ ] **Step 2: Run the documentation test and verify failure**

Run: `Rscript tests/testthat.R tests/testthat/test-documentation.R`

Expected: FAIL because the two examples and guides do not exist.

- [ ] **Step 3: Write Dockstore metadata and user guides**

```yaml
version: 1.2
workflows:
  - subclass: WDL
    primaryDescriptorPath: /workflows/cell_type_deconvolution.wdl
    testParameterFiles:
      - /tests/fixtures/dtangle.inputs.json
    name: cell-type-deconvolution
```

README and Terra guide sections must cover: positive normalized CPM, required GTF, no gene-type filter, exact ID matching, user-supplied linear LM22, dtangle and restart modes, marker method `ratio`, ten-group mapping, `0.0001` mean filter, full-matrix TCA, `log2_cpm` outputs, 9,000-sample run-time defaults, HDF5 layout, GitHub Actions smoke testing, and LM22 licensing.

The limitations section must state that LM22 estimates relative immune proportions, was derived from microarray data, does not model erythrocyte or platelet expression, and includes macrophage states in the monocyte/myeloid group. It must also state that the default path has no QN, TCA outputs are statistical estimates, zeros stop CPM preparation, ID matching is exact after trimming, retained groups can differ by cohort, and very small proportions can be unstable.

- [ ] **Step 4: Write the data dictionary and exact examples**

For every input and output, record WDL name, type, required mode, default, scale, validation, and description. State that `dtangle_shared_bulk` is LM22-limited and that `tca_expression` is the full mapped matrix after constant-gene removal. Make both JSON namespaces match `cell_type_deconvolution` exactly.

- [ ] **Step 5: Run checks and commit**

Run: `Rscript tests/testthat.R tests/testthat/test-documentation.R`

Expected: PASS.

Run: `miniwdl check workflows/cell_type_deconvolution.wdl`

Expected: PASS.

Run: `git diff --check`

Expected: no output.

```bash
git add README.md .dockstore.yml docs/terra.md docs/data-dictionary.md examples tests/testthat/test-documentation.R
git commit -m "docs: add Terra guidance for CPM and GTF inputs"
```

---

### Task 10: Final Verification and GitHub Actions Evidence

**Files:**
- Modify only files required by verified failures.

**Interfaces:**
- Consumes: complete implementation branch.
- Produces: clean local checks and passing GitHub Actions evidence for image construction and both workflow modes.

- [ ] **Step 1: Run all non-container local checks**

Run: `Rscript tests/testthat.R`

Expected: all tests pass when the pinned R environment is active.

Run: `Rscript tools/lint_r.R`

Expected: `R lint ok`.

Run: `miniwdl check workflows/cell_type_deconvolution.wdl`

Expected: no errors.

Run: `python3 tools/check_wdl_logging.py workflows`

Expected: all command blocks pass.

Run: `git diff --check`

Expected: no output.

- [ ] **Step 2: Confirm no LM22 data or secrets are tracked**

Run: `git ls-files | rg -i 'lm22[.]txt|token|credential|secret'`

Expected: no LM22 matrix, token, credential, or secret file. Documentation references are acceptable when the file path does not match `LM22.txt`.

- [ ] **Step 3: Push the branch and wait for GitHub Actions**

Run: `git push -u origin feat/dtangle-tca-pipeline`

Run: `gh run list --workflow pipeline-ci.yml --branch feat/dtangle-tca-pipeline --limit 1`

Expected: one completed successful run.

- [ ] **Step 4: Inspect and record the workflow evidence**

Run: `gh run view "$(gh run list --workflow pipeline-ci.yml --branch feat/dtangle-tca-pipeline --limit 1 --json databaseId --jq '.[0].databaseId')" --log-failed`

Expected: no failed logs. The summary must show the image build, R tests, R lint, WDL validation, dtangle-mode smoke, restart-mode smoke, and output assertions.

Record the Actions URL and tested image identifier in the final handoff. Do not claim the pipeline passed until this run is successful.

- [ ] **Step 5: Commit only verified corrections when needed**

If a verification failure requires a correction, rerun the focused check and every check in Step 1, then commit the exact changed files:

```bash
git add R scripts workflows tests envs .github README.md docs examples .dockstore.yml tools
git commit -m "fix: resolve pipeline verification failures"
```

If verification changes no file, do not create a commit.
