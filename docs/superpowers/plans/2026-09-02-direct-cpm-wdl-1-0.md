# Direct CPM WDL 1.0 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove the expression-preparation task and run dtangle, TCA, and BED export directly from one coordinate-preserving linear-CPM BED under WDL 1.0.

**Architecture:** Shared R readers validate the BED once per consuming task. `RunDtangle` maps gene IDs through the GTF and creates its gene-symbol log2-CPM view. `FitTca` creates the full gene-ID log2-CPM view, and `ExportTcaBeds` reads coordinates from the same BED. The workflow has one default GitHub image and one global setting each for preemptible attempts and retries.

**Tech Stack:** WDL 1.0, MiniWDL 1.15.0, R 4.3, tidyverse, dtangle 2.0.10, TCA 1.2.1, testthat, GitHub Actions, Terra/Cromwell

**Spec:** `docs/superpowers/specs/2026-09-02-direct-cpm-wdl-1-0-design.md`

## Global Constraints

- Every active workflow and task WDL file must declare `version 1.0`.
- Input expression is a BED or BED.GZ with `#chr`, `start`, `end`, `gene_id`, then sample columns.
- Expression values are finite, strictly positive linear CPM values.
- Apply `log2()` once without a pseudocount inside dtangle and TCA consumers.
- TCA uses all valid nonconstant gene-ID rows; GTF mapping only controls dtangle inputs.
- Default image is `ghcr.io/aou-multiomics-analysis/celltypedeconvolution:latest`.
- Workflow defaults are `preemptible_attempts = 2` and `max_retries = 2`, applied to every task.
- Keep task-specific CPU, memory, and disk inputs.
- Remove the public `pipeline_version` input and manifest field.
- Keep start, dimension, output, completion, and error logging in every WDL task.
- Do not build Docker locally. GitHub Actions performs the image and workflow smoke tests.

---

### Task 1: Define direct CPM views in the shared R layer

**Files:**
- Modify: `R/expression.R`
- Modify: `R/expression_bed.R`
- Modify: `tests/testthat/test-expression.R`

**Interfaces:**
- Consumes: `read_expression_bed(path) -> list(coordinates, cpm)` and `read_gtf_gene_annotation(path) -> tibble`
- Produces: `make_dtangle_expression(expression, annotation) -> list(log_expression, mapping_report)` and `make_tca_expression(expression) -> numeric matrix`

- [ ] **Step 1: Write failing direct-view tests**

Add tests that prove TCA keeps all BED genes while dtangle maps and aggregates symbols:

```r
testthat::test_that("direct CPM views separate TCA genes from dtangle symbols", {
  expression <- list(
    coordinates = tibble::tibble(
      `#chr` = c("chr1", "chr1", "chr2"),
      start = c(10L, 20L, 30L),
      end = c(11L, 21L, 31L),
      gene_id = c("ENSG1", "ENSG2", "ENSG3")
    ),
    cpm = matrix(
      c(4, 8, 16, 32, 64, 128),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(c("ENSG1", "ENSG2", "ENSG3"), c("S1", "S2"))
    )
  )
  annotation <- tibble::tibble(
    gene_id = c("ENSG1", "ENSG2"),
    gene_name = c("MARKER", "MARKER"),
    gene_type = c("protein_coding", "lncRNA")
  )

  tca_expression <- make_tca_expression(expression)
  dtangle_expression <- make_dtangle_expression(expression, annotation)

  testthat::expect_identical(rownames(tca_expression), c("ENSG1", "ENSG2", "ENSG3"))
  testthat::expect_equal(unname(tca_expression[, "S1"]), log2(c(4, 16, 64)))
  testthat::expect_identical(rownames(dtangle_expression$log_expression), "MARKER")
  testthat::expect_equal(unname(dtangle_expression$log_expression[1, ]), log2(c(20, 40)))
})
```

Add one test that confirms zero, negative, missing, and nonfinite BED values still fail through `read_expression_bed()`.

- [ ] **Step 2: Run the focused test and verify RED**

Run:

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-expression.R")'
```

Expected: FAIL because `make_dtangle_expression()` and `make_tca_expression()` do not exist.

- [ ] **Step 3: Implement the direct-view functions**

Replace preparation-specific composition with these focused interfaces:

```r
make_tca_expression <- function(expression) {
  if (!is.list(expression) || !all(c("coordinates", "cpm") %in% names(expression))) {
    stop("expression must contain coordinates and cpm", call. = FALSE)
  }
  cpm <- validate_cpm_matrix(expression$cpm)
  if (!identical(rownames(cpm), expression$coordinates$gene_id)) {
    stop("Expression genes and BED coordinates must match exactly", call. = FALSE)
  }
  log2(cpm)
}

make_dtangle_expression <- function(expression, annotation) {
  if (!is.list(expression) || !all(c("coordinates", "cpm") %in% names(expression))) {
    stop("expression must contain coordinates and cpm", call. = FALSE)
  }
  cpm <- validate_cpm_matrix(expression$cpm)
  annotation <- validate_gtf_gene_annotation(annotation)
  mapping_report <- make_cpm_mapping_report(cpm, annotation)
  symbol_cpm <- collapse_cpm_to_gene_names(cpm, annotation)
  if (nrow(symbol_cpm) == 0L || any(symbol_cpm <= 0)) {
    stop("No positive CPM genes map to a usable gene symbol", call. = FALSE)
  }
  list(log_expression = log2(symbol_cpm), mapping_report = mapping_report)
}
```

Keep `prepare_expression()` and `prepare_expression_bed()` temporarily because the existing WDL still calls the preparation script. Task 5 removes both functions with the obsolete task and script. Keep `read_expression_bed()`, GTF parsing, mapping, validation, and duplicate-symbol aggregation.

- [ ] **Step 4: Run the focused tests and verify GREEN**

Run:

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-expression.R")'
```

Expected: PASS with no failures.

- [ ] **Step 5: Commit the shared data contract**

```bash
git add R/expression.R R/expression_bed.R tests/testthat/test-expression.R
git commit -m "refactor: define direct CPM expression views"
```

---

### Task 2: Make dtangle consume the CPM BED and GTF directly

**Files:**
- Modify: `scripts/run_dtangle.R`
- Modify: `tests/testthat/test-dtangle.R`

**Interfaces:**
- Consumes: expression BED plus GTF or the temporary prepared-log compatibility input, LM22, and `make_dtangle_expression()` from Task 1
- Produces: the existing dtangle proportions, markers, metadata, overlap report, transformed LM22, shared bulk, and task log; Task 5 removes the temporary compatibility input

- [ ] **Step 1: Write failing CLI and WDL tests**

Add assertions that the CLI accepts `--expression` and `--gtf`:

```r
testthat::test_that("dtangle CLI reads direct CPM BED and GTF", {
  text <- paste(readLines("scripts/run_dtangle.R", warn = FALSE), collapse = "\n")
  testthat::expect_match(text, '"--expression"', fixed = TRUE)
  testthat::expect_match(text, '"--gtf"', fixed = TRUE)
})
```

- [ ] **Step 2: Run tests and verify RED**

Run:

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-dtangle.R")'
```

Expected: FAIL because the CLI does not accept the direct inputs.

- [ ] **Step 3: Update the dtangle CLI**

Add these options alongside the temporary `--bulk-log` compatibility option:

```r
optparse::make_option(
  "--expression",
  type = "character",
  help = "Coordinate-preserving BED of positive linear CPM values."
)
optparse::make_option(
  "--gtf",
  type = "character",
  help = "GTF used to map expression gene_id values to LM22 gene symbols."
)
```

Require either `--bulk-log` or the pair `--expression` and `--gtf`. Read and transform direct inputs with:

```r
expression <- read_expression_bed(options$expression)
annotation <- read_gtf_gene_annotation(options$gtf)
dtangle_expression <- make_dtangle_expression(expression, annotation)
bulk_log <- dtangle_expression$log_expression
```

Keep the existing LM22 validation, overlap threshold, marker selection, dtangle fit, outputs, and logging.

- [ ] **Step 4: Run the focused tests**

Run:

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-dtangle.R")'
```

Expected: PASS.

- [ ] **Step 5: Commit direct dtangle input**

```bash
git add scripts/run_dtangle.R tests/testthat/test-dtangle.R
git commit -m "refactor: run dtangle from the CPM BED"
```

---

### Task 3: Make TCA fit and export consume the CPM BED directly

**Files:**
- Modify: `scripts/fit_tca.R`
- Modify: `scripts/export_tca_beds.R`
- Modify: `tests/testthat/test-tca.R`
- Modify: `tests/testthat/test-bed-outputs.R`

**Interfaces:**
- Consumes: expression BED or temporary prepared-file compatibility inputs, TCA weights, optional covariates, and model settings
- Produces: unchanged TCA model, model log, filtered log2-CPM expression, excluded genes, cell-type BEDs, reconstruction QC, inventory, plots, and logs; Task 5 removes the temporary compatibility inputs

- [ ] **Step 1: Write failing direct TCA tests**

Add a fit CLI contract test:

```r
testthat::test_that("TCA fit reads a direct CPM BED", {
  text <- paste(readLines("scripts/fit_tca.R", warn = FALSE), collapse = "\n")
  testthat::expect_match(text, '"--expression"', fixed = TRUE)
  testthat::expect_match(text, "make_tca_expression(expression)", fixed = TRUE)
})
```

Add an export test that writes a small BED, reads its coordinates, aligns a subset of modeled genes, and expects the original coordinate order.

- [ ] **Step 2: Run focused tests and verify RED**

Run:

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-tca.R"); testthat::test_file("tests/testthat/test-bed-outputs.R")'
```

Expected: FAIL because neither CLI accepts the direct BED.

- [ ] **Step 3: Update the TCA fit CLI**

Add the direct option while retaining `--expression-log` temporarily:

```r
optparse::make_option(
  "--expression",
  type = "character",
  help = "Coordinate-preserving BED of positive linear CPM values."
)
```

Require exactly one expression source. Create `X` from the direct source with:

```r
expression <- read_expression_bed(options$expression)
X <- make_tca_expression(expression)
```

Keep `fit_tca_stage()` responsible for sample alignment, constant-gene exclusion, fitting, and model logging.

- [ ] **Step 4: Update direct coordinate extraction for export**

Add `--expression` while retaining `--coordinates` temporarily. Require one coordinate source. Read the direct BED and align coordinates to the filtered TCA expression:

```r
expression <- read_expression_bed(options$expression)
coordinates <- expression$coordinates
coordinate_index <- match(rownames(X), coordinates$gene_id)
if (anyNA(coordinate_index)) {
  stop("Every modeled gene_id must have BED coordinates", call. = FALSE)
}
coordinates <- coordinates[coordinate_index, , drop = FALSE]
```

Task 5 removes the duplicate `read_bed_coordinates()` implementation after the WDL switches to direct input.

- [ ] **Step 5: Run focused tests**

Run:

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-tca.R"); testthat::test_file("tests/testthat/test-bed-outputs.R")'
```

Expected: PASS.

- [ ] **Step 6: Commit direct TCA input**

```bash
git add scripts/fit_tca.R scripts/export_tca_beds.R tests/testthat/test-tca.R tests/testthat/test-bed-outputs.R
git commit -m "refactor: run TCA from the CPM BED"
```

---

### Task 4: Make preparation provenance optional for the interface transition

**Files:**
- Modify: `R/qc.R`
- Modify: `scripts/build_manifest.R`
- Modify: `tests/testthat/test-bed-outputs.R`

**Interfaces:**
- Consumes: export QC, proportions, TCA model and log, optional preparation mapping, optional pipeline version, optional dtangle metadata, effective parameters, and container image
- Produces: a manifest and QC summary that omit preparation fields when their optional inputs are absent while remaining compatible with the current WDL until Task 5

- [ ] **Step 1: Write failing manifest tests**

Add calls that omit both optional provenance inputs:

```r
testthat::expect_false("pipeline_version" %in% names(manifest))
testthat::expect_false("duplicate_gene_symbol_input_row_count" %in% qc_summary$metric)
testthat::expect_false("duplicate_gene_symbol_count" %in% qc_summary$metric)
```

Keep one existing compatibility call with mapping data and a pipeline version. It must retain the current fields until Task 5 switches the WDL atomically.

- [ ] **Step 2: Run focused tests and verify RED**

Run:

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-bed-outputs.R")'
```

Expected: FAIL because the manifest and QC functions still require preparation metadata and pipeline version.

- [ ] **Step 3: Simplify manifest functions**

Change the manifest signature to:

```r
build_output_manifest <- function(
    outputs,
    tca_version,
    parameters,
    container_image,
    pipeline_version = NULL)
```

Return `schema_version`, `created_utc`, `tca_version`, `software_versions`, `parameters`, `container_image`, and `outputs`. Add `pipeline_version` only when the optional value is nonempty. Task 5 removes the final caller.

Change the QC signature to:

```r
build_pipeline_qc_summary <- function(
    export_summary,
    original_proportions,
    combined_proportions,
    tca_weights,
    filter_report,
    tca_model,
    tca_log_lines,
    dtangle_metadata = NULL,
    mapping_report = NULL)
```

Add the two duplicate-symbol metrics only when `mapping_report` is supplied. Keep proportion, zero-adjustment, reconstruction, LM22, TCA iteration, convergence, and `tau_hat` metrics in both modes.

- [ ] **Step 4: Make the manifest CLI inputs optional**

Give `--mapping-report` and `--pipeline-version` `NULL` defaults in `scripts/build_manifest.R`, remove them from `required_options`, and read them only when supplied. Keep the current WDL task compatible until Task 5.

Allow the GitHub default image in `validate_container_image()`:

```r
default_image <- "ghcr.io/aou-multiomics-analysis/celltypedeconvolution:latest"
valid <- is.character(container_image) && length(container_image) == 1L &&
  !is.na(container_image) &&
  (identical(container_image, local_smoke_image) ||
    identical(container_image, default_image) ||
    grepl(digest_pattern, container_image))
```

- [ ] **Step 5: Run focused tests**

Run:

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-bed-outputs.R")'
```

Expected: PASS.

- [ ] **Step 6: Commit the provenance simplification**

```bash
git add R/qc.R scripts/build_manifest.R tests/testthat/test-bed-outputs.R
git commit -m "refactor: remove expression preparation provenance"
```

---

### Task 5: Rewrite the top workflow as WDL 1.0 with global runtime settings

**Files:**
- Modify: `workflows/cell_type_deconvolution.wdl`
- Modify: `workflows/tasks/dtangle.wdl`
- Modify: `workflows/tasks/proportions.wdl`
- Modify: `workflows/tasks/tca.wdl`
- Modify: `workflows/tasks/qc.wdl`
- Modify: `R/expression.R`
- Modify: `R/qc.R`
- Modify: `scripts/build_manifest.R`
- Modify: `scripts/run_dtangle.R`
- Modify: `scripts/fit_tca.R`
- Modify: `scripts/export_tca_beds.R`
- Delete: `scripts/prepare_expression.R`
- Delete: `workflows/tasks/expression.wdl`
- Modify: `tests/testthat/test-wdl-contract.R`

**Interfaces:**
- Consumes: one CPM BED, one GTF, exactly one proportion source, optional covariates, scientific parameters, task resources, and two global retry controls
- Produces: proportion, model, cell-type BED, QC, manifest, inventory, and log outputs

- [ ] **Step 1: Write failing WDL 1.0 and global-setting tests**

Add this contract test:

```r
testthat::test_that("all active WDL files use WDL 1.0", {
  files <- list.files("workflows", pattern = "[.]wdl$", recursive = TRUE, full.names = TRUE)
  declarations <- vapply(files, function(path) readLines(path, n = 1L, warn = FALSE), character(1))
  testthat::expect_true(all(declarations == "version 1.0"), info = paste(files, declarations))
})
```

Add assertions for:

```r
testthat::expect_match(text, 'String docker_image = "ghcr.io/aou-multiomics-analysis/celltypedeconvolution:latest"', fixed = TRUE)
testthat::expect_match(text, "Int preemptible_attempts = 2", fixed = TRUE)
testthat::expect_match(text, "Int max_retries = 2", fixed = TRUE)
testthat::expect_false(grepl("prepare_preemptible_attempts|fit_max_retries|pipeline_version", text))
testthat::expect_false(grepl("PrepareExpression|expression_tasks", text))
```

- [ ] **Step 2: Run the WDL contract test and verify RED**

Run:

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-wdl-contract.R")'
```

Expected: FAIL on WDL 1.1, preparation references, required Docker input, and task-specific retry controls.

- [ ] **Step 3: Rewrite the top workflow data flow**

Remove the expression-task import, call, and outputs. Use these call bindings:

```wdl
expression = expression,
gtf = gtf
```

for `RunDtangle`, and:

```wdl
expression = expression
```

for `FitTca` and `ExportTcaBeds`.

Remove `mapping_report`, `excluded_genes`, and `pipeline_version` from `BuildManifest`.

Change `RunDtangle` to require `File expression` and `File gtf`, then invoke `--expression` and `--gtf`. Change `FitTca` to require `File expression`. Change `ExportTcaBeds` to require `File expression` instead of `File coordinates`.

Delete `scripts/prepare_expression.R` and `workflows/tasks/expression.wdl`. Remove `prepare_expression()` and `prepare_expression_bed()` from `R/expression.R`. Remove the transitional `--bulk-log`, `--expression-log`, and `--coordinates` options from their CLIs. Remove the duplicate `read_bed_coordinates()` function. Remove `--mapping-report` and `--pipeline-version` from `scripts/build_manifest.R`. Remove the transitional `pipeline_version` and `mapping_report` parameters and conditional fields from `R/qc.R`. Keep the shared BED, GTF, mapping, validation, and direct-view functions.

- [ ] **Step 4: Define the global runtime inputs**

Use these workflow inputs:

```wdl
String docker_image = "ghcr.io/aou-multiomics-analysis/celltypedeconvolution:latest"
Int preemptible_attempts = 2
Int max_retries = 2
```

Pass the same values to every task call. Remove all task-prefixed preemptible and retry inputs from the workflow interface. Keep CPU, memory, and disk inputs unchanged.

- [ ] **Step 5: Convert all active WDL files to 1.0**

Set the first line of the main workflow and the four remaining task files to:

```wdl
version 1.0
```

Run MiniWDL after the change. If MiniWDL identifies 1.1-only syntax, replace only that syntax with its WDL 1.0 equivalent and preserve the task behavior.

- [ ] **Step 6: Run WDL tests and validation**

Run:

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-wdl-contract.R"); testthat::test_file("tests/testthat/test-wdl-manifest-boundary.R"); testthat::test_file("tests/testthat/test-workflow-validation.R")'
miniwdl check workflows/cell_type_deconvolution.wdl
Rscript tools/check_wdl_logging.R
```

Expected: PASS.

- [ ] **Step 7: Commit the WDL 1.0 workflow**

```bash
git add workflows scripts R/expression.R R/qc.R tests/testthat/test-wdl-contract.R tests/testthat/test-wdl-manifest-boundary.R tests/testthat/test-bed-outputs.R
git commit -m "refactor: simplify the WDL 1.0 interface"
```

---

### Task 6: Update smoke fixtures and correct the inventory assertion

**Files:**
- Modify: `scripts/generate_synthetic_fixture.R`
- Modify: `tests/fixtures/dtangle.inputs.json`
- Modify: `tests/fixtures/restart.inputs.json`
- Modify: `tests/smoke/assert_outputs.R`
- Modify: `tests/testthat/test-wdl-contract.R`
- Modify: `tests/testthat/test-documentation.R`

**Interfaces:**
- Consumes: the simplified workflow input schema
- Produces: deterministic dtangle and precomputed smoke JSON files and value-based output assertions

- [ ] **Step 1: Write failing fixture-contract tests**

Require both fixture JSON files to omit all preparation, pipeline-version, task-specific preemptible, and task-specific retry inputs. Require each to contain only the global controls:

```r
testthat::expect_identical(inputs$cell_type_deconvolution.preemptible_attempts, 2L)
testthat::expect_identical(inputs$cell_type_deconvolution.max_retries, 2L)
testthat::expect_false("cell_type_deconvolution.pipeline_version" %in% names(inputs))
```

Change the smoke-contract string from `identical(public_inventory, inventory)` to `all.equal()` with parser attributes ignored.

- [ ] **Step 2: Run fixture tests and verify RED**

Run:

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-wdl-contract.R"); testthat::test_file("tests/testthat/test-documentation.R")'
```

Expected: FAIL because fixtures still use the old public inputs and the smoke test compares readr metadata.

- [ ] **Step 3: Update fixture generation and checked-in JSON**

Generate these common runtime entries:

```r
"cell_type_deconvolution.preemptible_attempts" = 2L
"cell_type_deconvolution.max_retries" = 2L
```

Keep `celltype-deconvolution:test` in smoke fixtures because GitHub Actions builds that local CI tag. Remove `pipeline_version` and every task-prefixed preemptible or retry key.

- [ ] **Step 4: Correct the output-inventory comparison**

Replace parser-metadata identity with schema and value checks:

```r
stopifnot(
  identical(names(public_inventory), names(inventory)),
  isTRUE(all.equal(
    as.data.frame(public_inventory),
    as.data.frame(inventory),
    check.attributes = FALSE
  ))
)
```

This compares the published values and ignores source-specific readr parser attributes.

- [ ] **Step 5: Regenerate fixtures and run tests**

Run:

```bash
Rscript scripts/generate_synthetic_fixture.R tests/fixtures
Rscript scripts/generate_synthetic_fixture.R tests/fixtures
Rscript -e 'testthat::test_file("tests/testthat/test-wdl-contract.R"); testthat::test_file("tests/testthat/test-documentation.R")'
```

Expected: both generator runs succeed, the focused tests pass, and the deterministic fixture test confirms identical output from independent generations.

- [ ] **Step 6: Commit fixtures and smoke correction**

```bash
git add scripts/generate_synthetic_fixture.R tests/fixtures tests/smoke/assert_outputs.R tests/testthat/test-wdl-contract.R tests/testthat/test-documentation.R
git commit -m "test: update direct CPM workflow smoke inputs"
```

---

### Task 7: Update user documentation and examples

**Files:**
- Modify: `README.md`
- Modify: `docs/terra.md`
- Modify: `docs/data-dictionary.md`
- Modify: `examples/bed.inputs.json`
- Modify: `examples/precomputed-proportions.inputs.json`
- Modify: `tests/testthat/test-documentation.R`

**Interfaces:**
- Consumes: the final public WDL input and output schema
- Produces: Terra-ready examples and documentation that match the workflow

- [ ] **Step 1: Write failing documentation tests**

Require the examples to use the default image implicitly, omit `pipeline_version`, and expose only global retry controls:

```r
testthat::expect_false("cell_type_deconvolution.docker_image" %in% names(inputs))
testthat::expect_false("cell_type_deconvolution.pipeline_version" %in% names(inputs))
testthat::expect_true(all(c(
  "cell_type_deconvolution.preemptible_attempts",
  "cell_type_deconvolution.max_retries"
) %in% names(inputs)))
```

- [ ] **Step 2: Run the documentation test and verify RED**

Run:

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-documentation.R")'
```

Expected: FAIL because docs and examples describe task-specific controls and a required image.

- [ ] **Step 3: Update the documentation**

Document these facts in direct language:

- The BED is already normalized linear CPM.
- The workflow does not have an expression-preparation task.
- dtangle maps gene IDs with the GTF and applies `log2()` once.
- TCA reads every valid nonconstant gene from the BED and applies `log2()` once.
- The GitHub `latest` image is the default.
- Preemptible attempts and retries are global.
- All WDL sources use WDL 1.0.
- The removed preparation and pipeline-version outputs are no longer available.

Remove the Docker image field from user examples so the workflow default applies. Keep cloud paths and scientific settings in the examples.

- [ ] **Step 4: Run documentation and JSON checks**

Run:

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-documentation.R")'
python3 -m json.tool examples/bed.inputs.json
python3 -m json.tool examples/precomputed-proportions.inputs.json
```

Expected: PASS.

- [ ] **Step 5: Commit documentation**

```bash
git add README.md docs/terra.md docs/data-dictionary.md examples tests/testthat/test-documentation.R
git commit -m "docs: describe the direct CPM WDL 1.0 workflow"
```

---

### Task 8: Run final verification and update the pull request

**Files:**
- Verify: all changed source, workflow, test, fixture, example, and documentation files
- Modify only if a verification failure identifies a defect in the approved design

**Interfaces:**
- Consumes: completed Tasks 1 through 7
- Produces: a reviewed branch with local static evidence and GitHub Actions smoke evidence

- [ ] **Step 1: Run the full R test suite**

```bash
Rscript tests/testthat.R
```

Expected: zero failures. Dependency-based skips are acceptable only when the output names the unavailable package.

- [ ] **Step 2: Run R lint**

```bash
Rscript tools/lint_r.R
```

Expected: zero lint failures.

- [ ] **Step 3: Validate WDL and logging**

```bash
miniwdl check workflows/cell_type_deconvolution.wdl
Rscript tools/check_wdl_logging.R
```

Expected: both commands exit zero and every active WDL begins with `version 1.0`.

- [ ] **Step 4: Run repository integrity checks**

```bash
git diff --check
rg -n "PrepareExpression|expression_tasks|pipeline_version|version 1[.]1" workflows scripts R README.md docs/terra.md docs/data-dictionary.md examples
```

Expected: `git diff --check` exits zero. The active-source scan returns no obsolete contract references. Historical specifications and plans are excluded from this scan.

- [ ] **Step 5: Request code review**

Review the complete implementation range against the approved specification. Resolve critical and important findings with focused red-green tests. Record minor findings that do not block the requested workflow.

- [ ] **Step 6: Push the branch and monitor GitHub Actions**

```bash
git push origin feat/dtangle-tca-pipeline
gh run list --branch feat/dtangle-tca-pipeline --limit 5
```

Expected: the branch push succeeds and Docker Image CI starts for the existing pull request.

- [ ] **Step 7: Verify GitHub smoke workflows**

Resolve the GitHub run identifier from Step 6 and wait for it:

```bash
run_id=$(gh run list --branch feat/dtangle-tca-pipeline --limit 1 --json databaseId --jq '.[0].databaseId')
gh run watch "$run_id" --exit-status
```

Expected: the image build, R tests, lint, WDL checks, dtangle smoke workflow, precomputed smoke workflow, and output assertions pass.

- [ ] **Step 8: Commit any verification-only correction**

If Step 1 through Step 7 required a correction, run the affected focused test and the full verification set again, then commit the exact corrected files:

```bash
git add -u
git commit -m "fix: resolve direct CPM workflow verification"
```

If no correction was required, do not create an empty commit.
