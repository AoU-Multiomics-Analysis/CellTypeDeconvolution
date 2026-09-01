# Whole-Blood CIBERSORTx and TCA Pipeline Design

Date: 2026-09-01

Status: Approved design

## Purpose

This project will provide a reproducible Terra workflow for human whole-blood RNA-seq data. The workflow will use CIBERSORTx Fractions with a user-supplied LM22 signature matrix to estimate cell-type proportions. It will combine the LM22 cell types into major groups. It will then use Tensor Composition Analysis (TCA) to estimate a gene-by-sample expression matrix for each retained group.

The target cohort has approximately 9,000 samples.

## Goals

The workflow will:

- Run on Terra with WDL.
- Accept raw counts or TPM as input.
- Accept the CIBERSORTx email and token at run time.
- Accept LM22 at run time. The repository will not distribute LM22.
- Estimate LM22 relative proportions with CIBERSORTx Fractions.
- Combine the 22 LM22 cell types into major whole-blood groups.
- Remove groups with a cohort mean proportion below 0.01%.
- Fit one cohort-wide TCA model.
- Produce one sample-specific gene-expression matrix for each retained group.
- Support a restart from precomputed fractions.
- Record validation, quality-control, and provenance results.

## Non-goals

The first version will not:

- Run CIBERSORTx HiRes.
- Estimate absolute blood-cell counts.
- Model erythrocytes or platelets unless a future signature includes them.
- Require CD4, CD8, or any other group to remain after filtering.
- Run source-specific association tests with `tcareg()`.
- Publish CIBERSORTx credentials, LM22, or licensed CIBERSORTx assets.

## Selected Approach

Three approaches were considered.

1. CIBERSORTx Fractions only. This approach estimates composition, but it does not produce sample-specific cell-type expression.
2. CIBERSORTx HiRes. This approach can estimate cell-type expression directly, but it creates a larger licensed workflow and was not selected for the first version.
3. CIBERSORTx Fractions followed by TCA. This approach separates proportion estimation from cell-type expression inference. It also permits a restart from saved fractions.

The workflow will use approach 3.

## Scientific Model

CIBERSORTx will receive linear TPM values. Quantile normalization will be off because the mixture data are RNA-seq. B-mode batch correction will be on by default and will be configurable.

TCA will receive `log2(TPM + 1)` values. Therefore, the TCA results are inferred log-expression values. They are not counts, raw TPM, or absolute expression measurements.

For each gene and sample, TCA models the bulk value as a weighted mixture of cell-type-specific values. The sample weights are the combined CIBERSORTx proportions. The first version will set `refit_W = FALSE`, so TCA will not change these proportions.

When the user supplies technical covariates, the workflow will pass them to TCA as mixture-level `C2` covariates. The first version will not use source-specific `C1` covariates.

TCA requires all retained weights to be positive and each sample row to sum to one. After the cohort filter, the workflow will replace exact zero values in retained groups with `1e-6`. It will then normalize each sample row to one. The zero floor will be a configurable input.

## Inputs

### Required inputs

- A gene-by-sample expression matrix in TSV or compressed TSV format.
- `expression_type`, with value `counts` or `tpm`.

The user must also select one fraction-input mode:

1. CIBERSORTx mode requires a user-supplied LM22 signature matrix, a CIBERSORTx account email, and a CIBERSORTx token file.
2. Restart mode requires a precomputed CIBERSORTx fraction table with one sample identifier column and all 22 LM22 fraction columns.

For count input, a gene annotation table is also required. It must contain:

- `gene_id`
- `gene_symbol`
- `gene_length_bp`

For TPM input, a gene annotation table is optional when the matrix already uses unique gene symbols.

### Optional inputs

- A sample-by-covariate table for TCA technical covariates.
- CIBERSORTx permutation count.
- B-mode selection.
- The minimum LM22 gene-overlap fraction. The default is `0.80`.
- The cell-group mean threshold. The default is `0.0001`, which equals 0.01%.
- The exact-zero floor. The default is `1e-6`.
- The number of genes in each tensor-extraction shard. The initial default is 500.
- TCA and WDL run-time settings, such as CPU, memory, disk, and retry values.
- A switch that requests compressed TSV copies of the primary HDF5 matrices.

## Expression Preparation

The expression matrix will use genes as rows and samples as columns. The workflow will preserve the sample order after it validates the identifiers.

For count input, the workflow will:

1. Divide each gene count by its length in kilobases.
2. Normalize the rates in each sample to one million.
3. Map gene identifiers to gene symbols.
4. Sum TPM values when more than one gene identifier maps to the same symbol.
5. Remove unmapped genes and normalize the remaining values in each sample to one million.

For TPM input, the workflow will validate nonnegative values. If annotation is supplied, it will map identifiers to symbols and sum duplicate symbols. After symbol resolution, it will normalize each sample to one million.

The workflow will make two matrices:

- Linear TPM for CIBERSORTx.
- `log2(TPM + 1)` for TCA.

Genes that are constant across all samples will not enter TCA. The workflow will record each removed gene and its removal reason.

## LM22 Grouping

The workflow will combine LM22 proportions as follows:

| Output group | LM22 members |
| --- | --- |
| B cells | B cells naive; B cells memory; Plasma cells |
| CD4 T cells | T cells CD4 naive; T cells CD4 memory resting; T cells CD4 memory activated; T cells follicular helper; T cells regulatory (Tregs) |
| CD8 T cells | T cells CD8 |
| NK cells | NK cells resting; NK cells activated |
| Monocyte/myeloid | Monocytes; Macrophages M0; Macrophages M1; Macrophages M2 |
| Neutrophils | Neutrophils |
| Eosinophils | Eosinophils |
| Dendritic cells | Dendritic cells resting; Dendritic cells activated |
| Other LM22 | T cells gamma delta; Mast cells resting; Mast cells activated |

The `Other LM22` name is intentional. This group contains only the listed LM22 members. It does not represent all cell types that LM22 does not model.

After grouping, the workflow will calculate the cohort mean for each group. It will remove a group when its mean is less than `0.0001`. This rule applies to all groups. No group is mandatory.

The workflow will preserve three fraction tables:

1. The original 22-type CIBERSORTx table.
2. The combined table before filtering and zero adjustment.
3. The filtered, adjusted, and normalized table that TCA uses.

## Terra Workflow Architecture

The repository will have one top-level workflow:

`workflows/cell_type_deconvolution.wdl`

The top-level workflow will import separate task definitions for expression preparation, CIBERSORTx, TCA, matrix assembly, and QC.

The workflow will run these stages:

1. Validate and prepare expression data.
2. Run CIBERSORTx Fractions, or accept precomputed fractions.
3. Combine and filter the LM22 proportions.
4. Fit one TCA model across all samples and all eligible genes.
5. Divide the eligible genes into deterministic shards.
6. Extract the TCA tensor for each gene shard.
7. Assemble one matrix for each retained group.
8. Calculate QC results and create the output manifest.

The workflow will not divide the cohort by sample. A sample split would fit different models to different parts of the cohort and would make the inferred matrices less comparable.

The global TCA fit will estimate one common noise parameter and gene-specific source parameters. Each tensor-extraction task will use `tcasub()` to select the parameters for its genes. It will then use `tensor()` with the same cohort-wide model. This design limits tensor memory use without changing the fitted cohort model.

All shards will use a stable gene order and a fixed random seed when a package operation uses randomness. A failed or missing shard will stop matrix assembly.

## Containers

The workflow will use two container sources:

- The official CIBERSORTx Fractions container for proportion estimation.
- A project container for expression preparation, TCA, assembly, and QC.

The project container will use micromamba as its base. It will install pinned packages from `conda-forge` and `bioconda` where packages are available. It will include pinned versions of R, TCA, tidyverse packages, HDF5 support, and the test dependencies.

The repository will record container tags and immutable image identifiers in the output provenance. The workflow will never copy the CIBERSORTx token into an output.

## Outputs

### Primary outputs

- Original LM22 proportions.
- Combined proportions before filtering.
- Filtered and adjusted proportions used by TCA.
- One HDF5 expression matrix for each retained group.
- A machine-readable output manifest.
- A QC summary.

Each HDF5 output will contain:

- An `expression` dataset with genes as rows and samples as columns.
- Gene identifiers.
- Sample identifiers.
- The cell-group name.
- The scale value `log2_tpm_plus_1`.
- Pipeline and model version attributes.

Compressed TSV matrices will be optional because they can use substantially more storage and transfer time than HDF5 at this cohort size.

### Supporting outputs

- Prepared linear TPM.
- Prepared `log2(TPM + 1)`.
- Gene mapping and duplicate-resolution reports.
- Excluded-gene report.
- Cell-group filtering report.
- TCA model object and model log.
- Per-task logs.
- Checksums and software version records.

## Validation and Failure Rules

The workflow will stop for these conditions:

- Duplicate sample identifiers.
- Missing, negative, or nonnumeric expression values.
- Missing or invalid gene lengths for count input.
- Sample mismatches between expression, fractions, and covariates.
- Missing LM22 columns that prevent the approved group mapping.
- An LM22 gene-overlap fraction below the configurable threshold. The default threshold is `0.80`.
- Negative or nonfinite fraction values.
- Invalid fraction row sums after normalization.
- Fewer than two groups after the cohort-mean filter.
- Nonpositive values in the final TCA weight matrix.
- Nonfinite TCA estimates.
- Missing, duplicated, or incomplete gene shards.
- Output matrix dimensions or identifiers that do not match the manifest.

The workflow will not use a fixed reconstruction threshold as an automatic failure condition in the first version. Suitable thresholds can depend on the cohort.

## Quality Control

The QC results will include:

- Number and fraction of LM22 genes found in the mixture matrix.
- Duplicate-gene counts and resolution actions.
- CIBERSORTx correlation, RMSE, and empirical p-value when available.
- Cohort mean, retention status, and filter reason for each combined group.
- Number of exact zero values adjusted for each retained group.
- Maximum fraction row-sum error before and after normalization.
- TCA convergence information and fitted model parameters.
- Counts of excluded and processed genes.
- Matrix dimensions for each retained group.
- Per-sample bulk reconstruction correlation and RMSE in log-expression space.
- Summary distributions of reconstruction metrics.
- Software versions, container identifiers, parameters, and file checksums.

The reconstructed bulk value for QC will be the proportion-weighted sum of the inferred cell-type values. When technical covariates are supplied, the reconstruction will also include the fitted mixture-level covariate term. This comparison will use the same `log2(TPM + 1)` scale as the TCA input.

Plots will use a clean, minimal style. They will not contain titles or subtitles.

## Logging and Credential Handling

Every WDL command block will log:

- Stage name.
- Start time.
- Completion time.
- Input dimensions or input file names when safe.
- Output paths.
- Clear error context.

Commands will not enable shell tracing when they handle credentials. The CIBERSORTx token will be a localized input file. The command will read the token at run time and will not print it. The token file and its contents will not be workflow outputs.

The CIBERSORTx email is a normal Terra string input. Documentation will tell users that Terra metadata can record normal string inputs.

## Test Strategy

### R tests

Unit tests will cover:

- Count-to-TPM conversion.
- Gene mapping and duplicate handling.
- LM22 group aggregation.
- The cohort-mean group filter.
- Zero adjustment and row normalization.
- Sample-order validation.
- Gene-shard construction and assembly.
- HDF5 dimensions, identifiers, and attributes.
- Reconstruction QC calculations.

R code will use tidyverse syntax where it improves clarity.

### Workflow tests

A small deterministic synthetic data set will exercise the complete open part of the workflow. It will use precomputed fractions, so pull-request tests will not require CIBERSORTx credentials. The test will verify:

- Successful WDL execution.
- Expected retained and removed groups.
- Expected matrix dimensions and identifiers.
- Positive normalized TCA weights.
- Complete shard assembly.
- Finite tensor and QC results.

A separate manual integration test can run CIBERSORTx with repository secrets. It will not run for untrusted pull requests.

### GitHub Actions

GitHub Actions will:

- Validate WDL files.
- Lint and test R scripts.
- Build the project container.
- Run the synthetic workflow smoke test in the built container.
- Push versioned project images to GitHub Container Registry from the main branch.

A local Docker build will not be required for the standard development smoke test.

## Terra and Dockstore Delivery

The repository will contain:

- Dockstore metadata for the top-level WDL workflow.
- Example JSON input files for count, TPM, and precomputed-fraction modes.
- A Terra run guide.
- A data dictionary for all input and output tables.
- Guidance for setting CPU, memory, disk, retry, and preemptible values.

Initial run-time defaults will target a 9,000-sample cohort. Users will be able to override the defaults without editing WDL source files.

## Known Limitations

- LM22 estimates relative proportions among its modeled immune cell types. It does not give absolute cell counts.
- LM22 does not model all whole-blood components. In particular, the selected grouping does not add erythrocyte or platelet expression.
- The `Other LM22` group is not a complete unmodeled-cell category.
- The monocyte/myeloid group includes LM22 macrophage states, although circulating whole blood usually has few macrophages.
- TCA results are statistical estimates. They are not measurements from sorted cells.
- TCA will operate in log-expression space. Its results must not be labeled as TPM.
- Group retention is cohort-specific. Two cohorts can produce different sets of output matrices.
- Very small retained proportions can produce unstable estimates even after exact zero adjustment.
- CIBERSORTx use is subject to its current license and service terms.

## References

- [CIBERSORT](https://www.nature.com/articles/nmeth.3337)
- [CIBERSORTx](https://www.nature.com/articles/s41587-019-0114-2)
- [CIBERSORTx protocol](https://pmc.ncbi.nlm.nih.gov/articles/PMC7695353/)
- [TCA paper](https://www.nature.com/articles/s41467-019-11052-9)
- [TCA package](https://github.com/cozygene/TCA)
- [TCA function documentation](https://rdrr.io/cran/TCA/man/tca.html)
