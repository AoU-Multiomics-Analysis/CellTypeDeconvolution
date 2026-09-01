# Whole-Blood dtangle and TCA Pipeline Design

Date: 2026-09-01

Status: Approved design

## Purpose

This project will provide a reproducible Terra workflow for human whole-blood RNA-seq data. The workflow will use dtangle with a user-supplied LM22 signature matrix to estimate cell-type proportions. It will combine the 22 LM22 cell types into major groups. It will then use Tensor Composition Analysis (TCA) to estimate a gene-by-sample expression matrix for each retained group.

The target cohort has approximately 9,000 samples.

## Goals

The workflow will:

- Run on Terra with WDL.
- Accept raw counts or TPM as input.
- Accept LM22 at run time. The repository will not distribute LM22.
- Estimate 22 relative cell-type proportions with the open-source dtangle R package.
- Combine the 22 proportions into major whole-blood groups.
- Remove groups with a cohort mean proportion below 0.01%.
- Fit one cohort-wide TCA model.
- Produce one sample-specific gene-expression matrix for each retained group.
- Support a restart from precomputed proportions.
- Record validation, quality-control, and provenance results.

## Non-goals

The first version will not:

- Run CIBERSORT or CIBERSORTx.
- Require a CIBERSORT account, email, or token.
- Estimate absolute blood-cell counts.
- Model erythrocytes or platelets unless a future signature includes them.
- Require CD4, CD8, or any other group to remain after filtering.
- Run source-specific association tests with `tcareg()`.
- Redistribute LM22.

## Selected Approach

Three dtangle integration approaches were considered.

1. Run dtangle with all 22 LM22 profiles, then combine the estimated proportions. This preserves the full reference during marker selection and permits the approved major-group mapping.
2. Combine the LM22 reference profiles before dtangle. This reduces the number of source types, but it can blur distinct reference profiles and change marker selection.
3. Replace LM22 with an RNA-seq reference. This can reduce the platform difference, but it changes the cell-type definitions and requires a separate reference-validation project.

The workflow will use approach 1. It will run `dtangle()` at the original 22-type resolution and combine proportions after deconvolution.

## Scientific Model

### dtangle inputs and scale

The standard LM22 file contains non-log, linear expression values. It contains 547 gene rows and 22 cell-type columns. dtangle expects base-2 log-scale expression for both the mixtures and references.

The workflow will therefore use:

- `log2(TPM + 1)` for the bulk RNA-seq mixtures.
- `log2(LM22)` for the standard LM22 reference.

The workflow will not add a pseudocount to LM22 because the standard matrix contains positive values. It will stop if an LM22 value is zero or negative.

The workflow will not apply quantile normalization by default. The dtangle stage will use TPM normalization and the scale conversions above. Quantile normalization will be an optional sensitivity setting and will be off by default. When enabled, the workflow will join the log-scale reference and mixture profiles by their shared genes and use `limma::normalizeBetweenArrays()` once on the combined matrix. This setting will apply only to the dtangle input. It will never change the TCA input matrix.

The workflow will run `dtangle()` with:

- One LM22 reference profile for each of the 22 cell types.
- `data_type = "rna-seq"` because the mixture samples are RNA-seq.
- `marker_method = "ratio"`.
- The top 10% of ranked markers for each cell type.

The workflow will save the marker genes, number of markers, and fitted `gamma` value returned by dtangle.

### TCA model

TCA will receive the original `log2(TPM + 1)` matrix. Therefore, the TCA results are inferred log-expression values. They are not counts, raw TPM, or absolute expression measurements.

For each gene and sample, TCA models the bulk value as a weighted mixture of cell-type-specific values. The sample weights are the combined dtangle proportions. The first version will set `refit_W = FALSE`, so TCA will not change these proportions.

When the user supplies technical covariates, the workflow will pass them to TCA as mixture-level `C2` covariates. The first version will not use source-specific `C1` covariates.

TCA requires all retained weights to be positive and each sample row to sum to one. After the cohort filter, the workflow will replace exact zero values in retained groups with `1e-6`. It will then normalize each sample row to one. The zero floor will be a configurable input.

## Inputs

### Required inputs

- A gene-by-sample expression matrix in TSV or compressed TSV format.
- `expression_type`, with value `counts` or `tpm`.

The user must also select one proportion-input mode:

1. dtangle mode requires a user-supplied LM22 signature matrix.
2. Restart mode requires a precomputed 22-type dtangle proportion table with one sample identifier column and all 22 LM22 proportion columns.

For count input, a gene annotation table is also required. It must contain:

- `gene_id`
- `gene_symbol`
- `gene_length_bp`

For TPM input, a gene annotation table is optional when the matrix already uses unique gene symbols.

### Optional inputs

- A sample-by-covariate table for TCA technical covariates.
- The minimum LM22 gene-overlap fraction. The default is `0.80`.
- The dtangle marker fraction. The default is `0.10`.
- The dtangle marker-ranking method. The default and supported first-version value is `ratio`.
- A switch for joint quantile normalization of the dtangle reference and mixtures. The default is `false`.
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

The workflow will make one primary transformed matrix:

- `log2(TPM + 1)` for dtangle mixtures and TCA.

The workflow will keep linear TPM as a supporting output. It will create a separate dtangle working matrix that contains only genes shared with LM22. Optional quantile normalization, when enabled, will change only this working matrix.

Genes that are constant across all samples will not enter TCA. The workflow will record each removed gene and its removal reason.

## LM22 Validation and Preparation

The workflow will expect the standard LM22 structure:

- One gene-symbol column.
- All 22 standard LM22 cell-type columns.
- Unique gene symbols.
- Finite, positive expression values.

The workflow will not infer the LM22 scale from its numeric range. It will treat the input as the standard linear LM22 matrix and apply `log2()` once. This explicit rule prevents an ambiguous automatic scale decision.

The workflow will match LM22 and bulk genes by exact gene symbol after trimming surrounding white space. It will report missing and matched LM22 genes. It will stop when the matched fraction is below the configured minimum.

The dtangle task will transpose the matrices to the package interface:

- Mixtures: samples by shared genes.
- References: 22 cell types by shared genes.

## LM22 Grouping

The workflow will combine the 22 dtangle proportions as follows:

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

The workflow will preserve three proportion tables:

1. The original 22-type dtangle table.
2. The combined table before filtering and zero adjustment.
3. The filtered, adjusted, and normalized table that TCA uses.

## Terra Workflow Architecture

The repository will have one top-level workflow:

`workflows/cell_type_deconvolution.wdl`

The top-level workflow will import separate task definitions for expression preparation, dtangle, proportion processing, TCA, matrix assembly, and quality control.

The workflow will run these stages:

1. Validate and prepare expression data.
2. Validate and transform LM22, then run dtangle, or accept precomputed proportions.
3. Combine and filter the LM22 proportions.
4. Fit one TCA model across all samples and all eligible genes.
5. Divide the eligible genes into deterministic shards.
6. Extract the TCA tensor for each gene shard.
7. Assemble one matrix for each retained group.
8. Calculate quality-control results and create the output manifest.

The dtangle task will process the complete cohort in one task. It uses only the LM22-overlap genes, so the working matrix is small relative to the complete transcriptome.

The workflow will not divide the cohort by sample for TCA. A sample split would fit different models to different parts of the cohort and would make the inferred matrices less comparable.

The global TCA fit will estimate one common noise parameter and gene-specific source parameters. Each tensor-extraction task will use `tcasub()` to select the parameters for its genes. It will then use `tensor()` with the same cohort-wide model. This design limits tensor memory use without changing the fitted cohort model.

All shards will use a stable gene order and a fixed random seed when a package operation uses randomness. A failed or missing shard will stop matrix assembly.

## Container

The workflow will use one project container for expression preparation, dtangle, proportion processing, TCA, assembly, and quality control.

The project container will:

- Use micromamba as its base.
- Install pinned packages from `conda-forge` and `bioconda` where packages are available.
- Pin `r-dtangle` to version `2.0.10`.
- Pin R, TCA, tidyverse packages, HDF5 support, and test dependencies.
- Use an exact package version or source revision for any R package that is not available from the selected Conda channels.

The environment definition and container digest will be part of the output provenance. GitHub Actions, not a required local Docker build, will build and smoke-test the image.

## Outputs

### Primary outputs

- Original 22-type dtangle proportions.
- Combined proportions before filtering.
- Filtered and adjusted proportions used by TCA.
- One HDF5 expression matrix for each retained group.
- A machine-readable output manifest.
- A quality-control summary.

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
- The transformed LM22 matrix used by dtangle.
- The dtangle shared-gene matrix.
- The dtangle marker table and model metadata.
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
- Sample mismatches between expression, precomputed proportions, and covariates.
- Missing LM22 columns that prevent the approved 22-type model or group mapping.
- Duplicate LM22 gene symbols.
- Zero, negative, missing, or nonfinite LM22 expression values.
- An LM22 gene-overlap fraction below the configurable threshold. The default threshold is `0.80`.
- A cell type with no selected dtangle marker genes.
- Negative, nonfinite, or invalid dtangle proportion values.
- dtangle proportion row sums outside numeric tolerance.
- Fewer than two groups after the cohort-mean filter.
- Nonpositive values in the final TCA weight matrix.
- Nonfinite TCA estimates.
- Missing, duplicated, or incomplete gene shards.
- Output matrix dimensions or identifiers that do not match the manifest.

The workflow will not use a fixed reconstruction threshold as an automatic failure condition in the first version. Suitable thresholds can depend on the cohort.

## Quality Control

The quality-control results will include:

- Number and fraction of LM22 genes found in the mixture matrix.
- LM22 dimensions, value range, and validation status.
- Duplicate-gene counts and resolution actions.
- dtangle package version, `gamma`, marker method, marker fraction, and marker count by cell type.
- Maximum proportion row-sum error.
- Cohort mean, retention status, and filter reason for each combined group.
- Number of exact zero values adjusted for each retained group.
- Maximum fraction row-sum error before and after final normalization.
- TCA convergence information and fitted model parameters.
- Counts of excluded and processed genes.
- Matrix dimensions for each retained group.
- Per-sample bulk reconstruction correlation and RMSE in log-expression space.
- Summary distributions of reconstruction metrics.
- Software versions, container identifiers, parameters, and file checksums.

The reconstructed bulk value for TCA quality control will be the proportion-weighted sum of the inferred cell-type values. When technical covariates are supplied, the reconstruction will also include the fitted mixture-level covariate term. This comparison will use the same `log2(TPM + 1)` scale as the TCA input.

Plots will use a clean, minimal style. They will not contain titles or subtitles.

## Logging

Every WDL command block will log:

- Stage name.
- Start time.
- Completion time.
- Input dimensions or input file names when safe.
- Output paths.
- Clear error context.

The dtangle log will also record the shared-gene count, marker count for each cell type, package version, and selected preprocessing settings.

## Test Strategy

### R tests

Unit tests will cover:

- Count-to-TPM conversion.
- Gene mapping and duplicate handling.
- LM22 structure and scale validation.
- The `log2(LM22)` conversion.
- Bulk and LM22 gene matching.
- dtangle matrix orientation.
- Marker selection and proportion extraction.
- LM22 group aggregation.
- The cohort-mean group filter.
- Zero adjustment and row normalization.
- Sample-order validation.
- Gene-shard construction and assembly.
- HDF5 dimensions, identifiers, and attributes.
- Reconstruction quality-control calculations.

R code will use tidyverse syntax where it improves clarity.

### Workflow tests

A small deterministic synthetic data set will exercise the complete workflow. It will include a synthetic LM22-shaped reference with the 22 required column names. It will not contain or derive from the licensed LM22 values.

The smoke test will verify:

- Successful WDL execution.
- Finite, nonnegative 22-type dtangle proportions whose rows sum to one.
- Expected retained and removed groups.
- Expected matrix dimensions and identifiers.
- Positive normalized TCA weights.
- Complete shard assembly.
- Finite tensor and quality-control results.

A second test will start from precomputed proportions and verify the restart path.

### GitHub Actions

GitHub Actions will:

- Validate WDL files.
- Lint and test R scripts.
- Build the micromamba project container.
- Run the synthetic end-to-end workflow smoke test in the built container.
- Run the restart-path smoke test.
- Push versioned project images to GitHub Container Registry from the main branch.

A local Docker build will not be required for the standard development smoke test.

## Terra and Dockstore Delivery

The repository will contain:

- Dockstore metadata for the top-level WDL workflow.
- Example JSON input files for count, TPM, and precomputed-proportion modes.
- A Terra run guide.
- A data dictionary for all input and output tables.
- Guidance for setting CPU, memory, disk, retry, and preemptible values.

Initial run-time defaults will target a 9,000-sample cohort. Users will be able to override the defaults without editing WDL source files.

## Known Limitations

- LM22 estimates relative proportions among its modeled immune cell types. It does not give absolute cell counts.
- LM22 was developed from microarray data, while the target mixtures are RNA-seq. The log-scale conversion does not remove every cross-platform effect.
- The default workflow does not apply quantile normalization or another batch-correction method between LM22 and the mixtures.
- LM22 does not model all whole-blood components. In particular, the selected grouping does not add erythrocyte or platelet expression.
- The `Other LM22` group is not a complete unmodeled-cell category.
- The monocyte/myeloid group includes LM22 macrophage states, although circulating whole blood usually has few macrophages.
- dtangle estimates proportions only. TCA performs the separate cell-type-specific expression inference.
- TCA results are statistical estimates. They are not measurements from sorted cells.
- TCA will operate in log-expression space. Its results must not be labeled as TPM.
- Group retention is cohort-specific. Two cohorts can produce different sets of output matrices.
- Very small retained proportions can produce unstable estimates even after exact zero adjustment.

## References

- [dtangle paper](https://academic.oup.com/bioinformatics/article/35/12/2093/5165376)
- [dtangle source repository](https://github.com/gjhunt/dtangle)
- [dtangle package documentation](https://search.r-project.org/CRAN/refmans/dtangle/html/dtangle.html)
- [dtangle cross-platform vignette](https://gjhunt.github.io/dtangle/vign/sc_vignette.html)
- [LM22 and CIBERSORT paper](https://www.nature.com/articles/nmeth.3337)
- [CIBERSORT protocol with LM22 scale guidance](https://pmc.ncbi.nlm.nih.gov/articles/PMC5895181/)
- [TCA paper](https://www.nature.com/articles/s41467-019-11052-9)
- [TCA package](https://github.com/cozygene/TCA)
- [TCA function documentation](https://rdrr.io/cran/TCA/man/tca.html)
