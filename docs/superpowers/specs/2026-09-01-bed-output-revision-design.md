# Direct Cell-Type BED Output Design

**Date:** 2026-09-01

**Status:** Approved in conversation; written review pending

## Scope

This design revises the whole-blood dtangle and TCA pipeline. It supersedes the input-matrix, tensor-extraction, HDF5, and final-output parts of `2026-09-01-whole-blood-dtangle-tca-design.md`.

The workflow will accept a BED-like CPM matrix. It will preserve genomic coordinates. It will fit TCA once on the full prepared matrix. It will extract the full tensor once. It will write one compressed BED file for each retained major cell group.

The workflow will not create HDF5 files. It will not split genes into shards.

## Goals

- Preserve the genomic coordinates and gene IDs from the input expression BED.
- Keep original gene-ID rows for TCA.
- Use gene-symbol aggregation only for dtangle and LM22 alignment.
- Write one direct, analysis-ready BED.GZ matrix for each retained major group.
- Keep both the dtangle and precomputed-proportion workflow modes.
- Keep cohort-level TCA fitting, proportion outputs, exclusion reports, QC results, and provenance.
- Keep resource settings suitable as a starting profile for about 9,000 samples.

## Non-goals

- Do not create final or intermediate HDF5 files.
- Do not shard the tensor by gene.
- Do not write one file for each of the 22 LM22 subtypes.
- Do not recover cell-type values for genes that TCA cannot model.
- Do not change the ten-group LM22 mapping or the cohort-mean filter rule.

## Input Expression Contract

The expression input must be a tab-delimited BED or BED.GZ file. It must have these leading columns in this order:

```text
#chr    start    end    gene_id    sample_1    sample_2    ...
```

The sample columns contain normalized linear CPM values. Values must be finite and strictly positive. The workflow will not normalize values. It will not add a pseudocount.

The workflow will validate these rules:

- `#chr`, `start`, `end`, and `gene_id` exist in the required order.
- Chromosome values and gene IDs are not empty.
- Start and end values are valid integer coordinates.
- Each start value is less than its end value.
- Gene IDs are unique after surrounding white space is removed.
- Sample IDs are non-empty and unique.
- Sample values are finite and strictly positive.

The workflow will preserve chromosome, start, end, gene ID, gene order, and sample order in each final BED output. It will not convert coordinate systems.

## GTF Mapping and Two Expression Views

The GTF or GTF.GZ input remains required. The workflow will read `gene` records. It will map `gene_id` to `gene_name` by exact text after it removes surrounding white space. It will not filter on `gene_type` or `gene_biotype`.

Expression preparation will create two distinct views.

### Full TCA view

The TCA view will keep one row for each mapped input `gene_id`. It will preserve input BED order and coordinate metadata. It will not collapse rows that share one `gene_name`.

The workflow will calculate `log2(CPM)` without a pseudocount. It will remove constant genes before TCA fitting. It will record unmapped and constant genes in separate exclusion reports.

### dtangle view

The dtangle view will map rows to `gene_name`. It will sum linear CPM values when multiple input gene IDs map to the same gene name. It will then calculate `log2(CPM)` without a pseudocount.

Only this gene-symbol view will intersect with the user-supplied linear LM22 signature. LM22 order and the marker method `ratio` will remain unchanged.

This separation prevents gene-symbol aggregation from losing the coordinate identity that the TCA BED outputs need.

## Proportion Processing

The workflow will keep both proportion modes.

- In dtangle mode, it will estimate the standard 22 LM22 proportions.
- In restart mode, it will accept a precomputed sample-by-22 proportion table.

The workflow will combine the 22 estimates into these ten major groups:

- B cells
- CD4 T cells
- CD8 T cells
- Gamma-delta T cells
- NK cells
- Monocyte/myeloid
- Neutrophils
- Eosinophils
- Dendritic cells
- Mast cells

T follicular helper and regulatory T cells will remain in CD4 T cells. Gamma-delta T cells will remain separate.

The workflow will remove a group when its cohort mean is below `0.0001`. It will not require any named group to remain. It must retain at least two groups for TCA.

## TCA Fit and Direct Tensor Extraction

The workflow will fit one cohort-level TCA model. The model input will be the full mapped, nonconstant gene-ID matrix on the `log2(CPM)` scale. The weights will be the retained major-group proportions.

The workflow will extract the full cell-type tensor in one operation. It will not scatter extraction across gene shards. The implementation will assign tensor list names from `colnames(W)`, where `W` is the TCA weight matrix.

Before file output, it will validate each cell-group matrix:

- The list contains each retained cell group exactly once.
- Each matrix has genes in rows and samples in columns.
- Gene IDs match the prepared TCA matrix in exact order.
- Sample IDs match the input BED in exact order.
- Dimensions match the expected gene and sample counts.
- All inferred values are finite.

The inferred values remain on the `log2(CPM)` model scale. They are statistical cell-type-specific estimates. They are not linear CPM values or directly measured expression values.

## BED Output Contract

The workflow will write one `.bed.gz` file for each retained major group. File names will use stable group slugs. Examples include `cd4_t_cells.bed.gz` and `cd8_t_cells.bed.gz`.

Each file will have this schema:

```text
#chr    start    end    gene_id    sample_1    sample_2    ...
```

Rows will contain only genes that map through the GTF and remain nonconstant. Rows will follow their original input order. The first four fields will come from the input BED. Sample columns will follow the input order. Sample cells will contain the inferred values for that major group.

The workflow will also write a cell-group output table. Each row will contain the exact cell-group name, stable slug, BED file path, gene count, sample count, and scale `log2_cpm`.

## Supporting Outputs

The workflow will keep these supporting results:

- Prepared linear CPM and `log2(CPM)` matrices needed for method stages.
- GTF mapping and exclusion reports.
- Original 22-type proportions.
- Combined ten-group proportions.
- Adjusted TCA weights.
- Cell-group filter report.
- TCA model RDS.
- TCA expression and constant-gene exclusion report.
- Per-sample reconstruction metrics.
- QC summary and clean minimal QC plots without titles or subtitles.
- Effective parameters and software versions.
- A machine-readable output manifest with checksums.
- Task logs with stage, UTC start, completion, dimensions, and output paths.

The output manifest will list each cell-group BED as a primary output. It will record the group name, dimensions, scale, checksum, pipeline version, TCA version, and container image.

## Workflow Structure

The WDL will keep separate tasks for expression preparation, dtangle, proportion processing, TCA fitting, direct tensor export, and manifest creation.

The direct tensor-export task will load the fitted model and the full prepared TCA matrix. It will extract the complete tensor, calculate reconstruction QC, and write all retained cell-group BED files in one task.

The workflow will remove:

- Gene-shard creation and shard manifests.
- Tensor-shard extraction tasks.
- HDF5 assembly tasks.
- `hdf5r` from the pinned environment.
- HDF5 and optional-TSV workflow outputs.
- Shard size and shard resource inputs.
- HDF5 assembly settings and documentation.

The direct export task will have configurable CPU, memory, disk, retry, and preemptible settings. Existing high-memory assembly defaults can provide the starting values for the 9,000-sample profile.

## Error Handling

Every validation failure will stop the task with a nonzero exit status. Errors will identify the failed stage and rule. The workflow will stop for:

- A malformed input BED or invalid coordinate fields.
- Duplicate gene IDs or sample IDs.
- Nonfinite, zero, or negative CPM values.
- A malformed GTF or no mapped genes.
- Insufficient LM22 overlap in dtangle mode.
- A malformed precomputed proportion table in restart mode.
- Fewer than two retained major groups.
- A TCA fit failure.
- A tensor group, dimension, identifier, order, or finite-value mismatch.
- A BED output or manifest write failure.

## Testing

Unit tests will cover:

- BED and BED.GZ parsing.
- Exact preservation of coordinates, gene IDs, gene order, and sample order.
- Invalid coordinates, duplicate IDs, and invalid CPM values.
- Duplicate gene names that collapse in the dtangle view but remain separate in the TCA view.
- Unmapped and constant-gene exclusion.
- Tensor list naming from the TCA weight columns.
- Direct BED output schema, dimensions, values, and group file names.
- Reconstruction QC from the in-memory tensor.
- Output-manifest BED metadata.
- Removal of active HDF5 and sharding contracts.

GitHub Actions will build the pinned Micromamba image. It will run the full R test suite, R lint, WDL validation, logging checks, and two complete smoke workflows. One smoke will use dtangle mode. One smoke will use precomputed proportions. Both smoke runs will verify the BED coordinate and sample contracts.

A local Docker build will not be required.

## Documentation and Migration

README, Terra guidance, Dockstore examples, input examples, the data dictionary, and smoke fixtures will use the BED input contract.

The documentation will state that:

- The pipeline no longer accepts a gene-ID-first matrix without BED coordinates.
- The primary expression outputs are one BED.GZ file per retained major group.
- Output values use the `log2_cpm` model scale.
- The workflow does not produce HDF5 or tensor-shard files.
- The set of output groups can differ between cohorts.
- LM22 remains a user-supplied licensed input.

## Acceptance Criteria

- Both workflow modes accept the defined BED expression input.
- TCA keeps distinct input gene IDs even when gene symbols are duplicated.
- Each retained group has exactly one coordinate-preserving BED.GZ output.
- Every BED output preserves the original modeled-gene and sample order.
- BED values match the corresponding in-memory TCA matrix.
- HDF5 and sharding are absent from the active workflow, image dependency set, tests, and user documentation.
- Local static checks pass.
- The GitHub Actions image build and both smoke workflows pass before release.
