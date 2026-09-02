# Cell type deconvolution for human whole blood

<!-- workflow-badges:start -->
[![Docker Image CI](https://github.com/AoU-Multiomics-Analysis/CellTypeDeconvolution/actions/workflows/docker-image.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/CellTypeDeconvolution/actions/workflows/docker-image.yml)
[![Pipeline CI](https://github.com/AoU-Multiomics-Analysis/CellTypeDeconvolution/actions/workflows/pipeline-ci.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/CellTypeDeconvolution/actions/workflows/pipeline-ci.yml)
[![Update README workflow badges](https://github.com/AoU-Multiomics-Analysis/CellTypeDeconvolution/actions/workflows/update-readme-badges.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/CellTypeDeconvolution/actions/workflows/update-readme-badges.yml)
<!-- workflow-badges:end -->


This Terra workflow estimates relative immune-cell proportions with dtangle. It then uses Tensor Composition Analysis (TCA) to estimate a cell-group-specific expression matrix for each retained group. It is for human whole-blood RNA-seq cohorts. The default resources support a cohort of about 9,000 samples.

## Input data

Provide a gene-by-sample normalized CPM matrix. The first column must be `gene_id`. Values must be finite and strictly positive. The workflow does not calculate, renormalize, or add a pseudocount to CPM values. A zero or negative CPM value stops expression preparation.

Provide a GTF or GTF.GZ file. The workflow uses `gene` records with `gene_id` and `gene_name`. It does not filter genes by `gene_type` or `gene_biotype`. It matches identifiers by exact text after it trims surrounding white space. It removes expression rows with no usable `gene_name` and sums CPM values for duplicated gene names.

Select one proportion mode:

- **dtangle mode:** Provide a user-supplied, standard linear LM22 signature matrix. The repository does not distribute LM22. The workflow requires the 22 standard LM22 columns and positive finite values. It transforms the matrix with `log2()` once.
- **Restart mode:** Provide a precomputed 22-type LM22 proportion table. It must have a sample identifier column and every standard LM22 cell-type column.

Both modes require `expression`, `gtf`, and `docker_image`. See [the Terra guide](docs/terra.md) and the two input examples in [examples](examples).

## Method

The workflow creates `log2(CPM)` without a pseudocount. In dtangle mode, it intersects bulk genes with LM22 genes in LM22 order. dtangle uses RNA-seq mixtures, `marker_method = "ratio"`, and a default marker fraction of `0.10`. Joint quantile normalization is off by default. When enabled, it applies only to the joined dtangle reference and mixture matrix. It never changes the TCA matrix.

The workflow combines the 22 LM22 estimates into ten groups: B cells, CD4 T cells, CD8 T cells, Gamma-delta T cells, NK cells, Monocyte/myeloid, Neutrophils, Eosinophils, Dendritic cells, and Mast cells. T follicular helper and regulatory T cells are part of CD4 T cells. Gamma-delta T cells remain separate. It removes any group with a cohort mean below `0.0001` (0.01%). No group is required to remain.

TCA fits one model across the cohort. It uses the complete mapped `log2(CPM)` matrix after removal of constant genes. It does not limit TCA to LM22-shared genes. It uses the processed dtangle weights with `refit_W = FALSE`.

## Outputs

The primary outputs are the original 22-type proportions, combined proportions, TCA weights, cell-group-specific HDF5 matrices, quality control, and an output manifest. Each matrix has scale `log2_cpm`. It is an inferred log-expression value. It is not a count, linear CPM, or absolute cell-type expression measurement.

Each group HDF5 file contains:

- `expression`: genes in rows and samples in columns.
- `gene_id`: gene identifiers in matrix order.
- `sample_id`: sample identifiers in matrix order.
- attributes for the cell-group name, `log2_cpm` scale, pipeline version, and model version.

Set `write_tsv` to `true` only when you need compressed TSV copies. HDF5 is the recommended output for large cohorts. [The data dictionary](docs/data-dictionary.md) lists every workflow input and output.

## Run and validation

Import `workflows/cell_type_deconvolution.wdl` into a Terra workspace. Use one of the JSON examples as a starting point. Replace the example cloud paths and select a published image tag for a production run. The local image tag `celltype-deconvolution:test` is only for the repository smoke fixtures.

GitHub Actions builds the project image and runs lint checks, R tests, WDL validation, and small dtangle and restart smoke workflows. The smoke workflows check the declared output structure. Do not require a local Docker build for validation.

## LM22 license

LM22 is user-supplied. Review and comply with the LM22 source license and terms before you download, store, or use it. This repository does not redistribute the signature matrix.

## Known limitations

- LM22 estimates relative immune proportions. It does not estimate absolute blood-cell counts.
- LM22 was derived from microarray data. Platform differences can affect estimates from RNA-seq mixtures.
- LM22 does not model erythrocyte or platelet expression.
- The Monocyte/myeloid group includes macrophage states. It is not a pure monocyte estimate.
- The default path does not use quantile normalization. Use it only as a sensitivity analysis.
- TCA outputs are statistical estimates. They are not directly measured cell-type-specific expression values.
- CPM zeros stop preparation. The workflow does not add a pseudocount.
- Identifier matching is exact after trimming. Different identifier systems can reduce gene overlap.
- Retained groups can differ by cohort because the threshold uses the cohort mean.
- Very small proportions can be unstable, even after the retained-group zero floor is applied.
