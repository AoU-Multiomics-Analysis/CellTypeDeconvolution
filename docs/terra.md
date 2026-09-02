# Run the workflow in Terra

Import `workflows/cell_type_deconvolution.wdl` into a Terra workspace. Use [cpm.inputs.json](../examples/cpm.inputs.json) for dtangle mode. Use [precomputed-proportions.inputs.json](../examples/precomputed-proportions.inputs.json) for restart mode. Replace every `gs://YOUR_BUCKET/...` path with a path that your Terra billing project can read.

## Required data

`expression` is a gene-by-sample normalized CPM matrix in TSV or compressed TSV format. Its first column is `gene_id`. Every CPM value must be finite and strictly positive. The workflow stops on zero or negative values. It does not calculate CPM, renormalize CPM, or add a pseudocount.

`gtf` is a GTF or GTF.GZ file. It must have `gene` records with `gene_id` and `gene_name` attributes. The workflow does not filter on `gene_type` or `gene_biotype`. It matches expression and GTF identifiers by exact text after trimming surrounding white space. It removes rows with no usable gene name and sums CPM values for duplicate gene names.

## Select a proportion mode

Use dtangle mode when you set `lm22` and omit `precomputed_proportions`. LM22 is user-supplied. It must be the standard linear LM22 matrix with one gene-symbol column and the 22 standard LM22 cell-type columns. Values must be finite and strictly positive. The workflow applies `log2()` once. It does not infer the LM22 scale.

Use restart mode when you set `precomputed_proportions` and omit `lm22`. The file must contain one sample identifier column and all 22 standard LM22 proportion columns. Its sample IDs must match the prepared expression matrix.

Both modes require `expression`, `gtf`, and `docker_image`. Use a published GHCR image tag in Terra. `celltype-deconvolution:test` is the local test image for repository smoke runs. It is not a production Terra image.

## Analysis settings

The workflow creates `log2(CPM)` without a pseudocount. In dtangle mode, it intersects the prepared bulk data with LM22 genes and preserves LM22 gene order. dtangle uses RNA-seq mixture data, `marker_method = "ratio"`, and a default marker fraction of `0.10`. The minimum LM22 overlap is `0.80`.

Joint quantile normalization is off by default. If you set `dtangle_quantile_normalize` to `true`, it applies once to the joined log-scale LM22 and bulk profiles. It applies only to the dtangle input. It does not change the expression matrix used by TCA.

The workflow combines the 22 LM22 types into ten groups:

| Group | LM22 members |
| --- | --- |
| B cells | B cells naive; B cells memory; Plasma cells |
| CD4 T cells | T cells CD4 naive; T cells CD4 memory resting; T cells CD4 memory activated; T cells follicular helper; T cells regulatory (Tregs) |
| CD8 T cells | T cells CD8 |
| Gamma-delta T cells | T cells gamma delta |
| NK cells | NK cells resting; NK cells activated |
| Monocyte/myeloid | Monocytes; Macrophages M0; Macrophages M1; Macrophages M2 |
| Neutrophils | Neutrophils |
| Eosinophils | Eosinophils |
| Dendritic cells | Dendritic cells resting; Dendritic cells activated |
| Mast cells | Mast cells resting; Mast cells activated |

T follicular helper and regulatory T cells are in CD4 T cells. Gamma-delta T cells remain separate. The workflow removes a group when its cohort mean is below `0.0001` (0.01%). No group is forced to remain. It replaces retained exact zeros with the configurable `zero_floor`, then normalizes each sample row to one.

TCA fits one cohort-wide model with `refit_W = FALSE`. It uses the full mapped `log2(CPM)` expression matrix after removal of constant genes. It does not use only LM22-shared genes. Optional covariates are mixture-level technical covariates.

## Resources and outputs

The default resources support an initial whole-blood run of about 9,000 samples. The TCA fitting task defaults to 16 CPU, 192 GB memory, and 750 GB disk. The assembly task defaults to 8 CPU, 128 GB memory, and 500 GB disk. The default TCA shard size is 500 genes. Adjust runtime inputs only after you inspect your cohort size and storage requirements.

Each retained group has one HDF5 matrix. The `expression` dataset has genes in rows and samples in columns. The file also contains `gene_name` and `sample_id` datasets. Attributes record the cell-group name, the `log2_cpm` scale, the pipeline version, and the TCA version. These matrices contain inferred log-expression estimates. They are not counts, linear CPM, or absolute expression values. Set `write_tsv` to `true` only when compressed TSV copies are required.

The workflow also writes proportions, filtering reports, the fitted TCA model, gene-shard records, reconstruction quality control, logs, an output inventory, and a machine-readable manifest. See [the data dictionary](data-dictionary.md) for exact names.

## Validation and support

GitHub Actions builds the pinned project image. It runs R lint and tests, validates WDL and command logging, and runs dtangle and restart smoke workflows. The smoke workflows assert the declared output structure. The project does not require a local Docker build for a smoke test.

LM22 is not distributed with this repository. Review and comply with the LM22 source license and terms before use.

## Known limitations

- LM22 estimates relative immune proportions. It does not estimate absolute blood-cell counts.
- LM22 was derived from microarray data. RNA-seq and reference platform differences can affect estimates.
- LM22 does not model erythrocyte or platelet expression.
- The Monocyte/myeloid group includes macrophage states. It is not a pure monocyte estimate.
- The default path has no quantile normalization. Use it only as a sensitivity analysis.
- TCA outputs are statistical estimates. They are not directly measured cell-type-specific expression values.
- CPM zeros stop preparation. The workflow does not add a pseudocount.
- Identifier matching is exact after trimming. Different identifier systems can reduce overlap.
- Retained groups can differ by cohort because the mean filter is cohort-specific.
- Very small proportions can be unstable, even after the zero floor is applied.
