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
    zero_floor = 1e-6
  )
}

slugify_cell_group <- function(x) {
  x |>
    stringr::str_to_lower() |>
    stringr::str_replace_all("[^a-z0-9]+", "_") |>
    stringr::str_remove("^_") |>
    stringr::str_remove("_$")
}
