qc_wdl_text <- function() {
  paste(
    readLines(testthat::test_path("../..", "workflows", "tasks", "qc.wdl")),
    collapse = "\n"
  )
}

testthat::test_that("AssembleTca exposes manifest inventory and QC plots", {
  text <- qc_wdl_text()

  testthat::expect_match(
    text, "File output_inventory = \"outputs/assembled_outputs.tsv\"",
    fixed = TRUE
  )
  testthat::expect_match(
    text, "File qc_plots = \"outputs/qc_plots.pdf\"", fixed = TRUE
  )
})

testthat::test_that("BuildManifest rewrites inventory paths to localized outputs", {
  text <- qc_wdl_text()

  purrr::walk(c(
    "localized_group_files", "stage_localized_file()", "basename",
    "ln -s \"$source_absolute\" \"$destination\"",
    "~{sep('\\n', group_hdf5)}", "~{sep('\\n', group_tsv)}",
    "rewritten_inventory=outputs/assembled_outputs.localized.tsv",
    "--outputs \"$rewritten_inventory\"",
    "if [[ ! -e \"$rewritten_path\" ]]",
    "File provenance = \"outputs/assembled_outputs.localized.tsv\""
  ), ~ testthat::expect_match(text, .x, fixed = TRUE))
})
