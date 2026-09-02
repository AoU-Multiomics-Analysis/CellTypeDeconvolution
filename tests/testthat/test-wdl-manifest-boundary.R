qc_wdl_text <- function() {
  paste(
    readLines(testthat::test_path("../..", "workflows", "tasks", "qc.wdl")),
    collapse = "\n"
  )
}

testthat::test_that("BuildManifest exposes direct BED inventory and QC plots", {
  text <- qc_wdl_text()

  testthat::expect_match(
    text, "File cell_type_bed_inventory", fixed = TRUE
  )
  testthat::expect_match(
    text, "File qc_summary = export_qc_summary", fixed = TRUE
  )
  testthat::expect_match(
    text, "File qc_plots = export_qc_plots", fixed = TRUE
  )
  testthat::expect_match(
    text, "assembled_outputs.localized.tsv", fixed = TRUE
  )
  testthat::expect_match(
    text, "Array[File] cell_type_beds", fixed = TRUE
  )
  testthat::expect_false(grepl("group_hdf5|group_tsv", text))
})

testthat::test_that("BuildManifest rewrites direct BED inventory paths", {
  text <- qc_wdl_text()

  purrr::walk(c(
    "localized_bed_files", "stage_localized_file()", "basename",
    "ln -s \"$source_absolute\" \"$destination\"",
    "~{sep('\\n', cell_type_beds)}",
    "rewritten_inventory=outputs/assembled_outputs.localized.tsv",
    "--outputs \"$rewritten_inventory\"",
    "if [[ ! -e \"$rewritten_path\" ]]",
    "File provenance = \"outputs/assembled_outputs.localized.tsv\""
  ), ~ testthat::expect_match(text, .x, fixed = TRUE))
})
