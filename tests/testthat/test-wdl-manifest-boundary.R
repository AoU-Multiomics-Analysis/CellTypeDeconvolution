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
    text, "File qc_summary = \"outputs/qc_summary.tsv\"", fixed = TRUE
  )
  testthat::expect_match(
    text, "File qc_plots = export_qc_plots", fixed = TRUE
  )
  testthat::expect_match(
    text, "checksum_inventory.localized.tsv", fixed = TRUE
  )
  testthat::expect_match(
    text, "Array[File] cell_type_beds", fixed = TRUE
  )
  testthat::expect_false(grepl("group_hdf5|group_tsv", text))
})

testthat::test_that("BuildManifest keeps localized paths internal", {
  text <- qc_wdl_text()

  purrr::walk(c(
    "localized_bed_files", "stage_localized_file()", "basename",
    "ln -s \"$source_absolute\" \"$destination\"",
    "~{sep='\\n' cell_type_beds}",
    "checksum_inventory=checksum_inventory.localized.tsv",
    "public_inventory=outputs/output_inventory.tsv",
    "--outputs \"$checksum_inventory\"",
    "if [[ ! -e \"$rewritten_path\" ]]",
    "File provenance = \"outputs/output_inventory.tsv\""
  ), ~ testthat::expect_match(text, .x, fixed = TRUE))
  testthat::expect_false(grepl(
    'File provenance = "checksum_inventory.localized.tsv"',
    text,
    fixed = TRUE
  ))
})

testthat::test_that("top workflow publishes the final cross-stage QC summary", {
  text <- paste(
    readLines(testthat::test_path(
      "../..", "workflows", "cell_type_deconvolution.wdl"
    )),
    collapse = "\n"
  )

  testthat::expect_match(
    text,
    "dtangle_metadata = RunDtangle.metadata",
    fixed = TRUE
  )
  testthat::expect_match(
    text,
    "File qc_summary = BuildManifest.qc_summary",
    fixed = TRUE
  )
})
