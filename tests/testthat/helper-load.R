pipeline_root <- normalizePath(file.path(testthat::test_path(), "../.."))
r_files <- file.path(pipeline_root, "R", c("constants.R", "io.R"))
invisible(lapply(r_files[file.exists(r_files)], source, local = .GlobalEnv))
