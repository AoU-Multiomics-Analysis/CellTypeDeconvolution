file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]))
pipeline_root <- normalizePath(Sys.getenv(
  "CELLTYPE_ROOT",
  unset = file.path(dirname(script_path), "..")
))
r_files <- list.files(file.path(pipeline_root, "R"), pattern = "[.]R$", full.names = TRUE)
invisible(lapply(sort(r_files), source))
