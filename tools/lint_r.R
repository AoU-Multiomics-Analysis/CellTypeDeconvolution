#!/usr/bin/env Rscript

if (!requireNamespace("lintr", quietly = TRUE)) {
  stop("The lintr package is required. Install it with install.packages('lintr').", call. = FALSE)
}

paths <- c("scripts", "R/scripts", "R/utils")
paths <- paths[dir.exists(paths)]

r_files <- unlist(
  lapply(paths, list.files, pattern = "[.]R$", full.names = TRUE, recursive = TRUE),
  use.names = FALSE
)

if (length(r_files) == 0) {
  message("No R files found to lint.")
  quit(status = 0)
}

linters <- list(
  T_and_F_symbol_linter = lintr::T_and_F_symbol_linter(),
  equals_na_linter = lintr::equals_na_linter(),
  missing_argument_linter = lintr::missing_argument_linter(),
  semicolon_linter = lintr::semicolon_linter(),
  seq_linter = lintr::seq_linter(),
  unreachable_code_linter = lintr::unreachable_code_linter(),
  vector_logic_linter = lintr::vector_logic_linter()
)

lints <- unlist(
  lapply(r_files, lintr::lint, linters = linters, parse_settings = FALSE),
  recursive = FALSE
)

if (length(lints) > 0) {
  print(lints)
  quit(status = 1)
}

message("R lint ok")
