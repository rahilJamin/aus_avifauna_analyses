library(testthat)
Sys.setenv(TESTTHAT_EDITION = "3")

find_repo_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    if (file.exists(file.path(current, "R", "00_main_config.R"))) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not find repository root containing R/00_main_config.R.", call. = FALSE)
    }

    current <- parent
  }
}

repo_root <- find_repo_root()
test_dir(file.path(repo_root, "tests", "testthat"))
