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

source_function <- function(file, envir = parent.frame()) {
  source(file.path(find_repo_root(), "R", "functions", file), local = envir)
}

source_config <- function(envir = parent.frame()) {
  source(file.path(find_repo_root(), "R", "00_main_config.R"), local = envir)
}
