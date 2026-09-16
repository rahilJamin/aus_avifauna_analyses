import::from(readr, write_csv)

#' Normalise a path for safety checks
#'
#' Resolves a file or directory path as far as possible, preserves not-yet-created
#' path segments, converts separators to forward slashes, and collapses simple
#' dot segments. This helper is used by output guards before generated files are
#' written.
#'
#' @param path A single file or directory path.
#'
#' @return A normalised character path.

normalise_path_for_check <- function(path) {
  if (length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("Expected one non-empty output path.", call. = FALSE)
  }

  # Generated files may not exist yet, so resolve the nearest existing parent
  # and then append the missing path segments.
  ancestor <- path.expand(path)
  suffix <- character()

  while (!file.exists(ancestor) && !dir.exists(ancestor)) {
    parent <- dirname(ancestor)
    if (identical(parent, ancestor)) break

    suffix <- c(basename(ancestor), suffix)
    ancestor <- parent
  }

  resolved <- normalizePath(ancestor, winslash = "/", mustWork = FALSE)

  if (length(suffix) > 0) {
    resolved <- do.call(file.path, c(list(resolved), as.list(suffix)))
  }

  resolved <- gsub("\\\\", "/", resolved)

  # Collapse dot segments while preserving Windows drives and UNC roots.
  prefix <- if (startsWith(resolved, "//")) {
    "//"
  } else if (startsWith(resolved, "/")) {
    "/"
  } else {
    ""
  }

  parts <- strsplit(sub("^/+", "", resolved), "/", fixed = TRUE)[[1]]
  root_depth <- if (prefix == "//") {
    2L
  } else if (grepl("^[A-Za-z]:$", parts[1])) {
    1L
  } else {
    0L
  }

  clean <- character()

  for (part in parts) {
    if (part %in% c("", ".")) next

    if (part == "..") {
      if (length(clean) > root_depth) {
        clean <- clean[-length(clean)]
      }
    } else {
      clean <- c(clean, part)
    }
  }

  paste0(prefix, paste(clean, collapse = "/"))
}

#' Test whether a path is inside a parent directory
#'
#' Compares normalised paths and handles case-insensitive comparison on Windows.
#' It is used to prevent generated outputs from being written into protected
#' input directories.
#'
#' @param path A single file or directory path to check.
#' @param parent A single parent directory path.
#'
#' @return A logical value.

path_is_inside <- function(path, parent) {
  path <- normalise_path_for_check(path)
  parent <- normalise_path_for_check(parent)

  if (.Platform$OS.type == "windows") {
    path <- tolower(path)
    parent <- tolower(parent)
  }

  path == parent || startsWith(path, paste0(sub("/+$", "", parent), "/"))
}

#' Refuse output paths inside protected input directories
#'
#' Stops if a generated output path points inside known raw-data/input folders.
#' When `cfg$raw_data_dir` is available, that configured directory is also
#' treated as protected.
#'
#' @param path A single output file or directory path.
#'
#' @return The input path, invisibly, if it passes the guard.

assert_not_raw_output_path <- function(path) {
  raw_dirs <- c("data_raw", "data", "spatial_grids")

  if (exists("cfg", inherits = TRUE) && !is.null(cfg$raw_data_dir)) {
    raw_dirs <- unique(c(raw_dirs, cfg$raw_data_dir))
  }

  writes_to_raw <- vapply(
    raw_dirs,
    function(raw_dir) path_is_inside(path, raw_dir),
    logical(1)
  )

  if (any(writes_to_raw)) {
    stop(
      "Refusing to write into a raw-data directory: ",
      normalise_path_for_check(path),
      call. = FALSE
    )
  }

  invisible(path)
}

#' Write a CSV after checking the output path
#'
#' Creates the parent directory if needed, writes a CSV with `readr::write_csv()`,
#' and refuses to write inside protected input directories.
#'
#' @param x A data frame or tibble to write.
#' @param path Output CSV path.
#'
#' @return The output path, invisibly.

write_csv_quiet <- function(x, path) {
  assert_not_raw_output_path(path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write_csv(x, path)
  invisible(path)
}

#' Save an R object after checking the output path
#'
#' Creates the parent directory if needed, saves an RDS object with `saveRDS()`,
#' and refuses to write inside protected input directories.
#'
#' @param x An R object to save.
#' @param path Output RDS path.
#'
#' @return The output path, invisibly.

save_rds_quiet <- function(x, path) {
  assert_not_raw_output_path(path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, path)
  invisible(path)
}

#' Write text lines after checking the output path
#'
#' Creates the parent directory if needed, writes lines with `writeLines()`, and
#' refuses to write inside protected input directories.
#'
#' @param text Character vector to write.
#' @param path Output text-file path.
#'
#' @return The output path, invisibly.
write_lines_quiet <- function(text, path) {
  assert_not_raw_output_path(path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(text, path)
  invisible(path)
}
