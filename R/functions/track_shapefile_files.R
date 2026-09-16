#' Track every file belonging to an ESRI shapefile
#'
#' A shapefile is a collection of files with the same basename. This helper
#' returns all existing components so `targets` notices changes to attributes,
#' geometry indexes, projections, and text encodings as well as the `.shp` file.
#'
#' @param shp_path Path to the main `.shp` file.
#'
#' @return Character vector containing all files with the same basename. The
#'   required `.shp`, `.shx`, and `.dbf` components must be present.
track_shapefile_files <- function(shp_path) {
  if (!is.character(shp_path) || length(shp_path) != 1L ||
      is.na(shp_path) || !nzchar(shp_path) ||
      tolower(tools::file_ext(shp_path)) != "shp") {
    stop("`shp_path` must be one path ending in .shp.", call. = FALSE)
  }

  directory <- dirname(shp_path)
  basename_without_extension <- tools::file_path_sans_ext(basename(shp_path))

  if (!dir.exists(directory)) {
    stop("Missing shapefile directory: ", directory, call. = FALSE)
  }

  candidates <- list.files(directory, full.names = TRUE, all.files = FALSE)
  candidate_basenames <- tools::file_path_sans_ext(basename(candidates))
  files <- candidates[
    tolower(candidate_basenames) == tolower(basename_without_extension)
  ]

  extensions <- tolower(tools::file_ext(files))
  required_extensions <- c("shp", "shx", "dbf")
  missing_extensions <- setdiff(required_extensions, extensions)

  if (length(missing_extensions) > 0L) {
    stop(
      "The shapefile is missing required component(s): ",
      paste(paste0(".", missing_extensions), collapse = ", "),
      call. = FALSE
    )
  }

  files[order(match(extensions, c("shp", "shx", "dbf", "prj", "cpg", "qpj")), files)]
}

#' Select the main geometry file from tracked shapefile components
#'
#' @param shapefile_files Character vector returned by
#'   `track_shapefile_files()`.
#'
#' @return The single path ending in `.shp`.
main_shapefile_path <- function(shapefile_files) {
  shp <- shapefile_files[
    tolower(tools::file_ext(shapefile_files)) == "shp"
  ]

  if (length(shp) != 1L) {
    stop("Expected exactly one .shp file in the tracked components.", call. = FALSE)
  }

  shp
}

#' Select an expected file from a tracked file target
#'
#' File targets may return their paths without preserving vector names. This
#' helper matches the configured path itself, so downstream targets do not rely
#' on names or positions while still declaring a dependency on the whole file
#' target.
#'
#' @param tracked_files Character vector of paths returned by a file target.
#' @param expected_path Configured path of the required file.
#'
#' @return The matching path from `tracked_files`.
tracked_file_path <- function(tracked_files, expected_path) {
  if (!is.character(tracked_files) || length(tracked_files) < 1L ||
      anyNA(tracked_files) || any(!nzchar(tracked_files))) {
    stop("`tracked_files` must contain one or more file paths.", call. = FALSE)
  }

  if (!is.character(expected_path) || length(expected_path) != 1L ||
      is.na(expected_path) || !nzchar(expected_path)) {
    stop("`expected_path` must be one non-empty file path.", call. = FALSE)
  }

  normalise <- function(path) {
    normalised <- normalizePath(
      path,
      winslash = "/",
      mustWork = FALSE
    )

    if (.Platform$OS.type == "windows") tolower(normalised) else normalised
  }

  matches <- which(normalise(tracked_files) == normalise(expected_path))

  if (length(matches) != 1L) {
    stop(
      "Expected exactly one tracked file matching: ",
      expected_path,
      call. = FALSE
    )
  }

  unname(tracked_files[matches])
}
