# Structural checks cannot establish which inputs or settings produced an object.
# Use trusted workflow files and confirm their provenance before a manual run.
#' Read an RDS object with simple structural checks
#'
#' Reads a required workflow checkpoint and gives clearer errors when the file is
#' missing, unreadable, has the wrong class, or lacks required top-level names.
#' This function checks structure only; it does not prove the object was created
#' with the current code or settings.
#'
#' @param path Path to an RDS file.
#' @param label Human-readable object label used in error messages.
#' @param required_class Optional class that the object must inherit from.
#' @param required_names Optional character vector of required top-level names.
#'
#' @return The object read from `path`.

read_rds_checked <- function(path, label, required_class = NULL,
                             required_names = character()) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("No input path configured for ", label, ".", call. = FALSE)
  }
  if (!file.exists(path) || dir.exists(path)) {
    stop("Missing ", label, ": ", path,
         "\nPrepare this object first or correct its path in R/00_main_config.R.",
         call. = FALSE)
  }

  object <- tryCatch(
    readRDS(path),
    error = function(e) {
      stop("Cannot read ", label, " from ", path, ": ", conditionMessage(e),
           call. = FALSE)
    }
  )
  if (!is.null(required_class) && !inherits(object, required_class)) {
    stop("Invalid ", label, " at ", path, ": expected class ",
         paste(required_class, collapse = " or "), ".", call. = FALSE)
  }
  missing_names <- setdiff(required_names, names(object))
  if (length(missing_names) > 0) {
    stop("Invalid ", label, " at ", path, ": missing ",
         paste(missing_names, collapse = ", "), ".", call. = FALSE)
  }
  object
}
