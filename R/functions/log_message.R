#' Log a timestamped message
#'
#' Prints a timestamped message to the console and, optionally, appends it to a
#' log file after checking that the log path is not inside a protected input
#' directory.
#'
#' @param ... Message components passed to `paste()`.
#' @param log_file Optional path to a log file.
#'
#' @return The timestamped message, invisibly.

log_message <- function(..., log_file = NULL) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ",
                paste(..., collapse = ""))
  message(msg)
  if (!is.null(log_file)) {
    assert_not_raw_output_path(log_file)
    dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
    cat(msg, "\n", file = log_file, append = TRUE)
  }
  invisible(msg)
}
