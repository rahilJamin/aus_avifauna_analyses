#' Define the two main GLLVM formulae
#'
#' Keeps the environmental and fourth-corner model formulae in one inspectable,
#' testable place. Selected dbMEM terms are appended to both formulae.
#'
#' @param selected_mem_names Character vector of selected dbMEM column names.
#'
#' @return A named list containing `environment` and `fourth_corner` formulae.
define_gllvm_formulas <- function(selected_mem_names = character()) {
  if (!is.character(selected_mem_names) || anyNA(selected_mem_names) ||
      any(!nzchar(selected_mem_names)) || anyDuplicated(selected_mem_names) ||
      any(make.names(selected_mem_names) != selected_mem_names)) {
    stop(
      "Selected dbMEM names must be unique, non-empty syntactic names.",
      call. = FALSE
    )
  }

  mem_term <- if (length(selected_mem_names) > 0L) {
    paste("+", paste(selected_mem_names, collapse = " + "))
  } else {
    ""
  }

  environment <- stats::as.formula(paste(
    "~ pland_agriculture_2500m +
       area_mn_woodland_2500m +
       cohesion_woodland_2500m +
       log_effort",
    mem_term
  ))

  fourth_corner <- stats::as.formula(paste(
    "~ (pland_agriculture_2500m +
         area_mn_woodland_2500m +
         cohesion_woodland_2500m) :
       (log_mass + HWI) +
       log_effort",
    mem_term
  ))

  list(environment = environment, fourth_corner = fourth_corner)
}
