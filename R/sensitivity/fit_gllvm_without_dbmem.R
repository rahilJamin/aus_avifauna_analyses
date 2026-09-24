#' Prepare exactly the main GLLVM response and predictor data
#' @param config Main configuration.
#' @return Aligned data produced by the same main preparation functions.
prepare_sensitivity_gllvm <- function(config) {
  sites <- assign_ibra_subregions(readRDS(config$model_df_path),
                                  config$ibra_path, config$model_coord_crs)
  prepare_gllvm_model_data(sites, readRDS(config$species_matrix_path),
    readr::read_csv(config$trait_path, show_col_types = FALSE),
    config$gllvm_base_predictors, config$occupancy_threshold, config$noisy_miner_species)
}

#' Fit the matched environmental GLLVM with dbMEM predictors removed
#'
#' Called in a fresh R process. The response, base covariates, latent dimension,
#' random intercept, variational method and all three starts remain unchanged.
#' @param prepared Main GLLVM preparation bundle.
#' @param main_model_path Completed environmental model, read only for checks.
#' @param config,settings Main and sensitivity settings.
#' @return The dedicated sensitivity model path.
fit_gllvm_without_dbmem <- function(prepared, main_model_path, config, settings) {
  reference <- readRDS(main_model_path)
  assert_same_response(as.matrix(reference$y), as.matrix(prepared$Y))
  if (is.null(reference$X) || !all(config$gllvm_base_predictors %in% colnames(reference$X)) ||
      !isTRUE(all.equal(unname(as.matrix(reference$X[, config$gllvm_base_predictors, drop = FALSE])),
                        unname(as.matrix(prepared$X)), tolerance = 1e-12))) {
    stop("Prepared base predictors differ from the completed environmental GLLVM.", call. = FALSE)
  }
  # Release the multi-gigabyte main fit BEFORE constructing a new TMB objective.
  rm(reference)
  invisible(gc())
  TMB::openmp(n = config$model_threads)
  model <- gllvm::gllvm(
    y = as.matrix(prepared$Y), X = prepared$X,
    formula = define_gllvm_formulas()$environment,
    family = stats::binomial(link = "probit"), num.lv = config$gllvm_env_num_lv,
    row.eff = ~ (1 | ibra_sub), studyDesign = data.frame(ibra_sub = prepared$gllvm_df$ibra_sub),
    method = "VA", starting.val = "res", n.init = config$gllvm_n_init,
    seed = config$gllvm_seed, trace = TRUE
  )
  assert_same_response(as.matrix(prepared$Y), as.matrix(model$y))
  path <- sensitivity_path(settings$model_dir, "gllvm_environment_without_dbmem.rds")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(model, path)
  path
}

#' Run one memory-intensive sensitivity task in a fresh R process
#' @param task Name of the fitting or diagnostic function to call.
#' @param args Explicit arguments (only compact data or paths).
#' @param source_files Tracked R source paths loaded in the child process.
#' @return Task return value; child memory is released on completion.
sensitivity_worker <- function(task, args, source_files) {
  callr::r(function(task, args, files) {
    for (file in files) sys.source(file, envir = .GlobalEnv)
    do.call(get(task, envir = .GlobalEnv), args)
  }, args = list(task = task, args = args, files = source_files),
  libpath = .libPaths(), wd = getwd(), show = TRUE,
  system_profile = FALSE, user_profile = FALSE)
}
