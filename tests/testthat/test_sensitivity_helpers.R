test_that("spatial diagnostic plots retain signed statistics and missing bands", {
  source_sensitivity_function("produce_sensitivity_outputs.R")
  dat <- data.frame(series = "total", lower_km = c(0, 2.5, 5),
    upper_km = c(2.5, 5, 7.5), midpoint_km = c(1.25, 3.75, 6.25),
    status = c("ok", "no_pairs", "ok"), moran_i = c(-0.1, 0.9, 0.2),
    expected_i = -0.01, p_bh = c(0.01, NA_real_, 0.2),
    model = "with_gp", focal_gradient = "agriculture")
  plots <- plot_sensitivity_spatial_diagnostics(data.frame(), dat)
  expect_equal(plots$lcbd$data$moran_i, c(-0.1, NA_real_, 0.2))
  expect_identical(plots$lcbd$data$detection,
    c("BH p < 0.05", NA_character_, "BH p >= 0.05"))
  lines <- ggplot2::ggplot_build(plots$lcbd)$data[[2]]
  expect_equal(lines$y, c(-0.1, NA_real_, 0.2))
  dat$model <- "environment_with_dbmem"
  plots <- plot_sensitivity_spatial_diagnostics(dat, data.frame())
  expect_equal(plots$gllvm$data$moran_i, c(-0.1, NA_real_, 0.2))
  expect_length(plot_sensitivity_spatial_diagnostics(data.frame(), data.frame()), 0L)
})

test_that("LCBD sensitivity plots retain interior curvature and saved uncertainty", {
  source_sensitivity_function("produce_sensitivity_outputs.R")
  curves <- data.frame(
    focal_gradient = "agriculture", focal_label = "Agricultural cover",
    scenario_id = "main_70", lcbd_component = "replacement",
    gradient_value = c(-1, 0, 1), mean_partial_effect = c(0.1, -0.2, 0.1),
    lower_95 = c(0, -0.3, 0), upper_95 = c(0.2, -0.1, 0.2)
  )
  manifest <- data.frame(scenario_id = c("main_70", "completeness_90"),
    scenario_label = c("Baseline", "90% completeness"),
    status = c("complete", "unsupported"))
  plots <- plot_sensitivity_lcbd_curves(curves, manifest)
  layers <- ggplot2::ggplot_build(plots$agriculture)$data
  # Equal endpoints must not erase the interior trough or its uncertainty.
  expect_equal(layers[[3]]$y, curves$mean_partial_effect)
  expect_equal(layers[[3]]$x, curves$gradient_value)
  expect_equal(layers[[2]]$ymin, curves$lower_95)
  expect_equal(layers[[2]]$ymax, curves$upper_95)
  expect_length(unique(layers[[3]]$PANEL), 1L)
  expect_length(plot_sensitivity_lcbd_curves(curves[FALSE, ], manifest), 0L)
  curves$scenario_id <- "completeness_90"
  expect_error(plot_sensitivity_lcbd_curves(curves, manifest), "completed")
})

test_that("sensitivity configuration declares the complete planned scenarios", {
  source_config()
  source_sensitivity_config()
  source_sensitivity_function("prepare_sampling_sensitivity.R")

  scenarios <- define_sensitivity_scenarios(cfg, sensitivity_cfg)

  expect_named(scenarios, c(
    "main_70",
    "completeness_75", "completeness_80", "completeness_85", "completeness_90",
    "effort_12", "daily_incidence", "structured_protocol", "pools_20", "pools_30"
  ))
  expect_equal(sensitivity_cfg$n_lcbd_allocations, 30L)
  expect_equal(sensitivity_cfg$allocation_seeds, 42:71)
  expect_identical(scenarios$main_70$type, "main")
  expect_true(all(vapply(scenarios, function(x) x$pool_k >= 20L, logical(1))))
  expect_equal(scenarios$effort_12$minimum_units, 12L)
})

test_that("sensitivity output paths cannot leave their dedicated directory", {
  source_sensitivity_function("validate_sensitivity_inputs.R")

  valid <- sensitivity_path(file.path("outputs", "sensitivity", "tables"), "summary.csv")
  expect_match(gsub("\\\\", "/", valid), "^outputs/sensitivity/tables/summary[.]csv$")
  expect_error(
    sensitivity_path(file.path("outputs", "sensitivity", "..", "tables"), "summary.csv"),
    "stay inside"
  )
  expect_error(
    sensitivity_path(file.path("outputs", "sensitivity", "tables"), "../summary.csv"),
    "stay inside"
  )
})

test_that("higher-completeness sensitivity preserves matrix row order", {
  source_config()
  source_function("build_matrix_bundle.R")
  source_sensitivity_function("validate_sensitivity_inputs.R")

  chao <- data.frame(
    hex_id = c("h3", "h1", "h2", "failed"),
    chao_status = c("ok", "ok", "ok", "failed"),
    sampling_units = c(8, 7, 6, 20),
    S_obs = c(6, 6, 6, 20),
    chao2_completeness = c(0.9, 0.8, 0.79, 1)
  )
  response <- matrix(
    c(1, 0, 0, 1, 1, 1),
    nrow = 3,
    dimnames = list(c("h1", "h2", "h3"), c("sp_a", "sp_b"))
  )

  selected <- sensitivity_site_ids(chao, response, threshold = 0.8, cfg)

  # Chao rows are deliberately not in matrix order. Matrix order governs every
  # downstream model alignment, so the sensitivity must retain h1 then h3.
  expect_identical(selected, c("h1", "h3"))
  expect_error(
    sensitivity_site_ids(chao, response, threshold = cfg$min_completeness, cfg),
    "must exceed"
  )
})

test_that("paired GLLVM comparisons require identical values as well as names", {
  source_sensitivity_function("validate_sensitivity_inputs.R")

  reference <- matrix(c(1, 0, 0, 1), 2,
    dimnames = list(c("h1", "h2"), c("sp_a", "sp_b")))
  expect_no_error(assert_same_response(reference, reference))

  changed <- reference
  changed[1, 1] <- 0
  expect_error(assert_same_response(reference, changed), "identical response values")
  reordered <- reference[c("h2", "h1"), , drop = FALSE]
  expect_error(assert_same_response(reference, reordered), "identical response values")
})

test_that("Noisy Miner count checks allow exact zero and one detection rates", {
  source_sensitivity_function("validate_sensitivity_inputs.R")

  valid <- data.frame(
    n_checklists = c(2, 3),
    n_miner_checklists = c(0, 3),
    miner_detection_rate = c(0, 1)
  )
  expect_no_error(validate_miner_counts(valid))

  invalid_count <- valid
  invalid_count$n_miner_checklists[1] <- 3
  expect_error(validate_miner_counts(invalid_count), "Noisy Miner counts")
  invalid_rate <- valid
  invalid_rate$miner_detection_rate[2] <- 0.5
  expect_error(validate_miner_counts(invalid_rate), "Noisy Miner counts")
})

test_that("daily incidence units retain calendar dates and reject missing dates", {
  source_sensitivity_function("prepare_sampling_sensitivity.R")

  records <- data.frame(date = as.Date(c("2020-01-01", "2020-01-01", "2020-01-02")))
  labelled <- add_daily_incidence_unit(records, "day")
  expect_identical(labelled$incidence_unit, c("2020-01-01", "2020-01-01", "2020-01-02"))
  expect_error(add_daily_incidence_unit(records, "month_year"), "daily units")
  expect_error(add_daily_incidence_unit(transform(records, date = as.Date(NA)), "day"),
               "complete dates")
})

test_that("GLLVM residual alignment preserves orientation but never relabels data", {
  source_sensitivity_function("validate_sensitivity_inputs.R")

  response <- matrix(c(1, 0, 0, 1), 2,
    dimnames = list(c("h1", "h2"), c("sp_a", "sp_b")))
  transposed <- matrix(c(0.2, -0.3, 0.5, 0.1), 2,
    dimnames = list(c("sp_a", "sp_b"), c("h1", "h2")))
  aligned <- align_gllvm_residuals(transposed, response)
  expect_identical(dimnames(aligned), dimnames(response))
  expect_equal(aligned["h1", "sp_a"], 0.2)

  mismatched <- transposed
  rownames(mismatched)[1] <- "wrong_species"
  expect_error(align_gllvm_residuals(mismatched, response), "names do not match")
})

test_that("spatial residual diagnostics reject malformed inputs before testing", {
  source_sensitivity_function("calculate_residual_moran.R")

  residuals <- matrix(c(0.1, -0.1), ncol = 1,
    dimnames = list(c("h1", "h2"), "response"))
  duplicated_coordinates <- matrix(c(0, 0, 0, 0), ncol = 2,
    dimnames = list(c("h1", "h2"), c("x", "y")))
  expect_error(
    calculate_residual_moran(residuals, duplicated_coordinates, c(0, 2500), "B", "two.sided"),
    "unique projected coordinates"
  )
  valid_coordinates <- matrix(c(0, 0, 2500, 0), ncol = 2,
    dimnames = list(c("h1", "h2"), c("x", "y")), byrow = TRUE)
  expect_error(
    calculate_residual_moran(residuals, valid_coordinates, c(2500, 0), "B", "two.sided"),
    "strictly increasing"
  )
})

test_that("supplementary tables can be exported as SVG", {
  source_sensitivity_function("produce_sensitivity_outputs.R")

  path <- tempfile(fileext = ".svg")
  on.exit(unlink(path), add = TRUE)
  save_sensitivity_table_svg(
    data.frame(Scenario = c("Baseline", "Alternative"), Sites = c(10, 8)),
    path,
    title = "Table S0. Test export",
    note = "Synthetic formatting check."
  )
  expect_true(file.exists(path))
  expect_gt(file.info(path)$size, 0)
})
