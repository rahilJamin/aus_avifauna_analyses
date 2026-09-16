make_test_distance_bundle <- function(ids = c("h1", "h2")) {
  total <- matrix(
    c(0, 1, 1, 0),
    nrow = 2,
    dimnames = list(ids, ids)
  )
  list(
    total = total,
    replacement = total,
    richness_difference = total * 0
  )
}

test_that("LCBD thresholds use every named gradient and return plain numbers", {
  source_function("calculate_lcbd_gradient_thresholds.R")

  site_data <- data.frame(
    agriculture = 1:8,
    woodland = seq(10, 80, by = 10),
    unused = 99
  )
  focal_predictors <- c(
    agriculture_gradient = "agriculture",
    woodland_gradient = "woodland"
  )

  thresholds <- calculate_lcbd_gradient_thresholds(
    site_data,
    focal_predictors = focal_predictors,
    lower_tail_probability = 0.25,
    upper_tail_probability = 0.75
  )

  expect_setequal(thresholds$focal_predictor, c("agriculture", "woodland"))
  expect_equal(
    thresholds$global_low_threshold[
      thresholds$focal_predictor == "agriculture"
    ],
    as.numeric(stats::quantile(1:8, 0.25)),
    ignore_attr = TRUE
  )
  expect_null(names(thresholds$global_low_threshold))
})

test_that("LCBD thresholds reject ambiguous or unusable gradients", {
  source_function("calculate_lcbd_gradient_thresholds.R")

  site_data <- data.frame(varying = 1:4, constant = 1)

  expect_error(
    calculate_lcbd_gradient_thresholds(
      site_data,
      focal_predictors = c(varying = "varying", constant = "constant"),
      lower_tail_probability = 0.25,
      upper_tail_probability = 0.75
    ),
    "no distinct low and high"
  )
  expect_error(
    calculate_lcbd_gradient_thresholds(
      site_data,
      focal_predictors = c(a = "varying", b = "varying"),
      lower_tail_probability = 0.25,
      upper_tail_probability = 0.75
    ),
    "uniquely named"
  )
  expect_error(
    calculate_lcbd_gradient_thresholds(
      site_data,
      focal_predictors = c(varying = "varying"),
      lower_tail_probability = 0.8,
      upper_tail_probability = 0.2
    ),
    "ordered"
  )
})

test_that("LCBD prediction grids use ordered shared limits", {
  source_function("make_lcbd_prediction_grids.R")

  limits <- data.frame(
    focal_gradient = c("agriculture", "woodland"),
    curve_lower = c(0, 10),
    curve_upper = c(1, 20)
  )

  grids <- make_lcbd_prediction_grids(limits, curve_grid_points = 5L)

  expect_named(grids, c("agriculture", "woodland"))
  expect_equal(grids$agriculture, seq(0, 1, length.out = 5L))
  expect_equal(grids$woodland, seq(10, 20, length.out = 5L))
  expect_error(
    make_lcbd_prediction_grids(transform(limits, curve_upper = curve_lower), 5L),
    "ordered"
  )
  expect_error(make_lcbd_prediction_grids(limits, 1L), "at least two")
})

test_that("LCBD formulae retain all scientific and adjustment terms", {
  source_function("define_lcbd_formulas.R")

  formulas <- define_lcbd_formulas()
  expect_named(formulas, c("total", "replacement", "richness_difference"))

  expected_responses <- c(
    total = "relative_lcbd_total",
    replacement = "relative_lcbd_replacement",
    richness_difference = "relative_lcbd_richness_difference"
  )
  expected_predictors <- c(
    "pland_agriculture_2500m",
    "area_mn_woodland_2500m",
    "cohesion_woodland_2500m",
    "miner_detection_rate",
    "log_effort",
    "x_km",
    "y_km",
    "reference_pool"
  )

  for (component in names(formulas)) {
    variables <- all.vars(formulas[[component]])
    expect_identical(variables[1], expected_responses[[component]])
    expect_setequal(variables[-1], expected_predictors)

    formula_text <- paste(deparse(formulas[[component]]), collapse = " ")
    expect_match(formula_text, 'bs = "gp"')
    expect_match(formula_text, 'bs = "re"')
    expect_match(formula_text, "k = 50")
  }
})

test_that("saved and rebuilt LCBD matrices are compared by site ID", {
  source_function("validate_lcbd_beta_matrices.R")

  rebuilt <- make_test_distance_bundle()
  reversed_ids <- c("h2", "h1")
  saved <- lapply(
    rebuilt,
    function(x) x[reversed_ids, reversed_ids, drop = FALSE]
  )

  validation <- compare_lcbd_beta_matrices(
    saved_beta_object = list(dist_mats = saved),
    rebuilt_distance_matrices = rebuilt,
    ordered_site_ids = c("h1", "h2")
  )
  expect_equal(validation$maximum_absolute_difference, c(0, 0, 0))

  changed <- saved
  changed$total["h1", "h2"] <- 0.9
  changed$total["h2", "h1"] <- 0.9
  expect_error(
    compare_lcbd_beta_matrices(
      list(dist_mats = changed),
      rebuilt,
      c("h1", "h2")
    ),
    "not identical"
  )
  expect_error(
    compare_lcbd_beta_matrices(
      list(dist_mats = saved["total"]),
      rebuilt,
      c("h1", "h2")
    ),
    "missing"
  )
  expect_error(
    compare_lcbd_beta_matrices(
      list(dist_mats = saved),
      rebuilt,
      c("h1", "h1")
    ),
    "unique"
  )
})
