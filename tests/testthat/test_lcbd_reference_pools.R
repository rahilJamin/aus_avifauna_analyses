test_that("reference-pool allocation is reproducible and preserves site IDs", {
  source_function("assert_unique_ids.R")
  source_function("assign_lcbd_reference_pools.R")

  site_data <- data.frame(
    hex_id = paste0("hex_", seq_len(8)),
    x_km = c(0, 0.1, 0.2, 0.3, 10, 10.1, 10.2, 10.3),
    y_km = c(0, 0.1, 0.2, 0.3, 10, 10.1, 10.2, 10.3)
  )

  allocation_a <- assign_lcbd_reference_pools(
    site_data = site_data,
    allocation = 1L,
    allocation_seed = 42L,
    reference_pool_k = 2L,
    reference_pool_nstart = 10L,
    minimum_pool_sites = 2L
  )
  allocation_b <- assign_lcbd_reference_pools(
    site_data = site_data,
    allocation = 1L,
    allocation_seed = 42L,
    reference_pool_k = 2L,
    reference_pool_nstart = 10L,
    minimum_pool_sites = 2L
  )

  expect_identical(
    allocation_a$allocation_membership$reference_pool,
    allocation_b$allocation_membership$reference_pool
  )
  expect_setequal(allocation_a$allocation_membership$hex_id, site_data$hex_id)
  expect_equal(nrow(allocation_a$allocation_membership), nrow(site_data))
  expect_true(all(allocation_a$pool_sizes$reference_pool_size >= 2L))
  expect_true(all(allocation_a$allocation_membership$complete_pool))
})

test_that("reference-pool allocation rejects impossible clustering settings", {
  source_function("assert_unique_ids.R")
  source_function("assign_lcbd_reference_pools.R")

  site_data <- data.frame(
    hex_id = c("a", "b", "c"),
    x_km = c(0, 0, 1),
    y_km = c(0, 0, 1)
  )

  expect_error(
    assign_lcbd_reference_pools(site_data, 1L, 42L, 3L, 5L, 1L),
    "distinct coordinates"
  )
  expect_error(
    assign_lcbd_reference_pools(site_data, 1L, 42L, 2L, 5L, 3L),
    "No reference pool"
  )
})

test_that("gradient support retains all sites from informative pools", {
  source_function("assert_unique_ids.R")
  source_function("calculate_lcbd_gradient_support.R")

  focal_predictors <- c(
    agriculture = "agriculture",
    patch_area = "patch_area",
    cohesion = "cohesion",
    noisy_miner = "miner_detection_rate"
  )
  gradient_values <- c(1, 2, 3, 7, 8, 9, 4, 4.5, 5, 5.5)

  allocation_lcbd_data <- data.frame(
    hex_id = paste0("h", seq_along(gradient_values)),
    reference_pool = factor(c(rep("1", 6), rep("2", 4))),
    agriculture = gradient_values,
    patch_area = gradient_values,
    cohesion = gradient_values,
    miner_detection_rate = gradient_values
  )
  thresholds <- data.frame(
    focal_predictor = unname(focal_predictors),
    global_low_threshold = 3,
    global_high_threshold = 7
  )

  support <- calculate_lcbd_gradient_support(
    allocation_lcbd_data = allocation_lcbd_data,
    focal_predictors = focal_predictors,
    global_gradient_thresholds = thresholds,
    allocation = 1L,
    allocation_seed = 42L,
    minimum_pool_sites = 4L,
    minimum_tail_sites = 2L
  )

  expect_named(
    support$gradient_pools,
    c("agriculture", "patch_area", "cohesion", "noisy_miner")
  )
  expect_identical(support$gradient_pools$noisy_miner, "1")
  expect_equal(nrow(support$gradient_data$noisy_miner), 6L)
  expect_true(all(
    support$gradient_data$noisy_miner$reference_pool == "1"
  ))
  expect_equal(
    sum(support$gradient_support$informative_pool),
    length(focal_predictors)
  )
})

test_that("LCBD pool validation enforces sums and retains negative diagnostics", {
  source_function("calculate_reference_pool_lcbd.R")

  checks <- data.frame(
    raw_lcbd_sum = c(1, 1 + 1e-8),
    n_values_below_negative_tolerance = c(0L, 0L)
  )
  expect_no_error(validate_lcbd_pool_checks(checks, 1L, 1e-10, 1e-6))

  bad_sum <- checks
  bad_sum$raw_lcbd_sum[2] <- 1.01
  expect_error(
    validate_lcbd_pool_checks(bad_sum, 1L, 1e-10, 1e-6),
    "did not sum to one"
  )

  bad_negative <- checks
  bad_negative$n_values_below_negative_tolerance[1] <- 1L
  expect_no_error(
    validate_lcbd_pool_checks(bad_negative, 1L, 1e-10, 1e-6)
  )
})

test_that("within-pool LCBD returns aligned values for every component", {
  source_function("assert_unique_ids.R")
  source_function("calculate_lcbd_beta_matrices.R")
  source_function("calculate_reference_pool_lcbd.R")

  species_matrix <- matrix(
    c(
      1, 0, 0,
      0, 1, 0,
      1, 0, 1,
      0, 1, 1
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("site_", 1:4), paste0("sp_", 1:3))
  )
  distance_matrices <- calculate_lcbd_beta_matrices(species_matrix)$dist_mats
  complete_pool_data <- data.frame(
    hex_id = rownames(species_matrix),
    reference_pool = factor("pool_1"),
    reference_pool_size = 4L
  )

  result <- calculate_reference_pool_lcbd(
    complete_pool_data = complete_pool_data,
    distance_matrices = distance_matrices,
    allocation = 1L,
    allocation_seed = 42L,
    negative_tolerance = 1e-10,
    lcbd_sum_tolerance = 1e-6
  )

  expect_identical(result$allocation_lcbd_data$hex_id, rownames(species_matrix))
  expect_equal(nrow(result$allocation_lcbd_data), 4L)
  expect_equal(result$pool_checks$raw_lcbd_sum, c(1, 1, 1), tolerance = 1e-6)
  expect_true(all(
    result$allocation_lcbd_data$relative_lcbd_total ==
      4 * result$allocation_lcbd_data$raw_lcbd_total
  ))
})

test_that("ensemble summaries require complete allocation and curve coverage", {
  source_function("summarise_lcbd_ensemble.R")

  components <- c("total", "replacement", "richness_difference")
  result_rows <- expand.grid(
    allocation = 1:2,
    lcbd_component = components,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  result_rows$focal_gradient <- "agriculture"
  result_rows$focal_label <- "Agricultural cover"
  result_rows$focal_variable <- "agriculture"
  result_rows$n_sites <- 100
  result_rows$n_informative_pools <- 10
  result_rows$endpoint_low <- 0
  result_rows$endpoint_high <- 1
  result_rows$endpoint_difference <- c(0.1, 0.2, -0.1, -0.2, 0.05, 0.15)
  result_rows$endpoint_se <- 0.05
  result_rows$edf <- 2
  result_rows$deviance_explained <- 0.2
  result_rows$model_converged <- TRUE

  curve_rows <- expand.grid(
    allocation = 1:2,
    lcbd_component = components,
    gradient_value = c(0, 0.5, 1),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  curve_rows$focal_gradient <- "agriculture"
  curve_rows$focal_label <- "Agricultural cover"
  curve_rows$focal_variable <- "agriculture"
  curve_rows$partial_effect <- curve_rows$gradient_value
  curve_rows$partial_effect_se <- 0.05

  allocation_summary <- data.frame(
    allocation = 1:2,
    partition_id = c(1L, 2L)
  )

  ensemble <- summarise_lcbd_ensemble(
    per_allocation_results = result_rows,
    per_allocation_curves = curve_rows,
    allocation_summary = allocation_summary,
    expected_model_count = 6L
  )
  expect_equal(nrow(ensemble$main_results), 3L)
  expect_true(all(ensemble$main_results$n_allocations == 2L))
  expect_true(all(ensemble$main_curves$n_allocations == 2L))

  incomplete_curves <- curve_rows[-1, ]
  expect_error(
    summarise_lcbd_ensemble(
      result_rows,
      incomplete_curves,
      allocation_summary,
      6L
    ),
    "prediction grids differ"
  )

  duplicated_results <- rbind(result_rows, result_rows[1, ])
  expect_error(
    summarise_lcbd_ensemble(
      duplicated_results,
      curve_rows,
      allocation_summary,
      7L
    ),
    "duplicated"
  )
})
