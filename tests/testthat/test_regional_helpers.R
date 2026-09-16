test_that("regional gradient definitions use all four configured comparisons", {
  source_function("define_regional_gradients.R")

  site_data <- data.frame(
    agricultural_cover = c(0, 10, 30, 90, 100),
    mean_woodland_patch_area = c(1, 2, 3, 4, 5),
    woodland_cohesion = c(10, 20, 30, 40, 50),
    noisy_miner_detection = c(0, 0.1, 0.2, 0.8, 1)
  )

  gradients <- define_regional_gradients(
    site_data,
    agriculture_low_maximum = 15,
    agriculture_high_minimum = 85,
    lower_tail_probability = 0.1,
    upper_tail_probability = 0.9
  )

  expect_identical(
    gradients$gradient,
    c("agriculture", "patch_area", "cohesion", "noisy_miner")
  )
  expect_identical(
    gradients$predictor[gradients$gradient == "noisy_miner"],
    "noisy_miner_detection"
  )
  expect_true(all(gradients$low_threshold < gradients$high_threshold))
  expect_match(gradients$threshold_rule[2], "10th and 90th")
  expect_match(gradients$low_endpoint_definition[1], "< 15%")
  expect_match(gradients$high_endpoint_definition[1], "> 85%")

  expect_error(
    define_regional_gradients(site_data, 15, 85, 0.9, 0.1),
    "ordered"
  )
  constant <- transform(site_data, woodland_cohesion = 1)
  expect_error(
    define_regional_gradients(constant, 15, 85, 0.1, 0.9),
    "distinct low and high"
  )
})

test_that("regional Jaccard components are symmetric and additive", {
  source_function("calculate_regional_jaccard_components.R")

  species_matrix <- matrix(
    c(
      1, 0, 0,
      1, 1, 0,
      0, 1, 1,
      0, 0, 1
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("h", 1:4), paste0("sp", 1:3))
  )

  components <- calculate_regional_jaccard_components(species_matrix)
  matrices <- components$component_matrices

  expect_named(
    matrices,
    c("Total Jaccard", "Replacement", "Richness difference")
  )
  expect_lt(components$component_sum_error, 1e-10)
  expect_equal(matrices[["Total Jaccard"]], t(matrices[["Total Jaccard"]]))
  expect_equal(
    matrices[["Total Jaccard"]],
    matrices[["Replacement"]] + matrices[["Richness difference"]]
  )
  expect_equal(unname(diag(matrices[["Total Jaccard"]])), rep(0, 4))

  empty_site <- species_matrix
  empty_site[1, ] <- 0
  expect_error(
    calculate_regional_jaccard_components(empty_site),
    "at least one species"
  )
})

test_that("regional endpoint analysis respects strict agriculture boundaries", {
  source_function("assert_unique_ids.R")
  source_function("analyse_regional_gradient.R")

  subregions <- rep(c("A", "B", "C"), each = 6)
  predictor <- rep(c(0, 1, 2, 8, 9, 10), times = 3)
  site_ids <- paste0("site_", seq_along(predictor))
  site_data <- data.frame(
    hex_id = site_ids,
    ibra_subregion = subregions,
    site_row = seq_along(site_ids),
    agricultural_cover = predictor
  )
  species_matrix <- cbind(
    sp_a = rep(1L, length(site_ids)),
    sp_b = as.integer(seq_along(site_ids) %% 2L == 0L),
    sp_c = as.integer(seq_along(site_ids) %% 3L == 0L)
  )
  rownames(species_matrix) <- site_ids

  distance_matrix <- outer(
    seq_along(site_ids),
    seq_along(site_ids),
    function(i, j) abs(i - j) / length(site_ids)
  )
  dimnames(distance_matrix) <- list(site_ids, site_ids)
  component_matrices <- list("Total Jaccard" = distance_matrix)

  gradient_definition <- data.frame(
    gradient = "agriculture",
    gradient_label = "Agricultural cover",
    predictor = "agricultural_cover",
    low_threshold = 2,
    high_threshold = 8,
    threshold_rule = "Fixed percentage of land cover",
    low_endpoint_definition = "Agricultural cover < 2%",
    high_endpoint_definition = "Agricultural cover > 8%"
  )

  result_a <- analyse_regional_gradient(
    gradient_definition = gradient_definition,
    site_data = site_data,
    species_matrix = species_matrix,
    component_matrices = component_matrices,
    gradient_index = 1L,
    minimum_sites_per_endpoint = 2L,
    minimum_retained_subregions = 3L,
    bootstrap_draws = 9L,
    bootstrap_seed = 123L
  )
  result_b <- analyse_regional_gradient(
    gradient_definition = gradient_definition,
    site_data = site_data,
    species_matrix = species_matrix,
    component_matrices = component_matrices,
    gradient_index = 1L,
    minimum_sites_per_endpoint = 2L,
    minimum_retained_subregions = 3L,
    bootstrap_draws = 9L,
    bootstrap_seed = 123L
  )

  expect_true(all(result_a$regional_support$retained))
  expect_true(all(result_a$regional_support$n_sites_low == 2L))
  expect_true(all(result_a$regional_support$n_sites_high == 2L))
  expect_equal(nrow(result_a$subregion_pair_means), 6L)
  expect_equal(nrow(result_a$bootstrap_results), 18L)
  expect_identical(result_a$bootstrap_results, result_b$bootstrap_results)
  expect_false(any(
    result_a$subregion_pair_means$n_sites_subregion_1 != 2L |
      result_a$subregion_pair_means$n_sites_subregion_2 != 2L
  ))

  misaligned <- species_matrix
  rownames(misaligned)[1] <- "wrong"
  expect_error(
    analyse_regional_gradient(
      gradient_definition,
      site_data,
      misaligned,
      component_matrices,
      1L,
      2L,
      3L,
      9L,
      123L
    ),
    "misaligned"
  )
})
