make_gllvm_test_inputs <- function() {
  model_df <- data.frame(
    hex_id = paste0("hex_", 1:5),
    ibra_sub = c("A", "A", "B", "B", "B"),
    x = 1:5,
    y = 11:15,
    agriculture = c(0, 1, 2, 3, NA),
    patch_area = 5:1,
    cohesion = c(10, 20, 30, 40, 50),
    effort = log(1:5)
  )
  species_matrix <- matrix(
    c(
      1, 1, 1, 1,
      1, 1, 1, 0,
      0, 1, 1, 0,
      0, 1, 0, 0,
      0, 0, 0, 1
    ),
    nrow = 5,
    byrow = TRUE,
    dimnames = list(
      model_df$hex_id,
      c("sp_equal", "sp_high", "sp_missing_trait", "sp_low")
    )
  )
  traits <- data.frame(
    species1 = c("sp_high", "sp_equal", "sp_low"),
    mass = c(20, 10, 30),
    hand_wing_index = c(2, 1, 3)
  )

  list(model_df = model_df, species_matrix = species_matrix, traits = traits)
}

test_that("GLLVM preparation filters complete sites before species occupancy", {
  source_function("assert_unique_ids.R")
  source_function("prepare_gllvm_model_data.R")

  inputs <- make_gllvm_test_inputs()
  prepared <- prepare_gllvm_model_data(
    model_df = inputs$model_df,
    species_matrix = inputs$species_matrix,
    traits_raw = inputs$traits,
    base_predictors = c("agriculture", "patch_area", "cohesion", "effort"),
    occupancy_threshold = 0.5,
    noisy_miner_species = "Manorina melanocephala"
  )

  expect_identical(rownames(prepared$Y), paste0("hex_", 1:4))
  expect_identical(colnames(prepared$Y), c("sp_equal", "sp_high"))
  expect_equal(colMeans(prepared$Y)[["sp_equal"]], 0.5)
  expect_identical(prepared$missing_traits, "sp_missing_trait")
  expect_identical(rownames(prepared$trait_matrix), colnames(prepared$Y))
  expect_false(anyNA(prepared$trait_matrix))
  expect_true(all(prepared$species_counts$n_presence >= 2L))
  expect_equal(prepared$minimum_required_presences, 2L)
  expect_equal(prepared$dimensions$n_sites, 4L)
})

test_that("GLLVM preparation rejects invalid traits and duplicate site IDs", {
  source_function("assert_unique_ids.R")
  source_function("prepare_gllvm_model_data.R")

  inputs <- make_gllvm_test_inputs()
  invalid_traits <- inputs$traits
  invalid_traits$mass[invalid_traits$species1 == "sp_equal"] <- 0
  expect_error(
    prepare_gllvm_model_data(
      inputs$model_df,
      inputs$species_matrix,
      invalid_traits,
      c("agriculture", "patch_area", "cohesion", "effort"),
      0.5,
      "Manorina melanocephala"
    ),
    "Mass must be positive"
  )

  duplicated_sites <- inputs$model_df
  duplicated_sites$hex_id[2] <- duplicated_sites$hex_id[1]
  expect_error(
    prepare_gllvm_model_data(
      duplicated_sites,
      inputs$species_matrix,
      inputs$traits,
      c("agriculture", "patch_area", "cohesion", "effort"),
      0.5,
      "Manorina melanocephala"
    ),
    "duplicated"
  )
})

test_that("GLLVM formulas retain environmental, trait, and selected spatial terms", {
  source_function("define_gllvm_formulas.R")

  formulas <- define_gllvm_formulas(c("MEM1", "MEM3"))
  environment_terms <- all.vars(formulas$environment)
  trait_terms <- all.vars(formulas$fourth_corner)

  expect_setequal(
    environment_terms,
    c(
      "pland_agriculture_2500m", "area_mn_woodland_2500m",
      "cohesion_woodland_2500m", "log_effort", "MEM1", "MEM3"
    )
  )
  expect_true(all(c("log_mass", "HWI", "MEM1", "MEM3") %in% trait_terms))
  expect_match(
    paste(deparse(formulas$fourth_corner), collapse = " "),
    "pland_agriculture_2500m.*log_mass"
  )
  expect_error(define_gllvm_formulas(c("MEM1", "MEM1")), "unique")
})

test_that("ecological group validation requires one group per modelled species", {
  source_function("validate_gllvm_species_groups.R")

  key <- data.frame(
    species = c("sp_a", "sp_b", "sp_c"),
    ecological_group = c("woodland", "edge", "introduced")
  )
  expect_no_error(
    validate_gllvm_species_groups(c("sp_a", "sp_b"), key)
  )
  expect_error(
    validate_gllvm_species_groups(c("sp_a", "missing"), key),
    "missing from"
  )
  expect_error(
    validate_gllvm_species_groups(
      c("sp_a"),
      rbind(key, data.frame(species = "sp_a", ecological_group = "edge"))
    ),
    "more than one"
  )
})
