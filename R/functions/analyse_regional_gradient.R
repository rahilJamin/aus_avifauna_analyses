import::from(dplyr, arrange, bind_rows, case_when, desc, filter, group_by,
             if_else, left_join, mutate, n, pull, select, summarise, transmute)
import::from(magrittr, `%>%`)
import::from(stats, quantile, sd)
import::from(tibble, as_tibble, tibble)
import::from(tidyr, pivot_wider)
import::from(utils, combn)

#' Analyse one regional convergence gradient
#'
#' Compares between-subregion beta diversity at the low and high endpoint of one
#' gradient. This asks whether bird assemblages from different IBRA7 subregions
#' are more similar in intact versus disturbed endpoint habitats. Retained IBRA7
#' subregions must have enough sites at both endpoints. Uncertainty is estimated
#' by bootstrapping sites within each retained subregion and endpoint.
#'
#' @param gradient_definition One-row gradient definition from
#'   `define_regional_gradients()`.
#' @param site_data Aligned regional site table with `ibra_subregion` and
#'   `site_row` columns.
#' @param species_matrix Binary site-by-species matrix aligned to `site_data`.
#' @param component_matrices List of Jaccard component matrices.
#' @param gradient_index Index used to offset the bootstrap seed.
#' @param minimum_sites_per_endpoint Minimum low-end and high-end sites required
#'   per retained subregion.
#' @param minimum_retained_subregions Minimum retained subregions required for
#'   the gradient.
#' @param bootstrap_draws Number of site-level bootstrap draws.
#' @param bootstrap_seed Base bootstrap seed.
#'
#' @return A list of regional support, pair means, plotting points, endpoint
#'   summaries, endpoint contrasts, bootstrap results, and dimensions.

analyse_regional_gradient <- function(gradient_definition,
                                      site_data,
                                      species_matrix,
                                      component_matrices,
                                      gradient_index,
                                      minimum_sites_per_endpoint,
                                      minimum_retained_subregions,
                                      bootstrap_draws,
                                      bootstrap_seed) {
  required_definition_columns <- c(
    "gradient", "gradient_label", "predictor", "low_threshold",
    "high_threshold", "threshold_rule", "low_endpoint_definition",
    "high_endpoint_definition"
  )
  if (!is.data.frame(gradient_definition) || nrow(gradient_definition) != 1L ||
      !all(required_definition_columns %in% names(gradient_definition)) ||
      anyNA(gradient_definition[required_definition_columns])) {
    stop("A complete one-row regional gradient definition is required.", call. = FALSE)
  }

  current_gradient <- as.character(gradient_definition$gradient)
  current_label <- as.character(gradient_definition$gradient_label)
  current_predictor <- as.character(gradient_definition$predictor)
  current_low_threshold <- gradient_definition$low_threshold
  current_high_threshold <- gradient_definition$high_threshold

  required_site_columns <- c("hex_id", "ibra_subregion", "site_row", current_predictor)
  missing_site_columns <- setdiff(required_site_columns, names(site_data))
  if (length(missing_site_columns) > 0L) {
    stop(
      "The regional site table lacks: ",
      paste(missing_site_columns, collapse = ", "),
      call. = FALSE
    )
  }
  assert_unique_ids(site_data, "hex_id", "regional site table")

  species_matrix <- as.matrix(species_matrix)
  if (is.null(rownames(species_matrix)) || is.null(colnames(species_matrix)) ||
      anyDuplicated(rownames(species_matrix)) ||
      anyDuplicated(colnames(species_matrix)) || anyNA(species_matrix) ||
      !is.numeric(species_matrix) || !all(species_matrix %in% c(0, 1)) ||
      !identical(rownames(species_matrix), as.character(site_data$hex_id)) ||
      !identical(site_data$site_row, seq_len(nrow(site_data)))) {
    stop("Regional site rows and the binary species matrix are misaligned.", call. = FALSE)
  }
  if (!is.numeric(site_data[[current_predictor]]) ||
      anyNA(site_data[[current_predictor]]) ||
      any(!is.finite(site_data[[current_predictor]])) ||
      anyNA(site_data$ibra_subregion)) {
    stop("Regional predictors and subregion assignments must be complete.", call. = FALSE)
  }
  if (!is.numeric(current_low_threshold) ||
      !is.numeric(current_high_threshold) ||
      length(current_low_threshold) != 1L || length(current_high_threshold) != 1L ||
      is.na(current_low_threshold) || is.na(current_high_threshold) ||
      !is.finite(current_low_threshold) || !is.finite(current_high_threshold) ||
      current_low_threshold >= current_high_threshold) {
    stop("Regional endpoint thresholds must be finite and ordered.", call. = FALSE)
  }

  valid_positive_whole <- function(value) {
    is.numeric(value) && length(value) == 1L && !is.na(value) &&
      is.finite(value) && value >= 1 && value == floor(value)
  }
  if (!valid_positive_whole(gradient_index) ||
      !valid_positive_whole(minimum_sites_per_endpoint) ||
      !valid_positive_whole(minimum_retained_subregions) ||
      minimum_retained_subregions < 2L ||
      !valid_positive_whole(bootstrap_draws) || bootstrap_draws < 2L ||
      !valid_positive_whole(bootstrap_seed)) {
    stop("Regional support, bootstrap, or seed settings are invalid.", call. = FALSE)
  }

  if (!is.list(component_matrices) || length(component_matrices) < 1L ||
      is.null(names(component_matrices)) || any(!nzchar(names(component_matrices))) ||
      anyDuplicated(names(component_matrices))) {
    stop("Regional component matrices must be a uniquely named list.", call. = FALSE)
  }
  valid_component <- vapply(
    component_matrices,
    function(component_matrix) {
      is.matrix(component_matrix) && is.numeric(component_matrix) &&
        identical(dim(component_matrix), c(nrow(species_matrix), nrow(species_matrix))) &&
        identical(rownames(component_matrix), rownames(species_matrix)) &&
        identical(colnames(component_matrix), rownames(species_matrix)) &&
        !anyNA(component_matrix) && all(is.finite(component_matrix)) &&
        isTRUE(all.equal(component_matrix, t(component_matrix), tolerance = 1e-12))
    },
    logical(1)
  )
  if (!all(valid_component)) {
    stop("Regional component matrices are malformed or misaligned.", call. = FALSE)
  }

  message("Analysing ", current_label)

  # Endpoint assignment is gradient-specific. Agriculture uses fixed strict
  # thresholds; the other gradients use inclusive study-wide quantile tails.
  current_site_data <- site_data %>%
    mutate(
      predictor_value = .data[[current_predictor]],
      endpoint = case_when(
        current_gradient == "agriculture" &
          predictor_value < current_low_threshold ~ "Low end",
        current_gradient == "agriculture" &
          predictor_value > current_high_threshold ~ "High end",
        current_gradient != "agriculture" &
          predictor_value <= current_low_threshold ~ "Low end",
        current_gradient != "agriculture" &
          predictor_value >= current_high_threshold ~ "High end",
        TRUE ~ NA_character_
      )
    )

  current_support <- current_site_data %>%
    group_by(ibra_subregion) %>%
    summarise(
      n_sites_all = n(),
      n_sites_low = sum(endpoint == "Low end", na.rm = TRUE),
      n_sites_high = sum(endpoint == "High end", na.rm = TRUE),
      minimum_predictor_value = min(predictor_value),
      maximum_predictor_value = max(predictor_value),
      .groups = "drop"
    ) %>%
    mutate(
      gradient = current_gradient,
      gradient_label = current_label,
      predictor = current_predictor,
      low_threshold = current_low_threshold,
      high_threshold = current_high_threshold,
      threshold_rule = gradient_definition$threshold_rule,
      retained =
        n_sites_low >= minimum_sites_per_endpoint &
        n_sites_high >= minimum_sites_per_endpoint,
      .before = ibra_subregion
    ) %>%
    arrange(desc(retained), ibra_subregion)

  retained_subregions <- current_support %>%
    filter(retained) %>%
    pull(ibra_subregion)

  if (length(retained_subregions) < minimum_retained_subregions) {
    print(current_support, n = Inf)
    stop(
      "Fewer than ",
      minimum_retained_subregions,
      " IBRA subregions contain at least ",
      minimum_sites_per_endpoint,
      " sites at both ends of ",
      current_label,
      ".",
      call. = FALSE
    )
  }

  current_endpoint_sites <- current_site_data %>%
    filter(
      ibra_subregion %in% retained_subregions,
      !is.na(endpoint)
    ) %>%
    mutate(
      gradient = current_gradient,
      gradient_label = current_label,
      predictor = current_predictor,
      low_threshold = current_low_threshold,
      high_threshold = current_high_threshold,
      .before = hex_id
    )

  low_rows_by_subregion <- split(
    current_endpoint_sites$site_row[
      current_endpoint_sites$endpoint == "Low end"
    ],
    current_endpoint_sites$ibra_subregion[
      current_endpoint_sites$endpoint == "Low end"
    ]
  )

  high_rows_by_subregion <- split(
    current_endpoint_sites$site_row[
      current_endpoint_sites$endpoint == "High end"
    ],
    current_endpoint_sites$ibra_subregion[
      current_endpoint_sites$endpoint == "High end"
    ]
  )

  low_rows_by_subregion <- low_rows_by_subregion[retained_subregions]
  high_rows_by_subregion <- high_rows_by_subregion[retained_subregions]

  subregion_pairs <- as_tibble(
    t(combn(retained_subregions, 2)),
    .name_repair = "minimal"
  )
  names(subregion_pairs) <- c("subregion_1", "subregion_2")

  # Calculate the directly observed mean dissimilarity for every pair of IBRA
  # subregions, separately for each endpoint and beta-diversity component.
  current_pair_rows <- list()
  output_row <- 1L

  for (component_name in names(component_matrices)) {
    current_matrix <- component_matrices[[component_name]]

    for (endpoint_name in c("Low end", "High end")) {
      if (endpoint_name == "Low end") {
        rows_by_subregion <- low_rows_by_subregion
      } else {
        rows_by_subregion <- high_rows_by_subregion
      }

      for (pair_index in seq_len(nrow(subregion_pairs))) {
        subregion_1 <- subregion_pairs$subregion_1[pair_index]
        subregion_2 <- subregion_pairs$subregion_2[pair_index]

        rows_1 <- rows_by_subregion[[subregion_1]]
        rows_2 <- rows_by_subregion[[subregion_2]]
        site_pair_values <- as.vector(current_matrix[rows_1, rows_2])

        current_pair_rows[[output_row]] <- tibble(
          gradient = current_gradient,
          gradient_label = current_label,
          endpoint = endpoint_name,
          component = component_name,
          subregion_1,
          subregion_2,
          n_sites_subregion_1 = length(rows_1),
          n_sites_subregion_2 = length(rows_2),
          n_site_pairs = length(site_pair_values),
          mean_dissimilarity = mean(site_pair_values),
          sd_dissimilarity = sd(site_pair_values)
        )

        output_row <- output_row + 1L
      }
    }
  }

  current_pair_means <- bind_rows(current_pair_rows) %>%
    arrange(component, endpoint, subregion_1, subregion_2)

  # One plotting point is one focal subregion's mean dissimilarity to every
  # other retained subregion. Every other subregion receives equal weight.
  current_subregion_points <- bind_rows(
    current_pair_means %>%
      transmute(
        gradient,
        gradient_label,
        endpoint,
        component,
        ibra_subregion = subregion_1,
        other_subregion = subregion_2,
        pair_mean = mean_dissimilarity
      ),
    current_pair_means %>%
      transmute(
        gradient,
        gradient_label,
        endpoint,
        component,
        ibra_subregion = subregion_2,
        other_subregion = subregion_1,
        pair_mean = mean_dissimilarity
      )
  ) %>%
    group_by(
      gradient,
      gradient_label,
      endpoint,
      component,
      ibra_subregion
    ) %>%
    summarise(
      n_other_subregions = n(),
      mean_dissimilarity_to_other_subregions = mean(pair_mean),
      sd_among_subregion_pairs = sd(pair_mean),
      .groups = "drop"
    ) %>%
    arrange(component, ibra_subregion, endpoint)

  current_subregion_contrasts <- current_subregion_points %>%
    select(
      gradient,
      gradient_label,
      component,
      ibra_subregion,
      endpoint,
      mean_dissimilarity_to_other_subregions
    ) %>%
    pivot_wider(
      names_from = endpoint,
      values_from = mean_dissimilarity_to_other_subregions
    ) %>%
    transmute(
      gradient,
      gradient_label,
      component,
      ibra_subregion,
      low_end_mean = `Low end`,
      high_end_mean = `High end`,
      high_minus_low = high_end_mean - low_end_mean
    ) %>%
    arrange(component, ibra_subregion)

  # Bootstrap sites within each retained subregion and endpoint. This preserves
  # the regional/end-point design while estimating uncertainty from site-level
  # sampling variation.
  set.seed(bootstrap_seed + gradient_index - 1L)
  current_bootstrap_rows <- vector("list", bootstrap_draws)

  for (bootstrap_index in seq_len(bootstrap_draws)) {
    sampled_low_rows <- lapply(
      low_rows_by_subregion,
      function(rows) sample(rows, length(rows), replace = TRUE)
    )
    sampled_high_rows <- lapply(
      high_rows_by_subregion,
      function(rows) sample(rows, length(rows), replace = TRUE)
    )

    current_draw_rows <- list()
    current_draw_row <- 1L

    for (component_name in names(component_matrices)) {
      current_matrix <- component_matrices[[component_name]]

      for (endpoint_name in c("Low end", "High end")) {
        if (endpoint_name == "Low end") {
          sampled_rows_by_subregion <- sampled_low_rows
        } else {
          sampled_rows_by_subregion <- sampled_high_rows
        }

        pair_means <- rep(NA_real_, nrow(subregion_pairs))

        for (pair_index in seq_len(nrow(subregion_pairs))) {
          subregion_1 <- subregion_pairs$subregion_1[pair_index]
          subregion_2 <- subregion_pairs$subregion_2[pair_index]

          rows_1 <- sampled_rows_by_subregion[[subregion_1]]
          rows_2 <- sampled_rows_by_subregion[[subregion_2]]

          pair_means[pair_index] <- mean(current_matrix[rows_1, rows_2])
        }

        current_draw_rows[[current_draw_row]] <- tibble(
          gradient = current_gradient,
          gradient_label = current_label,
          bootstrap_draw = bootstrap_index,
          endpoint = endpoint_name,
          component = component_name,
          estimate = mean(pair_means)
        )

        current_draw_row <- current_draw_row + 1L
      }
    }

    current_bootstrap_rows[[bootstrap_index]] <- bind_rows(current_draw_rows)

    if (bootstrap_index %% 100L == 0L) {
      message(
        current_label,
        ": bootstrap draw ",
        bootstrap_index,
        " of ",
        bootstrap_draws
      )
    }
  }

  current_bootstrap_results <- bind_rows(current_bootstrap_rows)

  observed_endpoint_summary <- current_pair_means %>%
    group_by(gradient, gradient_label, endpoint, component) %>%
    summarise(
      estimate = mean(mean_dissimilarity),
      sd_among_subregion_pair_means = sd(mean_dissimilarity),
      n_subregions = length(retained_subregions),
      n_subregion_pairs = n(),
      .groups = "drop"
    )

  bootstrap_endpoint_summary <- current_bootstrap_results %>%
    group_by(gradient, gradient_label, endpoint, component) %>%
    summarise(
      bootstrap_se = sd(estimate),
      lower_95 = quantile(estimate, 0.025),
      upper_95 = quantile(estimate, 0.975),
      .groups = "drop"
    )

  current_endpoint_summary <- observed_endpoint_summary %>%
    left_join(
      bootstrap_endpoint_summary,
      by = c("gradient", "gradient_label", "endpoint", "component")
    ) %>%
    mutate(
      endpoint_definition = if_else(
        endpoint == "Low end",
        gradient_definition$low_endpoint_definition,
        gradient_definition$high_endpoint_definition
      ),
      low_threshold = current_low_threshold,
      high_threshold = current_high_threshold,
      .after = endpoint
    ) %>%
    arrange(component, endpoint)

  observed_contrasts <- current_endpoint_summary %>%
    select(gradient, gradient_label, component, endpoint, estimate) %>%
    pivot_wider(names_from = endpoint, values_from = estimate) %>%
    transmute(
      gradient,
      gradient_label,
      component,
      low_end_mean = `Low end`,
      high_end_mean = `High end`,
      high_minus_low = high_end_mean - low_end_mean
    )

  bootstrap_contrasts <- current_bootstrap_results %>%
    select(
      gradient,
      gradient_label,
      bootstrap_draw,
      component,
      endpoint,
      estimate
    ) %>%
    pivot_wider(names_from = endpoint, values_from = estimate) %>%
    mutate(high_minus_low = `High end` - `Low end`)

  current_contrast_rows <- vector("list", length(component_matrices))

  for (component_index in seq_along(component_matrices)) {
    component_name <- names(component_matrices)[component_index]

    observed_change <- observed_contrasts$high_minus_low[
      observed_contrasts$component == component_name
    ]
    bootstrap_change <- bootstrap_contrasts$high_minus_low[
      bootstrap_contrasts$component == component_name
    ]
    centred_bootstrap_error <- bootstrap_change - mean(bootstrap_change)

    current_contrast_rows[[component_index]] <- tibble(
      component = component_name,
      difference_bootstrap_se = sd(bootstrap_change),
      difference_lower_95 = quantile(bootstrap_change, 0.025),
      difference_upper_95 = quantile(bootstrap_change, 0.975),
      bootstrap_two_sided_p_value = (
        sum(abs(centred_bootstrap_error) >= abs(observed_change)) + 1
      ) / (bootstrap_draws + 1)
    )
  }

  current_endpoint_contrasts <- observed_contrasts %>%
    left_join(bind_rows(current_contrast_rows), by = "component") %>%
    mutate(
      low_threshold = current_low_threshold,
      high_threshold = current_high_threshold,
      interpretation = case_when(
        difference_upper_95 < 0 ~
          "Assemblages are more similar at the high end",
        difference_lower_95 > 0 ~
          "Assemblages are more similar at the low end",
        TRUE ~ "No clear endpoint difference"
      ),
      .after = gradient_label
    )

  dimensions <- tibble(
    gradient = current_gradient,
    gradient_label = current_label,
    n_sites_all = nrow(species_matrix),
    n_species = ncol(species_matrix),
    n_endpoint_sites = nrow(current_endpoint_sites),
    n_retained_subregions = length(retained_subregions),
    n_subregion_pairs = nrow(subregion_pairs)
  )

  list(
    regional_support = current_support,
    subregion_pair_means = current_pair_means,
    subregion_plot_points = current_subregion_points,
    subregion_contrasts = current_subregion_contrasts,
    endpoint_summary = current_endpoint_summary,
    endpoint_contrasts = current_endpoint_contrasts,
    bootstrap_results = current_bootstrap_results,
    dimensions = dimensions
  )
}
