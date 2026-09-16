#' Produce GLLVM tables, diagnostic plots, and publication figures
#'
#' Reads the completed model checkpoint and reproduces the validated
#' publication tables and plots. Plot data, scales, colours, and export
#' settings are carried over from the original output script unchanged.
#'
#' @param cfg Main workflow configuration with input and output paths.
#' @return Paths to the generated files, for a `targets` file target.

import::from(magrittr, `%>%`)
import::from(dplyr, arrange, bind_rows, case_when, count, desc, distinct,
             filter, group_by, if_else, inner_join, left_join, mutate,
             n_distinct, ntile, pull, recode, relocate, rename, select,
             slice_head, summarise, transmute, ungroup)
import::from(forcats, fct_rev)
import::from(stringr, str_detect)
import::from(tibble, as_tibble, rownames_to_column, tibble)
import::from(tidyr, expand_grid, replace_na, separate)
import::from(stats, coef, quantile)
import::from(grid, unit)
import::from(scales, percent)
import::from(rcartocolor, carto_pal)
import::from(cowplot, get_legend)
import::from(patchwork, wrap_plots)
import::from(gllvm, gllvm)
import::from(graphics, par, plot)
import::from(grDevices, dev.off, pdf)
import::from(ggplot2, aes, coord_fixed, element_blank, element_line,
             element_rect, element_text, facet_grid, facet_wrap, geom_col,
             geom_errorbarh, geom_point, geom_text, geom_tile, geom_vline,
             ggplot, ggsave, guide_colourbar, guide_legend, guides, labs,
             margin, scale_alpha_manual, scale_color_manual,
             scale_fill_gradientn, scale_fill_manual, scale_x_continuous,
             scale_x_discrete, scale_y_discrete, theme, theme_bw,
             theme_minimal, vars)

produce_gllvm_outputs <- function(cfg) {
  message("Producing GLLVM/JSDM outputs")

  metadata <- read_rds_checked(
    cfg$gllvm_metadata_path,
    "GLLVM metadata",
    required_class = "list",
    required_names = c(
      "dimensions", "species_counts", "traits", "species_order",
      "missing_traits", "minimum_required_presences", "occupancy_threshold"
    )
  )
  fit_env <- read_rds_checked(cfg$gllvm_env_path, "environment GLLVM")
  fit_fourth <- read_rds_checked(cfg$gllvm_fourth_corner_path, "fourth-corner GLLVM")

  # Export model summary tables ----

  species_counts <- metadata$species_counts
  traits <- metadata$traits
  species_order <- metadata$species_order
  missing_traits <- metadata$missing_traits
  minimum_required_presences <- metadata$minimum_required_presences

  # Preserve a species-exclusion record because trait matching changes both
  # models' species sets. An empty file clears exclusions from an earlier run.
  write_lines_quiet(
    missing_traits,
    file.path(cfg$table_dir, "06_gllvm_species_missing_traits.txt")
  )

  write_csv_quiet(
    metadata$dimensions,
    file.path(cfg$table_dir, "06_gllvm_model_dimensions_with_dbmem.csv")
  )

  beta_matrix_species <- coef(fit_env)$Xcoef %>%
    as.data.frame() %>%
    rownames_to_column("species") %>%
    as_tibble() %>%
    arrange(match(species, species_order)) %>%
    left_join(species_counts, by = "species") %>%
    left_join(traits, by = "species") %>%
    relocate(n_presence, n_absence, occupancy, .after = species)

  stopifnot(setequal(beta_matrix_species$species, species_order))
  stopifnot(all(beta_matrix_species$n_presence >= minimum_required_presences))
  stopifnot(all(beta_matrix_species$occupancy >= metadata$occupancy_threshold))

  # One species coefficient table also contains final occurrence counts,
  # occupancy and traits, avoiding separate copies of these quantities.
  write_csv_quiet(
    beta_matrix_species,
    file.path(cfg$table_dir, "06_gllvm_species_coefficients.csv")
  )

  fourth_corner_coefficients <- coef(fit_fourth)$B %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    as_tibble()

  write_csv_quiet(
    fourth_corner_coefficients,
    file.path(cfg$table_dir, "06_fourth_corner_coefficients.csv")
  )

  # Render basic model diagnostics ----

  assert_not_raw_output_path(cfg$diagnostic_dir)
  dir.create(cfg$diagnostic_dir, recursive = TRUE, showWarnings = FALSE)

  pdf(
    file.path(cfg$diagnostic_dir, "06_gllvm_basic_diagnostics.pdf"),
    width = 10,
    height = 8
  )

  par(mfrow = c(2, 2))
  plot(fit_env)

  par(mfrow = c(2, 2))
  plot(fit_fourth)

  dev.off()

  # Render publication figures ----

  figure_dir <- file.path(cfg$figure_dir, "gllvm_plots")

  assert_not_raw_output_path(figure_dir)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  message("Preparing GLLVM coefficient and trait figures")


  # -----------------------------------------------------------------------------
  # 3. Define the three environmental predictors to plot
  # -----------------------------------------------------------------------------

  focal_covariates <- c(
    "pland_agriculture_2500m",
    "cohesion_woodland_2500m",
    "area_mn_woodland_2500m"
  )

  pretty_cov_names <- c(
    pland_agriculture_2500m = "Agriculture cover",
    cohesion_woodland_2500m = "Woodland cohesion",
    area_mn_woodland_2500m  = "Mean woodland patch area"
  )

  heatmap_cov_order <- c(
    "Agriculture cover",
    "Woodland cohesion",
    "Mean woodland patch area"
  )

  # -----------------------------------------------------------------------------
  # 4. Extract coefficients directly from summary(fit_env)
  # -----------------------------------------------------------------------------
  # Estimates, standard errors and significance markers come from the same
  # fitted summary. Coefficients remain on the probit scale.

  coef_table <- summary(fit_env)$Coef.tableX %>%
    as.data.frame(check.names = FALSE) %>%
    rownames_to_column("term") %>%
    as_tibble() %>%
    transmute(
      term,
      beta = Estimate,
      se = `Std. Error`,
      z_value = `z value`,
      p_value = `Pr(>|z|)`
    )

  # -----------------------------------------------------------------------------
  # 5. Keep only the three predictor-by-species coefficients
  # -----------------------------------------------------------------------------
  # The term names are expected to look like this:
  # pland_agriculture_2500m:Acanthiza ewingii

  alpha <- 0.05

  heatmap_df <- coef_table %>%
    separate(
      term,
      into = c("covariate", "species"),
      sep = ":",
      remove = FALSE,
      extra = "merge"
    ) %>%
    filter(covariate %in% focal_covariates) %>%
    mutate(
      sig = case_when(
        p_value < alpha & beta > 0 ~ "Positive",
        p_value < alpha & beta < 0 ~ "Negative",
        TRUE ~ "NS"
      ),
      star = if_else(p_value < alpha, "*", ""),
      predictor = factor(
        pretty_cov_names[covariate],
        levels = rev(heatmap_cov_order)
      )
    )

  # -----------------------------------------------------------------------------
  # 6. Cap beta values for plotting
  # -----------------------------------------------------------------------------

  # The symmetric colour limit uses the existing 97th percentile of absolute
  # coefficients. This clips colours only; coefficients and tests are unchanged.
  beta_cap <- quantile(abs(heatmap_df$beta), 0.97, na.rm = TRUE)

  heatmap_df <- heatmap_df %>%
    mutate(
      beta_capped = pmin(pmax(beta, -beta_cap), beta_cap)
    )

  # -----------------------------------------------------------------------------
  # 7. Select species for the grouped heatmap
  # -----------------------------------------------------------------------------

  # Preserve the author-specified ecological groups and preferred species
  # rankings. Up to five supported species per group enter the selected panel.
  woodland_dependent_species <- c(
    "Melithreptus brevirostris", "Melithreptus gularis", "Lichenostomus melanops",
    "Melithreptus validirostris", "Nesoptilotis leucotis", "Melithreptus lunatus",
    "Ptilotula fusca", "Hylacola cauta", "Hylacola pyrrhopygia", "Acanthiza apicalis",
    "Pyrrholaemus sagittatus", "Petroica boodang", "Petroica goodenovii",
    "Acanthiza uropygialis", "Aphelocephala leucopsis", "Melanodryas cucullata",
    "Eopsaltria australis", "Petroica rodinogaster", "Melanodryas vittata",
    "Pomatostomus superciliosus", "Oreoica gutturalis", "Daphoenositta chrysoptera",
    "Falcunculus frontatus", "Pomatostomus temporalis", "Cinclosoma punctatum",
    "Calyptorhynchus lathami", "Polytelis swainsonii", "Lathamus discolor",
    "Neophema pulchella", "Ptilotula ornata", "Climacteris picumnus",
    "Climacteris erythrops", "Stagonopleura guttata", "Stagonopleura bella",
    "Acanthorhynchus tenuirostris", "Caligavis chrysops", "Melithreptus albogularis",
    "Phylidonyris pyrrhopterus", "Myzomela sanguinolenta", "Plectorhyncha lanceolata",
    "Sericornis frontalis", "Sericornis humilis", "Neosericornis citreogularis",
    "Gerygone olivacea", "Gerygone fusca", "Acanthiza ewingii", "Gerygone mouki",
    "Acanthiza lineata", "Acanthiza nana", "Acanthiza reguloides", "Petroica phoenicea",
    "Pachycephala olivacea", "Pachycephala pectoralis", "Parvipsitta porphyrocephala",
    "Glossopsitta concinna", "Parvipsitta pusilla", "Platycercus elegans",
    "Alisterus scapularis", "Aprosmictus erythropterus", "Platycercus caledonicus",
    "Cacomantis flabelliformis", "Chalcites basalis", "Heteroscenes pallidus",
    "Chalcites lucidus", "Chalcites osculans", "Cacomantis variolosus",
    "Scythrops novaehollandiae", "Chalcites minutillus", "Todiramphus sanctus",
    "Ceyx azureus", "Malurus lamberti", "Alectura lathami", "Myiagra rubecula",
    "Rhipidura albiscapa", "Myiagra cyanoleuca", "Rhipidura rufifrons",
    "Dicrurus bracteatus", "Oriolus sagittatus", "Coracina papuensis",
    "Edolisoma tenuirostre", "Corcorax melanorhamphos", "Dicaeum hirundinaceum",
    "Smicrornis brevirostris", "Pachycephala rufiventris", "Microeca fascinans",
    "Cormobates leucophaea", "Zanda funerea", "Callocephalon fimbriatum",
    "Myiagra inquieta", "Artamus cyanopterus", "Petroica rosea",
    "Nesoptilotis flavicollis", "Melithreptus affinis"
  )

  matrix_edge_species <- c(
    "Anthochaera carunculata", "Anthochaera chrysoptera", "Entomyzon cyanotis",
    "Anthochaera paradoxa", "Manorina melanocephala", "Phylidonyris novaehollandiae",
    "Pardalotus striatus", "Pardalotus punctatus", "Acanthiza pusilla",
    "Acanthiza chrysorrhoa", "Colluricincla harmonica", "Eolophus roseicapilla",
    "Trichoglossus moluccanus", "Cacatua sanguinea", "Cacatua galerita",
    "Cacatua tenuirostris", "Platycercus eximius", "Eudynamys orientalis",
    "Manorina flavigula", "Dacelo novaeguineae", "Malurus cyaneus",
    "Grallina cyanoleuca", "Rhipidura leucophrys", "Coracina novaehollandiae",
    "Gymnorhina tibicen", "Cracticus torquatus", "Corvus coronoides",
    "Strepera versicolor", "Cracticus nigrogularis", "Strepera graculina",
    "Strepera fuliginosa", "Corvus mellori", "Corvus tasmanicus",
    "Corvus orru", "Neochmia temporalis", "Stizoptera bichenovii",
    "Zosterops lateralis",
    "Philemon citreogularis", "Philemon corniculatus", "Barnardius zonarius",
    "Platycercus adscitus", "Centropus phasianinus", "Eurystomus orientalis",
    "Struthidea cinerea"
  )

  open_country_species <- c(
    "Epthianura albifrons", "Ptilotula penicillata", "Nymphicus hollandicus",
    "Northiella haematogaster", "Psephotus haematonotus", "Neophema elegans",
    "Neophema chrysostoma", "Acanthagenys rufogularis", "Gavicalis virescens",
    "Malurus melanocephalus", "Artamus leucorynchus", "Taeniopygia guttata",
    "Cincloramphus mathewsi"
  )

  introduced_species <- c(
    "Passer domesticus", "Alauda arvensis", "Carduelis carduelis",
    "Chloris chloris", "Sturnus vulgaris", "Turdus merula",
    "Acridotheres tristis"
  )

  group_levels <- c(
    "Woodland-dependent species",
    "Matrix/edge native species",
    "Open-country native species",
    "Introduced species"
  )

  species_group_key_all <- bind_rows(
    tibble(
      species = woodland_dependent_species,
      ecological_group = "Woodland-dependent species"
    ),
    tibble(
      species = matrix_edge_species,
      ecological_group = "Matrix/edge native species"
    ),
    tibble(
      species = open_country_species,
      ecological_group = "Open-country native species"
    ),
    tibble(
      species = introduced_species,
      ecological_group = "Introduced species"
    )
  ) %>%
    filter(species %in% unique(heatmap_df$species)) %>%
    distinct(species, ecological_group)

  # Only species in the author-defined ecological groups enter grouped summaries.
  # Other modelled species remain visible in the all-species heatmap.

  selected_species_rank <- tibble(
    species = c(
      "Climacteris erythrops",
      "Pachycephala olivacea",
      "Eopsaltria australis",
      "Petroica rosea",
      "Rhipidura rufifrons",
      
      "Manorina melanocephala",
      "Cracticus nigrogularis",
      "Grallina cyanoleuca",
      "Gymnorhina tibicen",
      "Cacatua sanguinea",

      "Psephotus haematonotus",
      "Ptilotula penicillata",
      "Taeniopygia guttata",
      "Cincloramphus mathewsi",
      "Nymphicus hollandicus",
      
      "Acridotheres tristis",
      "Sturnus vulgaris",
      "Passer domesticus",
      "Carduelis carduelis",
      "Turdus merula"
    ),
    selection_rank = rep(1:5, times = 4)
  )

  species_group_key <- heatmap_df %>%
    inner_join(species_group_key_all, by = "species") %>%
    filter(sig != "NS") %>%
    distinct(species, ecological_group) %>%
    left_join(selected_species_rank, by = "species") %>%
    mutate(
      selection_rank = replace_na(selection_rank, 999),
      ecological_group = factor(ecological_group, levels = group_levels)
    ) %>%
    group_by(species, ecological_group, selection_rank) %>%
    summarise(
      n_supported = sum(heatmap_df$sig[match(species, heatmap_df$species)] != "NS", na.rm = TRUE),
      .groups = "drop"
    ) %>%
    select(species, ecological_group, selection_rank) %>%
    inner_join(
      heatmap_df %>%
        group_by(species) %>%
        summarise(
          n_supported = sum(sig != "NS"),
          max_abs_beta = max(abs(beta), na.rm = TRUE),
          .groups = "drop"
        ),
      by = "species"
    ) %>%
    group_by(ecological_group) %>%
    arrange(selection_rank, desc(n_supported), desc(max_abs_beta), species, .by_group = TRUE) %>%
    slice_head(n = 5) %>%
    ungroup() %>%
    select(species, ecological_group) %>%
    mutate(
      ecological_group = factor(ecological_group, levels = group_levels)
    )

  species_order <- heatmap_df %>%
    inner_join(species_group_key, by = "species") %>%
    mutate(
      response_component = case_when(
        ecological_group == "Woodland-dependent species" &
          covariate == "cohesion_woodland_2500m" ~ beta,
        ecological_group == "Woodland-dependent species" &
          covariate == "area_mn_woodland_2500m" ~ beta,
        ecological_group == "Woodland-dependent species" &
          covariate == "pland_agriculture_2500m" ~ -beta,
        
        ecological_group == "Matrix/edge native species" &
          covariate == "pland_agriculture_2500m" ~ beta,
        ecological_group == "Matrix/edge native species" &
          covariate == "cohesion_woodland_2500m" ~ -beta,
        ecological_group == "Matrix/edge native species" &
          covariate == "area_mn_woodland_2500m" ~ -beta,
        
        ecological_group == "Open-country native species" &
          covariate == "pland_agriculture_2500m" ~ beta,
        ecological_group == "Open-country native species" &
          covariate == "cohesion_woodland_2500m" ~ -beta,
        ecological_group == "Open-country native species" &
          covariate == "area_mn_woodland_2500m" ~ -beta,
        
        ecological_group == "Introduced species" &
          covariate == "pland_agriculture_2500m" ~ beta,
        ecological_group == "Introduced species" &
          covariate == "cohesion_woodland_2500m" ~ -beta,
        ecological_group == "Introduced species" &
          covariate == "area_mn_woodland_2500m" ~ -beta,
        
        TRUE ~ 0
      ),
      response_component = if_else(sig == "NS", 0, response_component)
    ) %>%
    group_by(species, ecological_group) %>%
    summarise(
      n_supported = sum(sig != "NS"),
      response_score = sum(response_component, na.rm = TRUE),
      max_abs_beta = max(abs(beta), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(
      ecological_group,
      desc(n_supported),
      desc(response_score),
      desc(max_abs_beta),
      species
    ) %>%
    pull(species)

  heatmap_selected_df <- heatmap_df %>%
    inner_join(species_group_key, by = "species") %>%
    mutate(
      ecological_group = factor(ecological_group, levels = group_levels),
      species = factor(species, levels = species_order)
    )

  # -----------------------------------------------------------------------------
  # 9. Plot selected species heatmap with group labels
  # -----------------------------------------------------------------------------

  p_heatmap_selected_with_legend <- ggplot(
    heatmap_selected_df,
    aes(x = species, y = predictor, fill = beta_capped)
  ) +
    geom_tile(
      colour = "grey90",
      linewidth = 0.025
    ) +
    geom_text(
      aes(label = star),
      colour = "grey10",
      size = 7.25,
      vjust = 0.78
    ) +
    facet_grid(
      . ~ ecological_group,
      scales = "free_x",
      space = "free_x"
    ) +
    scale_fill_gradientn(
      colours = rev(carto_pal(7, "Geyser")),
      limits = c(-beta_cap, beta_cap),
      breaks = c(-beta_cap, 0, beta_cap),
      labels = c(
        sprintf("\u2264 %.2f", -beta_cap),
        "0",
        sprintf("\u2265 %.2f", beta_cap)
      ),
      name = expression(beta~"(probit)"),
      guide = guide_colourbar(
        title.position = "top",
        title.hjust = 0.5,
        direction = "vertical",
        barwidth = unit(0.35, "cm"),
        barheight = unit(4.1, "cm"),
        ticks.colour = "black",
        frame.colour = "black"
      )
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(
      strip.background = element_rect(
        fill = NULL,
        colour = NULL,
        linewidth = NULL
      ),
      strip.text.x = element_blank(),
      panel.spacing.x = unit(0.7, "lines"),
      panel.grid = element_blank(),
      axis.text.x = element_text(
        face = "italic",
        size = 12,
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        colour = "grey20"
      ),
      axis.text.y = element_text(
        size = 12,
        colour = "grey20"
      ),
      axis.ticks = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8.5),
      plot.margin = margin(4, 8, 4, 8)
    )

  heatmap_legend_grob <- get_legend(p_heatmap_selected_with_legend)

  p_heatmap_selected <- p_heatmap_selected_with_legend + theme(legend.position = "none")


  # -----------------------------------------------------------------------------
  # 10. Prepare all species heatmap in four panels
  # -----------------------------------------------------------------------------

  all_species_order <- sort(unique(heatmap_df$species))

  species_panel_key <- tibble(
    species = all_species_order,
    species_number = seq_along(all_species_order),
    panel = paste0("Panel ", ntile(species_number, 4))
  )

  plot_list <- list()

  for (i in 1:4) {
    
    panel_name <- paste0("Panel ", i)
    
    panel_species <- species_panel_key %>%
      filter(panel == panel_name) %>%
      pull(species)
    
    panel_df <- heatmap_df %>%
      filter(species %in% panel_species) %>%
      mutate(
        species = factor(species, levels = panel_species)
      )
    
    plot_list[[i]] <- ggplot(
      panel_df,
      aes(x = species, y = predictor, fill = beta_capped)
    ) +
      geom_tile(colour = "grey93", linewidth = 0.05) +
      geom_text(
        aes(label = star),
        colour = "grey10",
        size = 5.5,
        vjust = 0.78
      ) +
      scale_fill_gradientn(
        colours = rev(carto_pal(7, "Geyser")),
        limits = c(-beta_cap, beta_cap),
        breaks = c(-beta_cap, 0, beta_cap),
        labels = c(
          sprintf("\u2264 %.2f", -beta_cap),
          "0",
          sprintf("\u2265 %.2f", beta_cap)
        ),
        name = expression(beta~"(probit)"),
        guide = guide_colourbar(
          title.position = "top",
          title.hjust = 0.5,
          direction = "horizontal",
          barwidth = unit(7, "cm"),
          barheight = unit(0.35, "cm"),
          ticks.colour = "black",
          frame.colour = "black"
        )
      ) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      coord_fixed(ratio = 1) +
      labs(x = NULL, y = NULL) +
      theme_minimal(base_size = 9) +
      theme(
        axis.text.x = element_text(
          face = "italic",
          size = 7.5,
          angle = 90,
          hjust = 1,
          vjust = 0.5,
          colour = "grey20"
        ),
        axis.text.y = element_text(
          size = 8.,
          hjust = 1,
          colour = "grey20"
        ),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8.5),
        plot.margin = margin(4, 8, 4, 8)
      )
  }

  p_heatmap_all <- wrap_plots(
    plot_list,
    ncol = 1,
    guides = "collect"
  ) &
    theme(legend.position = "none")

  # -----------------------------------------------------------------------------
  # 12. Save plots
  # -----------------------------------------------------------------------------

  ggsave(
    filename = file.path(figure_dir, "beta_heatmap_selected_grouped.svg"),
    plot = p_heatmap_selected,
    width = 14.5,
    height = 3.5,
    dpi = 800,
    units = "in"
  )

  # Retain one shared legend for assembling the selected and full heatmaps.
  ggsave(
    filename = file.path(figure_dir, paste0("heatmap_standalone_legend.svg")), 
    plot = heatmap_legend_grob, 
    width = 6,
    height = 1.5,
    dpi = 600,
    bg = "transparent"
  )

  ggsave(
    filename = file.path(figure_dir, paste0("beta_heatmap_all.svg")),
    plot = p_heatmap_all,
    width = 13,
    height = 12,
    dpi = 800,
    units = "in"
  )

  # Extract trait interactions from the saved fourth-corner fit. The existing
  # intervals use estimate +/- 1.96 standard errors on the coefficient scale.
  fit_trait <- fit_fourth

  # 1. Extract coefficients from the gllvm summary
  sum_fit_trait <- summary(fit_trait)
  coef_df <- as.data.frame(sum_fit_trait$Coef.tableX)
  coef_df$Term <- rownames(coef_df)

  # 2. Filter, clean, and calculate CIs
  plot_data <- coef_df %>%
    # Keep only the trait-environment interactions
    filter(str_detect(Term, ":")) %>%
    rename(
      Estimate = Estimate,
      StdError = `Std. Error`
    ) %>%
    # Calculate 95% Confidence Intervals
    mutate(
      CI_lower = Estimate - (1.96 * StdError),
      CI_upper = Estimate + (1.96 * StdError),
      # Determine significance for coloring
      Significance = case_when(
        CI_lower > 0 ~ "Positive",
        CI_upper < 0 ~ "Negative",
        TRUE ~ "Non-significant"
      )
    ) %>%
    # Split the interaction term into Landscape and Trait
    separate(Term, into = c("Landscape", "Trait"), sep = ":") %>%
    # Clean up labels for publication
    mutate(
      Landscape = case_when(
        Landscape == "pland_agriculture_2500m" ~ "Agricultural Cover",
        Landscape == "area_mn_woodland_2500m" ~ "Mean Woodland Area",
        Landscape == "cohesion_woodland_2500m" ~ "Woodland Cohesion"
      ),
      Trait_Clean = case_when(
        Trait == "log_mass" ~ "Log Body Mass",
        Trait == "HWI" ~ "Hand-Wing Index"
      )
    )

  # Ensure consistent ordering of traits on the y-axis
  plot_data$Trait_Clean <- factor(plot_data$Trait_Clean, levels = c("Log Body Mass", "Hand-Wing Index"))


  # 3. Build the trait coefficient plot
  p_morphology <- ggplot(plot_data, aes(x = Estimate, y = Trait_Clean, color = Significance)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0, linewidth = 0.8) +
    geom_point(size = 4) +
    facet_wrap(~ Landscape, scales = "free_x", ncol = 3) +
    scale_color_manual(
      values = c(
        "Positive" = "#008080",       # Teal
        "Negative" = "#ca562c",       # Orange
        "Non-significant" = "grey70"  # Grey
      )
    ) +
    labs(
      x = expression("" * beta * " Coefficient (probit)"),
      y = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(
      strip.text = element_text(size = 12),
      strip.background = element_rect(fill = "grey97", color = "grey55"),
      panel.spacing = unit(1.5, "lines"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", linetype = "dotted"),
      panel.grid.major.x = element_blank(),
      axis.text.y = element_text(color = "black", size = 12),
      axis.text.x = element_text(color = "black"),
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.margin = margin(t = 10, r = 15, b = 10, l = 10)
    )

  ggsave(
    filename = file.path(figure_dir, paste0("trait_heatmap_selected.svg")),
    plot = p_morphology,
    width = 12.5,
    height = 6.,
    dpi = 800,
    units = "in"
  )



  # -----------------------------------------------------------------------------
  # 3. Build species key
  # -----------------------------------------------------------------------------

  species_group_key <- bind_rows(
    tibble(
      species = woodland_dependent_species,
      ecological_group = "Woodland-dependent species"
    ),
    tibble(
      species = matrix_edge_species,
      ecological_group = "Edge-tolerant species"
    ),
    tibble(
      species = open_country_species,
      ecological_group = "Open-country species"
    ),
    tibble(
      species = introduced_species,
      ecological_group = "Introduced species"
    )
  ) %>%
    distinct(species, ecological_group)

  # -----------------------------------------------------------------------------
  # 4. Prepare model-response table
  # -----------------------------------------------------------------------------

  group_order <- c(
    "Woodland-dependent species",
    "Edge-tolerant species",
    "Open-country species",
    "Introduced species"
  )

  direction_levels <- c("Negative", "Positive")

  heatmap_habitat_df <- heatmap_df %>%
    left_join(species_group_key, by = "species") %>%
    mutate(
      ecological_group = factor(ecological_group, levels = group_order),
      
      predictor_label = recode(covariate, !!!pretty_cov_names),
      predictor_label = factor(
        predictor_label,
        levels = heatmap_cov_order
      ),
      
      coefficient_direction = case_when(
        beta < 0 ~ "Negative",
        beta > 0 ~ "Positive",
        TRUE ~ NA_character_
      ),
      
      supported_direction = case_when(
        sig %in% c("Negative", "Positive") ~ sig,
        TRUE ~ NA_character_
      )
    )

  # -----------------------------------------------------------------------------
  # 5. Denominators
  # -----------------------------------------------------------------------------

  # Use every modelled species in each ecological group as the denominator.
  # Light bars show all coefficient signs; dark bars show supported signs
  # under the existing alpha. Both layers therefore share the same denominator.
  group_denominators <- heatmap_habitat_df %>%
    distinct(species, ecological_group) %>%
    count(ecological_group, name = "n_group_species") %>%
    mutate(ecological_group = factor(ecological_group, levels = group_order))

  # -----------------------------------------------------------------------------
  # 6. Count all coefficient signs
  # -----------------------------------------------------------------------------

  all_direction_counts <- heatmap_habitat_df %>%
    filter(!is.na(coefficient_direction)) %>%
    group_by(ecological_group, predictor_label, coefficient_direction) %>%
    summarise(
      n_all_direction = n_distinct(species),
      .groups = "drop"
    )

  # -----------------------------------------------------------------------------
  # 7. Count supported coefficient signs
  # -----------------------------------------------------------------------------

  supported_direction_counts <- heatmap_habitat_df %>%
    filter(!is.na(supported_direction)) %>%
    transmute(
      ecological_group,
      predictor_label,
      coefficient_direction = supported_direction,
      species
    ) %>%
    group_by(ecological_group, predictor_label, coefficient_direction) %>%
    summarise(
      n_supported = n_distinct(species),
      .groups = "drop"
    )
  # -----------------------------------------------------------------------------
  # 8. Combine counts
  # -----------------------------------------------------------------------------

  direction_plot_df <- expand_grid(
    ecological_group = factor(group_order, levels = group_order),
    predictor_label = factor(heatmap_cov_order, levels = heatmap_cov_order),
    coefficient_direction = factor(direction_levels, levels = direction_levels)
  ) %>%
    left_join(
      all_direction_counts,
      by = c("ecological_group", "predictor_label", "coefficient_direction")
    ) %>%
    left_join(
      supported_direction_counts,
      by = c("ecological_group", "predictor_label", "coefficient_direction")
    ) %>%
    left_join(
      group_denominators,
      by = "ecological_group"
    ) %>%
    mutate(
      n_all_direction = replace_na(n_all_direction, 0L),
      n_supported = replace_na(n_supported, 0L),
      
      proportion_all_direction = n_all_direction / n_group_species,
      proportion_supported = n_supported / n_group_species,
      
      signed_all = if_else(
        coefficient_direction == "Negative",
        -proportion_all_direction,
        proportion_all_direction
      ),
      
      signed_supported = if_else(
        coefficient_direction == "Negative",
        -proportion_supported,
        proportion_supported
      ),
      
      ecological_group_short = case_when(
        ecological_group == "Open-country native species" ~ "Open-country natives",
        TRUE ~ as.character(ecological_group)
      ),
      
      ecological_group_label = paste0(
        ecological_group_short,
        "\n(n = ",
        n_group_species,
        ")"
      )
    )

  direction_plot_df <- direction_plot_df %>%
    mutate(
      ecological_group_label = factor(
        ecological_group_label,
        levels = group_denominators %>%
          arrange(ecological_group) %>%
          mutate(
            ecological_group_short = case_when(
              ecological_group == "Open-country native species" ~ "Open-country natives",
              TRUE ~ as.character(ecological_group)
            ),
            ecological_group_label = paste0(
              ecological_group_short,
              "\n(n = ",
              n_group_species,
              ")"
            )
          ) %>%
          pull(ecological_group_label)
      )
    )

  # -----------------------------------------------------------------------------
  # 9. Plotting layers
  # -----------------------------------------------------------------------------

  all_bar_df <- direction_plot_df %>%
    mutate(
      bar_value = signed_all,
      evidence = "All coefficients"
    )

  supported_bar_df <- direction_plot_df %>%
    mutate(
      bar_value = signed_supported,
      evidence = "Supported coefficients"
    )

  # -----------------------------------------------------------------------------
  # 10. Plot
  # -----------------------------------------------------------------------------

  p_habitat_direction_by_covariate <- ggplot() +
    
    geom_vline(
      xintercept = 0,
      linewidth = 0.35,
      colour = "grey35"
    ) +
    
    geom_col(
      data = all_bar_df,
      aes(
        x = bar_value,
        y = fct_rev(ecological_group_label),
        fill = coefficient_direction,
        alpha = evidence
      ),
      width = 0.68,
      colour = NA
    ) +
    
    geom_col(
      data = supported_bar_df,
      aes(
        x = bar_value,
        y = fct_rev(ecological_group_label),
        fill = coefficient_direction,
        alpha = evidence
      ),
      width = 0.48,
      colour = "grey25",
      linewidth = 0.12
    ) +
    facet_wrap(
      facets = vars(predictor_label),
      ncol = 1
    ) +
    
    scale_x_continuous(
      limits = c(-1.12, 1.12),
      breaks = seq(-1, 1, by = 0.5),
      labels = function(x) percent(abs(x), accuracy = 1),
      expand = c(0, 0)
    ) +
    
    scale_fill_manual(
      values = c(
        Negative = "#C95A38",
        Positive = "#007C7A"
      ),
      name = "Coefficient sign"
    ) +
    
    scale_alpha_manual(
      values = c(
        "All coefficients" = 0.23,
        "Supported coefficients" = 1
      ),
      name = "Bar shade"
    ) +
    
    guides(
      fill = guide_legend(
        order = 1,
        ncol = 1, # Stacks legends vertically
        override.aes = list(alpha = 1),
        title.position = "top"
      ),
      alpha = guide_legend(
        order = 2,
        ncol = 1, # Stacks legends vertically
        override.aes = list(
          fill = "grey35",
          colour = NA
        ),
        title.position = "top"
      )
    ) +
    
    labs(
      x = "Percentage of species by coefficient direction",
      y = NULL
    ) +
    
    theme_minimal(base_size = 10.5) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_blank(),
      
      panel.border = element_rect(
        colour = "grey82",
        linewidth = 0.30,
        fill = NA
      ),
      
      strip.background = element_rect(
        fill = "grey95", 
        colour = "grey82", 
        linewidth = 0.1
      ),
      strip.text.x = element_text(
        size = 9.5,
        colour = "grey10",
        face = "bold",
        margin = margin(t = 6, b = 6) 
      ),
      
      axis.text.y = element_text(
        colour = "grey20",
        size = 8.6,
        lineheight = 0.95
      ),
      axis.text.x = element_text(
        colour = "grey25",
        size = 8.4
      ),
      axis.title.x = element_text(
        colour = "grey10",
        size = 9.5,
        margin = margin(t = 8)
      ),
      axis.ticks = element_line(
        colour = "grey55",
        linewidth = 0.25
      ),
      
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.spacing.x = unit(5, "mm"),
      legend.title = element_text(size = 8.8),
      legend.text = element_text(size = 8.5),
      legend.key.height = unit(4, "mm"),
      legend.key.width = unit(8, "mm"),
      
      panel.spacing.y = unit(4, "mm"), 
      plot.margin = margin(10, 10, 10, 10)
    )

  ggsave(
    filename = file.path(figure_dir, "gllvm_habitat_direction.svg"),
    plot = p_habitat_direction_by_covariate,
    width = 10,
    height = 8.2,
    dpi = 800,
    units = "in"
  )

  message("Finished GLLVM coefficient and trait figures")


  message("Completed GLLVM publication outputs")

  output_files <- publication_output_paths(cfg)$gllvm
  missing_files <- output_files[!file.exists(output_files)]
  if (length(missing_files) > 0L) {
    stop("Publication files were not created: ",
         paste(missing_files, collapse = ", "), call. = FALSE)
  }
  output_files
}
