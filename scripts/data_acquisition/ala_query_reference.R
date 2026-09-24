# =============================================================================
# Script: scripts/data_acquisition/ala_query_reference.R
# Purpose: Preserve the supplied ALA acquisition helper for provenance review.
# Workflow stage: Data acquisition reference; outside the targets pipeline.
# Inputs: A species table with Species.Name and the original query extent.
# Outputs: None when sourced; the function is defined but never called here.
# Execution note: Reproduce the analysis using the frozen ala_df_raw.csv.
# =============================================================================
#
# The supplied download script called this helper with first_year = 1990 and
# batch_size = 10. It split its species table into six groups, combined the
# returned records, and contained additional interactive editing commands.
#
# USER CHECK: the exact taxon list, query boundary, acquisition date and any
# manual edits for the released snapshot are not fully established by those
# scripts. The helper below is preserved as supplied, not claimed to reconstruct
# the exact released CSV. In particular, its selected fields do not explicitly
# include basisOfRecord or coordinateUncertaintyInMeters, which the CSV contains.
#
# Calling this function requires galah to be configured by the caller. Package
# calls are namespace-qualified, so packages do not need to be attached. No
# credentials, network requests or package loading are performed by sourcing
# this reference. _targets.R does not source scripts/data_acquisition/ or submit
# a live query when reproducing the frozen snapshot.
#
# Supplied acquisition helper ----

get_ALA_data <- function(first_year = NULL, dTol = 10000, species_list, shp, output = NULL, batch_size = 50) {
  # Simplify shapefile and create convex hull
  shp_union <- sf::st_union(shp)
  shp_simplified <- sf::st_simplify(shp_union, dTolerance = dTol)
  hull <- sf::st_convex_hull(shp_union)
  
  sp_select <- species_list$Species.Name
  
  # Split species list into batches
  sp_batches <- split(sp_select, ceiling(seq_along(sp_select) / batch_size))
  
  # Initialize empty list to store results
  sp_occ_list <- list()
  
  # Loop over batches
  for (i in seq_along(sp_batches)) {
    sp_batch <- sp_batches[[i]]
    message("Processing batch ", i, " of ", length(sp_batches))
    
    # Build the galah call
    galah_call_obj <- galah::galah_call() |>
      galah::galah_identify(sp_batch) |>
      galah::galah_polygon(hull) |>
      galah::galah_select(
        decimalLongitude, decimalLatitude, family, genus, species,
        vernacularName, eventDate, dataResourceName, samplingProtocol, samplingEffort,
        occurrenceStatus, datasetName
      ) |>
      galah::galah_filter(year >= first_year)
    
    # Retrieve occurrences
    sp_occ <- galah::atlas_occurrences(galah_call_obj)
    
    # Append to list
    sp_occ_list[[i]] <- sp_occ
  }
  
  # Combine all occurrences
  sp_occ_all <- dplyr::bind_rows(sp_occ_list)
  
  # Proceed as before with sp_occ_all instead of sp_occ
  # Rename columns and filter out rows with missing coordinates
  sp_df <- sp_occ_all |>
    dplyr::rename(
      longitude = decimalLongitude,
      latitude = decimalLatitude,
      common_name = vernacularName,
      date = eventDate
    ) |>
    dplyr::filter(!is.na(longitude) & !is.na(latitude))
  
  # Convert to sf object and spatially filter
  sp_sf <- sf::st_as_sf(
    sp_df,
    coords = c("longitude", "latitude"),
    crs = sf::st_crs(shp)
  )
  
  # Spatial filtering using the simplified shapefile
  sp_sf <- sp_sf[shp_simplified, ]
  
  coords <- sf::st_coordinates(sp_sf)
  
  # Extract coordinates and drop geometry
  sp_sf <- sp_sf |>
    dplyr::mutate(
      longitude = coords[, 1],
      latitude = coords[, 2]
    ) |>
    sf::st_drop_geometry()
  
  # Write to CSV if output path is provided
  if (!is.null(output)) {
    write.csv(sp_sf, output, row.names = FALSE)
  }
  
  return(sp_sf)
}
