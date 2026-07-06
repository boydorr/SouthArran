# BTA coding 

# remove coding source, date etc

build_empty_trait_matrix <- function(
    survey_data, 
    biotic_path,
    output_dir        = "outputs",
    site_filter       = c("standard", "fetlar"),
    abund_col         = "value",
    min_site_prop     = 0.50,
    top_abund_prop    = 0.90,
    fetlar_abund_prop = 0.75
) {
  # need to fill in!
  #' What does function do in a line?
  #' 
  #' @description 
  #' @param 
  #' @param 
  #' 
  #' @return 
  
  site_filter <- match.arg(site_filter)
  
  # Trait Modalities - selected relevant to study
  trait_groups <- list(
    bioturbation = c("bio_none", "bio_biodiff", "bio_upward",
                     "bio_downward", "bio_bioirr", "bio_surf_mod"),
    size         = c("size<5", "size_5_10", "size_10_20",
                     "size_20_40", "size_40_80", "size>80"),
    longevity    = c("life<1", "life_1_3", "life_3_5",
                     "life_5_10", "life>10"),
    feeding      = c("feed_deposit", "feed_suspension", "feed_detritus",
                     "feed_scavenger", "feed_grazer", "feed_predator"),
    morphology   = c("morph_soft", "morph_tunic", "morph_exoskeleton",
                     "morph_crustose", "morph_cushion", "morph_stalked", "morph_endoskeleton"),
    mobility     = c("mob_fixed", "mob_limited", "mob_slow", "mob_mobile")
  )
  
  all_trait_cols <- unlist(trait_groups, use.names = FALSE)
  
  # Resolve bta_name 
  message("Step 1: Resolving BTA species names...")
  
  survey_data_bta <- survey_data_adjusted %>%
    mutate(bta_name = case_when(
      !is.na(accepted_name) & accepted_name != "" ~ accepted_name,
      TRUE ~ SpeciesName_clean
    ))
  
  # Apply inclusion criteria 
  message("Step 2: Applying '", site_filter, "' inclusion criteria...")
  
  species_summary <-  survey_data_bta%>%
    mutate(abundance = suppressWarnings(as.numeric(.data[[abund_col]]))) %>%
    filter(!is.na(abundance), abundance > 0) %>%
    group_by(bta_name) %>%
    summarise(
      n_sites     = n_distinct(GrabSite), #add _base????
      total_abund = sum(abundance, na.rm = TRUE),
      .groups     = "drop"
    )
  
  total_sites <- n_distinct( survey_data_bta$GrabSite) #add _base????
  
  cumulative <- species_summary %>%
    arrange(desc(total_abund)) %>%
    mutate(cum_prop = cumsum(total_abund) / sum(total_abund))
  
  bta_included <- if (site_filter == "standard") {
    cumulative %>%
      filter((n_sites / total_sites) >= min_site_prop | cum_prop <= top_abund_prop)
  } else {
    cumulative %>%
      filter(cum_prop <= fetlar_abund_prop)
  }
  
  message("  ", nrow(bta_included), " / ", nrow(species_summary),
          " species pass inclusion criteria")
  
  # Attach taxonomy 
  taxonomy <- survey_data_bta %>%
    distinct(bta_name, aphia_id, genus, family, order, class, phylum)
  
  bta_included <- bta_included %>%
    left_join(taxonomy, by = "bta_name")
  
  # Flag BIOTIC coverage 
  message("Step 3: Checking BIOTIC coverage...")
  
  biotic <- read_csv(biotic_path, show_col_types = FALSE) %>%
    mutate(biotic_name = str_squish(SpeciesName)) %>%
    select(biotic_name, Size, LifeSpan,        # Biotic columns (relevant traits) are selected
           feedingmethod, mobility, growthform,           
           Bioturbator, Habit, biozone)
  
  bta_included <- bta_included %>%
    left_join(biotic, by = c("bta_name" = "biotic_name")) %>%
    mutate(
      in_biotic = !is.na(feedingmethod),
      # Flag higher-taxon entries for family-level grouping
      taxon_level = case_when(
        str_detect(bta_name, " ")  ~ "species",   # has a space = binomial
        bta_name %in% taxonomy$family ~ "family",
        bta_name %in% taxonomy$order  ~ "order",
        bta_name %in% taxonomy$class  ~ "class",
        TRUE                           ~ "genus"
      )
    )
  
  n_biotic    <- sum(bta_included$in_biotic)
  n_not_biotic <- sum(!bta_included$in_biotic)
  
  message("  In BIOTIC:     ", n_biotic)
  message("  Not in BIOTIC: ", n_not_biotic, " (need manual research)")
  
  # Build empty matrix 
  message("Step 4: Building empty trait matrix...")
  
  empty_matrix <- bta_included %>%
    select(bta_name, aphia_id, taxon_level,
           genus, family, order, class, phylum,
           n_sites, total_abund, in_biotic,
           Size, LifeSpan, feedingmethod, # traits
           mobility, growthform, Bioturbator) %>%
    # Add all trait columns as NA (not 0 because it is not yet coded- needs to be filled in)
    bind_cols(
      matrix(NA_real_,
             nrow = nrow(bta_included),
             ncol = length(all_trait_cols),
             dimnames = list(NULL, all_trait_cols)) %>%
        as_tibble()
    )
  
  # Group summary for family-level coding
  # Shows which families have multiple species — candidates for grouped coding
  family_groups <- bta_included %>%
    filter(!in_biotic) %>%                      # focus on unmatched species
    group_by(phylum, family, order) %>%
    summarise(
      n_species        = n(),
      species          = paste(bta_name, collapse = "; "),
      any_in_biotic    = any(in_biotic),
      .groups          = "drop"
    ) %>%
    arrange(desc(n_species))
  
  # Outputs 
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Output 1. Empty matrix to fill in manually
  write_csv(empty_matrix,
            file.path(output_dir,
                      paste0("trait_matrix_TO_FILL_", site_filter, ".csv")))
  
  # Output 2. BIOTIC-matched species only (with raw BIOTIC fields alongside)
  write_csv(
    empty_matrix %>% filter(in_biotic),
    file.path(output_dir,
              paste0("biotic_matched_", site_filter, ".csv"))
  )
  
  # Output 3. Unmatched species needing manual research - sorted by family for grouping
  write_csv(
    empty_matrix %>% filter(!in_biotic) %>%
      arrange(phylum, family, bta_name),
    file.path(output_dir,
              paste0("manual_research_needed_", site_filter, ".csv"))
  )
  
  # Output 4. Family grouping summary
  write_csv(family_groups,
            file.path(output_dir,
                      paste0("family_groups_", site_filter, ".csv")))
  
  message("\nOutputs written to: ", output_dir)
  message("  trait_matrix_TO_FILL_", site_filter,
          ".csv       — full matrix to complete (", nrow(empty_matrix), " species)")
  message("  biotic_matched_", site_filter,
          ".csv            — ", n_biotic, " species with BIOTIC reference data")
  message("  manual_research_needed_", site_filter,
          ".csv    — ", n_not_biotic, " species needing literature search")
  message("  family_groups_", site_filter,
          ".csv           — family groupings for bulk coding")
  
  invisible(list(
    empty_matrix   = empty_matrix,
    family_groups  = family_groups,
    trait_groups   = trait_groups    # return structure for reference
  ))
}


# Function that checks the BTA trait table 

check_trait_matrix <- function(filled_matrix, trait_groups, target_sum = 3) {
  
  all_trait_cols <- unlist(trait_groups, use.names = FALSE)
  
  # 1. Check no trait columns got read in as character
  col_classes <- filled_matrix %>%
    select(all_of(all_trait_cols)) %>%
    summarise(across(everything(), class)) %>%
    unlist()
  
  non_numeric <- col_classes[col_classes != "numeric"]
  
  if (length(non_numeric) > 0) {
    message("ERROR: ", length(non_numeric), " trait column(s) not numeric:")
    print(non_numeric)
  } else {
    message("WORKS: all trait columns are numeric")
  }
  
  # 2. Check each trait group sums to target_sum (3 for fuzzy coding) or 0 (no information on trait) per species
  row_sum_issues <- purrr::map_dfr(trait_groups, function(cols) {
    filled_matrix %>%
      mutate(row_sum = rowSums(across(all_of(cols)), na.rm = TRUE)) %>%
      filter(!row_sum %in% c(target_sum, 0) & !is.na(row_sum)) %>%
      select(bta_name, row_sum)
  }, .id = "trait_group")
  
  if (nrow(row_sum_issues) > 0) {
    message("WARNING: ", nrow(row_sum_issues),
            " species/trait-group combinations don't sum to ", target_sum, " or 0")
  } else {
    message("GREAT: all trait groups sum to ", target_sum, " or 0 per species")
  }
  
  invisible(list(
    non_numeric_cols = non_numeric,
    row_sum_issues   = row_sum_issues
  ))
}