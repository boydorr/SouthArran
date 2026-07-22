
# Creating Read function
read_speciesdata <- function(
               # Arguments
                data_dir = "data/survey_data",   
                pattern = "species_data",   
                species_col = "Taxa",  
                resolve_worms = TRUE # using WoRMS database to fix taxonomy
                      ){
 #' Read and cleans survey data of species that were found at what grab site
 #' 
 #' @description This function reads in all the survey data, pivots into long format grabs, combines the surveys across the years, then cleans the species names to then be searched through the WoRMS database for accepted scienfic name and larger taxonomic resolution. 
 #' @param data_dir The folder containing with Excel files of survey report
 #' @param pattern The filename pattern to search for and match
 #' @param species_col The name of the species column in excel sheets
 #' @param resolve_worms Arguement so function will automatically validate and fix species names using the WoRMS database
 #' 
 #' @return A dataframe of combined survey data across the years with corrected species names
  
  #1. Read data
    # Retrieving files
    paths <- list.files(data_dir, pattern = pattern, full.names = TRUE)
  
  # Stop immediately if no files are found
  if (length(paths) == 0) stop("ERROR: No files matching '", pattern, "' in: ", data_dir) # returns error message
  message("Reading ", length(paths), " file(s)...")
  
  # Pull the year out of each filename
  years <- gsub(paste0(".*", pattern, "_(\\d{4}).*\\.xlsx$"), "\\1", basename(paths))
  
  #2. Prep for analysis- Process data
      # function in a function
  process_one_file <- function(path, year) {
    readxl::read_excel(path) %>%
      # Removes empty cells after grab columns
      select(-matches("^\\.\\.\\.\\d+$")) %>% 
      
      # Now pivoting data (not species) into long format
      mutate(across(-all_of(species_col), as.character)) %>%
      pivot_longer(
        cols      = -all_of(species_col),
        names_to  = "GrabSite",
        values_to = "value"
      ) %>%
      
      # Replace any NA values with "0"
      mutate(value = replace_na(value, "0")) %>%
      # Add a year column
      mutate(year = year) %>%
      
      # Species to SpeciesName and keep only the columns we need
      select(SpeciesName = all_of(species_col), GrabSite, value, year) %>%
      
      # Split GrabSite into base site and grab number (replicates eg bucket 1)
      mutate(
        GrabSite_base = str_replace(GrabSite, "-\\d+$", ""),
        GrabNumber    = case_when(
          str_detect(GrabSite, "-\\d+$") ~ str_extract(GrabSite, "\\d+$"),
          TRUE ~ "1"
        ),
        
        # Extract station code from GrabSite_base
        GrabSite_station = case_when(
          str_detect(GrabSite_base, "^[DT]\\d+G\\d+$") ~ str_extract(GrabSite_base, "^[DT]\\d+"),
          str_detect(GrabSite_base, "^[DT]\\d+_") ~ str_extract(GrabSite_base, "^[DT]\\d+"),
          str_detect(GrabSite_base, "^G\\d+$") ~ NA_character_,
          TRUE ~ NA_character_
        )
      )
  }
  
  # Combine into one data frame
      combined <- map2(paths, years, process_one_file) %>% # map2 runs through one path and year at a time
      list_rbind()  # list_rbind stacks all the resulting dataframes into one
  
  message("Combined: ", nrow(combined), " rows across years: ",
          paste(sort(unique(combined$year)), collapse = ", "))
  
  #3. Clean data 
      combined <- combined %>% # within new combined dataframe
      
        mutate(SpeciesName_clean = SpeciesName %>%
                 
                 str_trim() %>% # removes leading/trailing whitespace
                 
               # Replacing (juv.), (juv) and (Juv) with (juvenile)
                 str_replace_all("\\(juv\\.?\\)", "(juvenile)") %>%   # (juv) or (juv.)
                 str_replace_all("\\(Juv\\.?\\)", "(juvenile)") %>%   # (Juv) or (Juv.)
                 str_replace_all("\\bjuv\\.\\B",  "(juvenile)") %>%   # juv. 
                 str_replace_all("\\bJuv\\.\\B",  "(juvenile)") %>%   # Juv.
                 str_replace_all("\\bjuv\\b",     "(juvenile)") %>%   # juv 
                 str_replace_all("\\bJuv\\b",     "(juvenile)") %>%   # Juv
                 
              # Convert uppercase words to title case 
                 str_replace_all("\\b[A-Z]{2,}\\b", function(m) str_to_title(m)) %>% # keeping capital from second character ({2})
                 
               # Standardise "spp." variants - to become "sp."
                   str_replace_all("\\bspp\\.{0,2}", "sp.") %>%   # spp.
                   str_replace_all("\\bsp\\.{2,}", "sp.") %>%     # sp..
                   str_replace_all("\\bsp\\b(?!\\.)", "sp.") %>%  # sp
                 
              # Remove aggregate labels like "Lumbrineris agg.", keeping just the genus
                 str_replace_all("\\s+agg\\..*$", "") %>%
                 
               # For slash aggregates
                 #Harmothoe extenuata/fragilis --> make genus Harmothoe
                 str_replace_all("Harmothoe extenuata/fragilis", "Harmothoe") %>%
                 # Lumbrineris cingulata/gracilis --> Lumbrineridae (Family)
                 str_replace_all("Lumbrineris cingulata/gracilis", "Lumbrineridae") %>%
                #Goniadidae/Hesionidae --> Phyllodocida (Order) 
                 str_replace_all("Goniadidae/Hesionidae", "Phyllodocida") %>%
                 
                 str_trim() # Redo: removes leading/trailing whitespace again after substitutions
        ) %>%

      # Cleaning for WoRMS query - further stripped version of the clean name
      # Removes life-stage/sex qualifiers
      mutate(SpeciesName_worms = SpeciesName_clean %>%
               str_replace_all("\\s*\\(juvenile\\)\\s*", " ") %>%  # remove (juvenile)
               str_replace_all("\\s*\\(female\\)\\s*",   " ") %>%  # remove (female)
               str_replace_all("\\s*\\(male\\)\\s*",     " ") %>%  # remove (male)
               str_replace_all("\\s*\\(zoea larvae\\)\\s*",     " ") %>%  # remove (zoea larvae)
               str_replace_all("\\s*\\(larvae\\)\\s*",     " ") %>%  # remove (larvae)
               str_replace_all("\\s*\\(praniza\\)\\s*",     " ") %>%  # remove (praniza)
               str_replace_all("\\s*\\bsp\\.\\s*", " ") %>%  # remove sp.
               
               
               str_squish()  # collapses any leftover internal whitespace
      )
      
      # Common naming corrections 
      corrections_path <- file.path(dirname(data_dir), "unassigned species list.xlsx")  # finding file in data folder
      
      if (file.exists(corrections_path)) {
        corrections <- readxl::read_excel(corrections_path) %>%
          rename(
            SpeciesName_worms = 1,   # old/incorrect name in first column
            corrected_name    = 2    # accepted replacement in second column
          ) %>%
          mutate(across(everything(), str_trim)) %>%   # remove white space
          filter(!is.na(SpeciesName_worms) & SpeciesName_worms != "",
                 !is.na(corrected_name)    & corrected_name    != "")
        
        n_before <- sum(combined$SpeciesName_worms %in% corrections$SpeciesName_worms) # to verify how many names where correct
        combined <- combined %>%
          left_join(corrections, by = "SpeciesName_worms") %>%
          mutate(SpeciesName_worms = coalesce(corrected_name, SpeciesName_worms)) %>% # replacing old name with new corrected names where there is a match in clean name list
          select(-corrected_name)
        
        message("Applied ", n_before, " manual name correction(s) from: ", basename(corrections_path))
        
      } else {
        message("No name corrections file found at: ", corrections_path, " — skipping.")
      }
      
  #4. Resolving taxonomic records using WoRMS database
      if (resolve_worms) { # Only runs if resolve_worms is true
        
        # Creating unique cleaned names from data
        unique_names <- unique(combined$SpeciesName_worms)
        unique_names <- unique_names[
          !is.na(unique_names) & # removing any NA values
            unique_names != "" & # removing empty strings
            !grepl("^\\s*$", unique_names) # removes strings of whitespace 
            ]
        
        message("Querying WoRMS for ", length(unique_names), " unique names...") # Needs to pass above checks
        
        # Use SHARK4R to query on WoRMS database
          # Match_worms_taxa handles retries and strips special characters automatically
        worms_result <- tryCatch(
          SHARK4R::match_worms_taxa(unique_names,
                                    fuzzy=TRUE,
                                    best_match_only = TRUE,
                                    marine_only = TRUE,
                                    verbose = FALSE),
          error = function(e) {
            warning("Failed: ", e$message, " — skipping taxonomy resolution.")
            return(NULL)
          }
        )
        if (!is.null(worms_result)) { # if found message
          
          message(
            "WoRMS returned columns:\n",
            paste(names(worms_result), collapse = ", ")
          )
          
          # Pull columns that are wanted and rename them clearly
          worms_taxonomy <- worms_result %>%
            select(
              SpeciesName_worms = scientificname,  # the name that was queried
              accepted_name     = valid_name,      # WoRMS-accepted current name
              aphia_id          = AphiaID,         # ID
              worms_status      = status,          # "accepted", "synonym", "unaccepted" on WoRMs database
              genus, family, order, class, phylum  # taxonomic information
            ) %>%
            # Keep only the first match if the same name appears twice in WoRMS results 
            distinct(SpeciesName_worms, .keep_all = TRUE)
          
          # Combine the WoRMS search results onto data by the cleaned species name
          combined <- combined %>%
            left_join(worms_taxonomy, by = "SpeciesName_worms")
          
          #resend through Worms to make sure taxonomy is updated
          accepted_unique <- combined %>%
            dplyr::filter(
              !is.na(accepted_name)
            ) %>%
            dplyr::distinct(
              accepted_name
            ) %>%
            dplyr::pull()
          
          message(
            "Re-querying WoRMS for ",
            length(accepted_unique),
            " accepted names..."
          )
          
          accepted_tax <- tryCatch(
            
            SHARK4R::match_worms_taxa(
              accepted_unique,
              fuzzy = FALSE,
              best_match_only = TRUE,
              marine_only = TRUE,
              verbose = FALSE
            ),
            
            error = function(e) {
              warning(
                "Second WoRMS query failed: ",
                e$message
              )
              NULL
            }
          )
          
          if (!is.null(accepted_tax)) {
            
            accepted_taxonomy <- accepted_tax %>%
              dplyr::select(
                accepted_name = valid_name,
                genus_new = genus,
                family_new = family,
                order_new = order,
                class_new = class,
                phylum_new = phylum,
                aphia_id_new    = AphiaID,
                worms_status_new = status
                
              ) %>%
              dplyr::distinct(
                accepted_name,
                .keep_all = TRUE
              )
            
            combined <- combined %>%
              dplyr::left_join(
                accepted_taxonomy,
                by = "accepted_name"
              ) %>%
              dplyr::mutate(
                genus = dplyr::coalesce(
                  genus_new,
                  genus
                ),
                family = dplyr::coalesce(
                  family_new,
                  family
                ),
                order = dplyr::coalesce(
                  order_new,
                  order
                ),
                class = dplyr::coalesce(
                  class_new,
                  class
                ),
                phylum = dplyr::coalesce(
                  phylum_new,
                  phylum
                ),
               aphia_id = dplyr::coalesce(aphia_id_new, 
                                          aphia_id),
            worms_status = dplyr::coalesce(worms_status_new,
                                           worms_status)) %>%
              dplyr::select(
                -dplyr::ends_with("_new")
              )
}
          unresolved <- combined %>%
            dplyr::filter(
              is.na(accepted_name)
            ) %>%
            dplyr::distinct(
              SpeciesName_worms
            ) %>%
            dplyr::pull(
              SpeciesName_worms
            )
          
          if (length(unresolved) > 0) {
            
            message(
              length(unresolved),
              " unresolved name(s):\n",
              paste(
                " -",
                unresolved,
                collapse = "\n"
              )
            )
          }
        }
      }
    return(combined)
} 
