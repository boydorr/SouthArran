# Libraries 
library(tibble)
library(purrr)
library(stringr)

# Adjusting for fuzzy coding 
  # Continuous traits (size and longevity) in Biotic 
  # Marking boundaries of trait groups 

get_breaks <- function(cols) { 
  
  #' @description Function that takes columns of Biotic data with numeric range in their name and extracts lower and upper number boundary
  #'  
  #'  @param cols The column names that will become the upper and lower boundaries of the bins
  #'  @return Column of extracted upper and lower values 
  
  nums <- str_extract_all(cols, "\\d+\\.?\\d*") # extract all numbers in the column name (e.g. 10-20 --> 10, 20)
  tibble( col   = cols, 
          lower = map2_dbl(cols, nums, ~ if (str_detect(.x, "<")) 0 else as.numeric(.y[1])),
          upper = map2_dbl(cols, nums, ~ if (str_detect(.x, ">")) Inf else as.numeric(.y[length(.y)])))
}

# Parses free-text ranges from Biotic into numeric val_min / val_max and NA if nothing is reported
parse_range <- function(text) {
  #' @description Function that parses free-text ranges from Biotic into numeric val_min/val_max and NA if nothing is reported
  #' @param text Text from columns in Biotic
  #' 
  #' @return Tibble with one row per input string
  
  text_l <- str_to_lower(text)
  nums   <- str_extract_all(text_l, "\\d+\\.?\\d*") %>% map(as.numeric) # converting to a numeric
  
  map2_dfr(text_l, nums, function(t, n) { # a function to return a dataframe - looping over each text string and numeric vector (t and n)
    if (is.na(t) || length(n) == 0) return(tibble(val_min = NA_real_, val_max = NA_real_))
    if (str_detect(t, "<")) return(tibble(val_min = 0,      val_max = n[1]))
    if (str_detect(t, ">")) return(tibble(val_min = n[1],   val_max = Inf))
    tibble(val_min = min(n), val_max = max(n))
  })
  
}

# Given a species' (val_min, val_max) range and a trait group's bin boundaries, finds every bin the range overlaps and splits target_sum evenly across them
  # Bins with no overlap / no data stay NA 
code_numeric_range <- function(val_min, 
                               val_max,
                               breaks,
                               target_sum = 3) {   # target sum here is the fuzzy coding assignment (3)
  #' @description Function that takes a species' (val_min, val_max) range and a trait group's bin boundaries, finds every bin the range overlaps and splits target_sum evenly across them
  #' @param val_min lower range boundary
  #' @param val_max upper range boundary
  #' @param breaks the bin lookup
  #' 
  #' @return A row per species (with a complete column per bin)
  
  map2_dfr(val_min, val_max, function(mn, mx) {
    vals <- set_names(rep(NA_real_, nrow(breaks)), breaks$col)
    if (is.na(mn) || is.na(mx)) return(vals)
    spanned <- breaks$col[breaks$lower < mx & breaks$upper > mn] # finding overlapping bins
    if (length(spanned) > 0) vals[spanned] <- target_sum / length(spanned) # for every overlapping the fuzzy coding value is shared over the number of bins
    vals # returning the value
  })
}

# Categorical Traits (e.g. feeding type) from Biotic
# keyword_map is a named list: 
  # names = target trait columns
  # values = keyword to look for in the source text
# Any number of keywords can match at once (e.g. "deposit and suspension feeder" matches both feed_deposit and feed_suspension), and target_sum is split evenly across whatever matched
code_keyword_match <- function(text,
                               keyword_map,
                               target_sum = 3) {
  #' @description Function the uses the categorical trait from Biotic, matches keywords and assigns fuzzy coding to matrix
  #' @param text Text vectors of categorical traits descriptions
  #' @param target_sum Fuzzy code total for distribution amoung bins
  #' 
  #' @return  A row per species
  
  map_dfr(text, function(t) {
    vals <- set_names(rep(NA_real_, length(keyword_map)), names(keyword_map)) # using keyword map
    if (is.na(t)) return(vals)
    t_l <- str_to_lower(t)
    matched <- names(keyword_map)[map_lgl(keyword_map, ~ str_detect(t_l, .x))]
    if (length(matched) > 0) vals[matched] <- target_sum / length(matched)
    vals
  })
  
}
