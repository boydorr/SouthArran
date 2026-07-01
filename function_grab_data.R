# Function for reading in grabdata

read_grabdata <- function(path = "data/grab_data/grabsites.xlsx",
                          records = NULL){
  
  #'@description This function reads in all the grab environmental and spatial data, accounts for different date formats, and validates the grab information against the records from survey data.
  #' 
  #' @param path Setting the file location to read grab data 
  #' @param records Arguement to choose what records to pass through to validate data
  #' 
  #' @returns Cleaned grabsite survey dataframe
  
  #1. Read data
  grabs <- readxl::read_excel(
    path,
    col_types = c("guess", "text", rep("guess", 8)) # need to set a text for date column cause 2022 data only has year
  ) %>%
    select(-matches("\\.\\.\\.\\d+$")) %>% # drop stray empty columns like previous
    
  #2. Parse mixed date formats
    mutate(
      Date_raw = Date,
      Date = case_when(
        str_detect(Date, "^\\d{4}$") ~ NA_Date_,
        str_detect(Date, "/|-")      ~ dmy(Date, quiet = TRUE),
        TRUE ~ as.Date(as.integer(Date), origin = "1899-12-30")
      ),
      Year = case_when(
        str_detect(Date_raw, "^\\d{4}$") ~ as.integer(Date_raw),
        TRUE                             ~ year(Date)
      )
    ) %>%
    select(-Date_raw)%>% 
    
    # Split GrabSite into base site and grab replicate number
    mutate(
      GrabSite_base = str_replace(GrabSite, "-\\d+$", ""),
      GrabNumber    = case_when(
        str_detect(GrabSite, "-\\d+$") ~ str_extract(GrabSite, "\\d+$"),
        TRUE ~ "1"
      )
    )
  
  #3. Validate against records
  if (!is.null(records))
  {
    success = validate_records_and_grabs(records, grabs)
    if (success != TRUE)
    {
      cat(success)
      stop("Problem with grab sites")
    }
  }
  
 return(grabs)
}


# Function for validating 
validate_records_and_grabs <- function(records, grabs) {
  
  #' @description This function validates two dataframes, the species survey data (records) and the environmental and spatial information from grabsites (grabs)
  #' @param records Species survey data
  #' @param grabs The environmental and spatial information from grabsites
  #' 
  #' @return Results if there no mismatches are found in either direction it returns TRUE, meaning the two dataframes are in sync or FALSE if the two dataframes are not a match it stops and builds and error message of where to issue occurs
  
  # Key for each dataframe
  records_keys <- records %>%
    mutate(key = paste(GrabSite_base, .data$year, sep = "_")) %>%
    pull(key)
  
  grabs_keys <- grabs %>%
    mutate(key = paste(GrabSite_base, .data$Year, sep = "_")) %>%
    pull(key)
  
  # Check for grabs missing from records
  missing_from_records <- setdiff(grabs_keys, records_keys)
  
  # Check for records missing from grabs
  missing_from_grabs <- setdiff(records_keys, grabs_keys)
  
  if (length(missing_from_records) == 0 && length(missing_from_grabs) == 0) {
    return(TRUE)
  }
  
  # Error message
  msg <- ""
  if (length(missing_from_records) > 0)
    msg <- paste0(msg, "Grabs missing from records: ", paste(missing_from_records, collapse = ", "), "\n")
  if (length(missing_from_grabs) > 0)
    msg <- paste0(msg, "Records missing from grabs: ", paste(missing_from_grabs, collapse = ", "), "\n")
  
  return(msg)
}

