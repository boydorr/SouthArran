#Presence and Absence analysis 
  # crude analysis 
# Building a presence and absence abundance version 

# converting P and numbers to 1, keeping 0 at 0 
function_presence_and_absence <- function(survey_data,
                                 abund_col = "value"){
  # Building a presence and absence abundance version of survey data
  #' @description This function returns any presence value (P or p) and any value greater than 1 to all equal 1. Counts of 0 remain 0. This is to create a crude dataframe of only presence and absence counts. 
  #' @param survey_data Dataframe of the species survey records
  #' @param abund_col The name of the column that has the values (abundance count)
  #'
  #'
  #' @return Dataframe of survey data that are just presence and absence counts of 1 and 0 respectively
  
  survey_data_pa <- survey_data %>% # start pipeline
    mutate(!!abund_col := case_when( # over writing abundance column - value
      # conditions
      .data[[abund_col]] %in% c("P", "p")        ~ 1,   # presence-only records assigned to a value of 1
      suppressWarnings(as.numeric(.data[[abund_col]])) == 0 ~ 0, # 0 is set to 0 
      suppressWarnings(as.numeric(.data[[abund_col]])) > 0  ~ 1, # anything above 1 is set to 1
      TRUE ~ NA_real_))
  
  return(survey_data_pa)
}