# ReadMe Functions

## Overview

This document summarises the functions developed for this project and their usages.

### Reading in raw data 

[species_data.R](species_data.R)

Script that has **1** function: *read_speciesdata()*

If you need to clean imported long-term species survey data use *read_speciesdata()*

-	Returns a dataframe of combined survey data across the years with corrected species names using the WoRMS database


[grab_data.R](grab_data.R)

Script that has **2** functions *read_grabdata()* and *validate_records_and_grabs()*

*read_grabdata()*

If you need to read in function for grabsite data- environmental (PSA and sediment classification) and spatial (depth, latitude and longitude)

*validate_records_and_grabs()*
If you need to validate the grab site and date information in the grab dataframe and the species survey dataframe, use this function to look for errors


### Data Manipulation

[presence_absence.R](presence_absence.R)

Script that has **1** function: *function_presence_and_absence()*

If you need to create a dataframe of only presence and absence counts. Raw data changes any presence value (P or p) and any value greater than 1 to all equal 1. Counts of 0 remain 0. 


[spatial_protection_levels.R](spatial_protection_levels.R)
 
 This script has **2** functions: *add_protection_zones()* and *add_zones_protect_depth_sedi()*
 
 *add_protection_zones()* 
 If you need to assign spatial protection zones (by fishing gear restrictions) to grab sites based on an MPA management-zone shapefile, and then attach the results to grab dataframe by site
 
 *add_zones_protect_depth_sedi()*
 If you need to convert depth and sediment grain size into categorical variables 
 
 
### Biological Trait Analysis

[category_bins_fuzzycode.R](category_bins_fuzzycode.R)

Script that has **4** functions: *get_breaks()*, *parse_range()*, *code_numeric_range()*, and *code_keyword_match()*

These are different functions that relate to the BIOTIC databse in order to extract categories and create column bins for them 


*get_breaks()*

If you need to create columns of Biotic data with numeric range in their name and extracts lower and upper number boundary for size of species


*parse_range()*

If you need to parse free-text ranges from Biotic into numeric val_min/val_max and NA if nothing is reported


*code_numeric_range()*

If you need to take a species' (val_min, val_max) range and a trait group's bin boundaries, then across every bin the range overlaps and splits target_sum (3 for fuzzy coding) evenly between them
 
 
*code_keyword_match()*

If you need to the categorical traits from Biotic, matches keywords and assigns fuzzy coding to matrix
 
 
[BTA_matrix.R](BTA_matrix.R)

Script that has **2** functions: *build_empty_trait_matrix()* and *check_trait_matrix()*

*build_empty_trait_matrix()*

If you need to build an empty trait matrix for fuzzy coding for BTA


*check_trait_matrix()*

If you need to check the BTA trait table fuzzy coding adds up (either to 3 or 0). Input validating step.


[BTA_PCA.R](BTA_PCA.R)

Script that has **5** functions: *standardise_traits()*, *build_sample_trait_matrix()*, *run_bta_pca()*, *run_pca_lm()* and *run_bta_pca_workflow()*

The first 4 functions are aspects that come together for the full workflow analysis step in *run_bta_pca_workflow()* function

- *standardise_traits()* creates weighted values within columns for fuzzy coding (BTA)

- *build_sample_trait_matrix()* combines abundance data with standardised trait data to produce community-weighted mean traits per sample

- *run_bta_pca()* that runs a PCA on the sample x trait matrix

- *run_pca_lm()* that fits a linear model using the formula "PC1/PC2 ~ management period(or protection level) + phi + depth"


 *run_bta_pca_workflow()* function returns a list: sample-trait matrix, PCA full output, plots and linear model output
 
 
### Data Analysis 

[Analysis_CommunityStructure.R](Analysis_CommunityStructure.R)

Script that has **3** functions: *prepare_community_data()*, *run_nmds_analysis()* and *analyses_bio_patterns()*


*prepare_community_data()*

If you need to combine raw survey data (species and grabs) and create a clean abundance matrix for future analysis that returns a list of matrix, meta, sample_ids, and combined which is metadata + species matrix combined


*run_nmds_analysis()*

If you need to create an NMDS ordination and with support statistical tests such as PERMANOVA, betadisper, SIMPER


*analyses_bio_patterns()*

If you need to calculate biodiversity measurements (Richness, Shannon Diversity, Pielou's evenness) from a community matrix and test them against group covariates (period/protection levels) using linear models


[functional_redundancy.R](functional_redundancy.R)

Script that has **4** functions: *plot_fd_h_scatter()*, *test_slope_difference()*, *run_fr_lm()* and *run_fd_fr_workflow()*


The first 3 functions are aspects that come together for the full workflow analysis step in *run_fd_fr_workflow()* function

- *plot_fd_h_scatter()* helps visualise functional redundancy by plotting FD/H' (scatter plot)

- *test_slope_difference()* tests whether the FD~H' slope differ significantly between groups 

- *run_fr_lm()* linear models with funcitonal redundancy (FR) itself as the response variable, against period/protection + phi/depth (explanatory variables)


[SAC.R](SAC.R)

Script with **1** function: *run_station_sac()*

If you need to plot a species accumulation curve

