---
output:
  word_document: default
  html_document: default
---
# South Arran MSc Project README File
## Sarah Stevenson
#### Date:2026-06-11. Last Modified [19/08/2026]


### Project Description 

#### Overview of project
This project has multiple components that help evaluate ecological changes before and after management implementation within the South Arran Marine Protected Area (MPA).
The project combines all long-term benthic survey datasets collected within the MPA. These datasets are standardised and cleaned, with a particular focus on taxonomic records and survey metadata, to create a consistent dataset suitable for biological analyses.
The project also includes a function that helps build an empty Biological Trait Analysis (BTA) fuzzy coding matrix. It first retrieves and fills in species information from BIOTIC for the benthic species present in the dataset and then focuses only on the functional traits relevant to this study: bioturbation, body size, longevity, feeding mode, morphology, and mobility.

**Main features:**

- Community Structure and taxonomic diversity analyses 

- Construction of a functional trait matrix for species recorded in survey datasets

- Biological Trait Analysis (BTA)

- Evaluation of survey effort sufficiency and sampling completeness

- Map of South Arran MPA, fishing gear restrictions, survey stations and grabs samples visualised

### Usage 
Instructions for reproducing the data cleaning, processing and analysis workflow will be included

### Contents

#### Scripts and Functions

##### Scripts 
***project_script.Rmd***
- This is the completed workflow of the project, including instructions on how to use the functions I made. 
The main aspects of this script involve reading and cleaning the raw data files, creating the empty BTA fuzzy coding matrix, and conducting the analyses used to answer my thesis objectives. Initially, the script is applied to all data points before being filtered to include only T survey stations, which are more relevant to my study and contain complete data. 
Biodiversity metrics and NMDS plots were then generated to analyse shifts in community structure and address Objective 1. 
Next, the manually completed BTA matrix was imported and analysed to investigate changes in functional traits, addressing Objective 2. The outputs included a PCA, a biplot, and corresponding statistical analyses using linear models. Finally, a species accumulation curve based on Shannon diversity and grab samples was produced to assess current survey effort.

***map.R***
- Script that produces a map of the survey sites within South Arran MPA
- Fishing gear restrictions included


***LRT-model.R***
Script that runs Likelihood Ratio Tests (LRTs). Using most plausible complex model was plotted, and LRTs were used to assess the significance of interaction and main effects


##### Functions
***function_species_data.R***

Script that has **1** function: *read_speciesdata()*
If you need to clean imported long-term species survey data use *read_speciesdata()*
-	Returns a dataframe of combined survey data across the years with corrected species names using the WoRMS database

***function_grab_data.R***
Script that has **2** functions *read_grabdata()* and *validate_records_and_grabs()*

*read_grabdata()*
If you need to read in function for grabsite data- environmental (PSA and sediment classification) and spatial (depth, latitude and longitude)

*validate_records_and_grabs()*
If you need to validate the grab site and date information in the grab dataframe and the species survey dataframe, use this function to look for errors

***presence_absence.R***
Script that has **1** function: *function_presence_and_absence()*
If you need to create a dataframe of only presence and absence counts. Raw data changes any presence value (P or p) and any value greater than 1 to all equal 1. Counts of 0 remain 0. 

***category_bins_fuzzycode.R***
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
 
 ***spatial_protection_levels.R***
 This script has **2** functions: *add_protection_zones()* and *add_zones_protect_depth_sedi()*
 
 *add_protection_zones()* 
 If you need to assign spatial protection zones (by fishing gear restrictions) to grab sites based on an MPA management-zone shapefile, and then attach the results to grab dataframe by site
 
 *add_zones_protect_depth_sedi()*
 If you need to convert depth and sediment grain size into categorical variables 
 
***function_BTA_matrix.R***
Script that has **2** functions: *build_empty_trait_matrix()* and *check_trait_matrix()*

*build_empty_trait_matrix()*
If you need to build an empty trait matrix for fuzzy coding for BTA

*check_trait_matrix()*
If you need to check the BTA trait table fuzzy coding adds up (either to 3 or 0). Input validating step.


***Analysis_CommunityStructure.R***
Script that has **3** functions: *prepare_community_data()*, *run_nmds_analysis()* and *analyses_bio_patterns()*

*prepare_community_data()*
If you need to combine raw survey data (species and grabs) and create a clean abundance matrix for future analysis that returns a list of matrix, meta, sample_ids, and combined which is metadata + species matrix combined

*run_nmds_analysis()*
If you need to create an NMDS ordination and with support statistical tests such as PERMANOVA, betadisper, SIMPER


*analyses_bio_patterns()*
If you need to calculate biodiversity measurements (Richness, Shannon Diversity, Pielou's evenness) from a community matrix and test them against group covariates (period/protection levels) using linear models


***function_BTA_PCA.R***
Script that has **5** functions: *standardise_traits()*, *build_sample_trait_matrix()*, *run_bta_pca()*, *run_pca_lm()* and *run_bta_pca_workflow()*

The first 4 functions are aspects that come together for the full workflow analysis step in *run_bta_pca_workflow()* function

- *standardise_traits()* creates weighted values within columns for fuzzy coding (BTA)
- *build_sample_trait_matrix()* combines abundance data with standardised trait data to produce community-weighted mean traits per sample
- *run_bta_pca()* that runs a PCA on the sample x trait matrix
- *run_pca_lm()* that fits a linear model using the formula "PC1/PC2 ~ management period(or protection level) + phi + depth"

 *run_bta_pca_workflow()* function returns a list: sample-trait matrix, PCA full output, plots and linear model output
 

***function_functional_redundancy.R***
Script that has **4** functions: *plot_fd_h_scatter()*, *test_slope_difference()*, *run_fr_lm()* and *run_fd_fr_workflow()*

The first 3 functions are aspects that come together for the full workflow analysis step in *run_fd_fr_workflow()* function
- *plot_fd_h_scatter()* helps visualise functional redundancy by plotting FD/H' (scatter plot)
- *test_slope_difference()* tests whether the FD~H' slope differ significantly between groups 
- *run_fr_lm()* linear models with funcitonal redundancy (FR) itself as the response variable, against period/protection + phi/depth (explanatory variables)

***SAC.R***
Script with **1** function: *run_station_sac()*
If you need to plot a species accumulation curve


#### Data 
##### Survey and Grab Data
Survey Data was extracted from NatureScot reports, corresponds with excel files in survey_data and grab_data folders.
- Species abundance count and the grab site it was collected from

2012 
*"species_data_2012.xlsx"* 
https://www.nature.scot/sites/default/files/2025-06/naturescot-commissioned-report-539.pdf

2013 
*"species_data_2013.xlsx"* 
https://www.nature.scot/sites/default/files/2025-06/naturescot-commissioned-report-745.pdf

2015
- July 
*"species_data_2015_1.xlsx"*
https://www.nature.scot/doc/naturescot-commissioned-report-945-infaunal-and-psa-analyses-benthic-samples-collected-around-isle

- September 
*"species_data_2015_2.xlsx"*
https://www.nature.scot/doc/naturescot-commissioned-report-946-infaunal-and-psa-analyses-benthic-samples-collected-south-arran

2022
*"species_data_2022.xlsx"* 
Survey data provided by Dr. Kelly Saunders, NatureScot

##### Map details
Folder - *Fishing measures relevant to South Arran*
Shapefiles of fishing management boundaries

Folder - *Survey boxes - 2021 and 2016*
Shapefiles of the survey stations introduced from 2015 onwards

Provided by Dr. Kelly Saunders, NatureScot

##### Taxonomic Information 
Taxonomic names and classifications were validated using the World Register of Marine Species (WoRMS)
https://www.marinespecies.org/aphia.php


##### Trait Information

*"biotic.csv"* - An csv file that returns species taxonomic information and trait information. Accessed [16/06/2026]: https://api.mba.ac.uk/help_marlin 


*Websites for manual assignment* 
Databases to search species names for their traits to help with assignment of fuzzy coding.

Encyclopedia of Life website
https://eol.org/ 

Sea Life Base
https://www.sealifebase.ca/search.php 

World Register of Marine Species (WoRMS)
https://www.marinespecies.org/aphia.php

##### Additional 
**Folders** 

*data* - Folder containing the survey data- the species abundance count (/survey_data), environmental information (/grab_data), fishing regulation boundaries (/Fishing measures relevant to South Arran), biotic.csv, and name corrections to be read through before submitting names through the WoRMS database (unassigned species list.xlsx).

*outputs* - Folder containing the four outputs returned from the BTA function that creates a empty matrix for trait information and matches species names with the biotic dataframe. 

### Credits 

**Supervisors:**
Richard Reeve, University of Glasgow

Kelly Saunders, NatureScot

### Contact Information 
Sarah Stevenson 

3170686S@student.gla.ac.uk 

MSc Project, University of Glasgow 
