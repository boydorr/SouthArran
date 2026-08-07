---
output:
  word_document: default
  html_document: default
---
# South Arran MSc Project README File
## Sarah Stevenson
#### Date:2026-06-11. Last Modified [07/08/2026]


### Project Description 
#### Overview of project

This project combines long-term benthic survey datasets collected within the South Arran Marine Protected Area (MPA). The data are standardised and cleaned, with a particular focus on taxonomic records and survey metadata, to create a consistent dataset suitable for biological analyses.

**Main features:**

- Community Structure and taxonomic diversity analyses 

- Construction of a functional trait matrix for species recorded in survey datasets

- Biological Trait Analysis (BTA)

- Evaluation of survey effort sufficiency and sampling completeness

### Usage 
Instructions for reproducing the data cleaning, processing and analysis workflow will be included

### Contents

#### Scripts and Functions

##### Scripts 
*project_script.Rmd* 
- WHAT?

*map.R*
- Map of the survey sites across the years

##### Functions
***function_species_data.R***

Script that has 1 function *read_speciesdata()*
If you need to clean imported long-term species survey data use *read_speciesdata()*
-	Returns a dataframe of combined survey data across the years with corrected species names using the WoRMS database


***function_grab_data.R***
Script that has 2 functions *read_grabdata()* and *validate_records_and_grabs()*

*read_grabdata()*
If you need to read in function for grabsite data- environmental (PSA and sediment classification) and spatial (depth, latitude and longitude)

*validate_records_and_grabs()*
If you need to validate the grab site and date information in the grab dataframe and the species survey dataframe, use this function to look for errors

***function_BTA.R***
Script that has 2 functions *build_empty_trait_matrix()* and *check_trait_matrix()*

*build_empty_trait_matrix()*
If you need to build...

*check_trait_matrix()*
If you need to check the BTA trait table fuzzy coding adds up (either to 3 or 0). Input validating step.

***PCA.R***


***Analysis_CommunityStructure.R***


***SAC.R***


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
