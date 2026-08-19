# South Arran MSc Project README File

Sarah Stevenson

*Date:2026-06-11. Last Modified [19/08/2026]*


## Project Description 

### Overview of project

This project has multiple components that help evaluate ecological changes before and after management implementation within the South Arran Marine Protected Area (MPA).
The project combines all long-term benthic survey datasets collected within the MPA. These datasets are standardised and cleaned, with a particular focus on taxonomic records and survey metadata, to create a consistent dataset suitable for biological analyses.
The project also includes a function that helps build an empty Biological Trait Analysis (BTA) fuzzy coding matrix. It first retrieves and fills in species information from BIOTIC for the benthic species present in the dataset and then focuses only on the functional traits relevant to this study: bioturbation, body size, longevity, feeding mode, morphology, and mobility.

### Main features:

- Community Structure and taxonomic diversity analyses 

- Construction of a functional trait matrix for species recorded in survey datasets

- Biological Trait Analysis (BTA)

- Evaluation of survey effort sufficiency and sampling completeness

- Map of South Arran MPA, fishing gear restrictions, survey stations and grabs samples visualised

## Usage 

Instructions for reproducing the data cleaning, processing and analysis workflow will be included

## Contents


### Scripts

[project_script.Rmd](project_script.Rmd)

This is the completed workflow of the project, including instructions on how to use the functions I made. 
The main aspects of this script involve reading and cleaning the raw data files, creating the empty BTA fuzzy coding matrix, and conducting the analyses used to answer my thesis objectives. Initially, the script is applied to all data points before being filtered to include only T survey stations, which are more relevant to my study and contain complete data. 
Biodiversity metrics and NMDS plots were then generated to analyse shifts in community structure and address Objective 1. 
Next, the manually completed BTA matrix was imported and analysed to investigate changes in functional traits, addressing Objective 2. The outputs included a PCA, a biplot, and corresponding statistical analyses using linear models. Finally, a species accumulation curve based on Shannon diversity and grab samples was produced to assess current survey effort.

[map.R](map.R)

Script that produces a map of the survey sites within South Arran MPA (Fishing gear restrictions included)


[LRT-model.R](LRT-model.R)

Script that runs Likelihood Ratio Tests (LRTs). Using most plausible complex model was plotted, and LRTs were used to assess the significance of interaction and main effects


### Functions

[All functions created for this project are located and described in R folder](ReadMe_functions.md)


###  Data 

#### Folders

[Survey data](data/survey_data)
The species abundance count

[Environmental information](data/grab_data)
The grab samples with location coordinates, depth (m), sediment type and mean phi (grain size)

[Fishing regulation boundaries](<data/Fishing measures relevant to South Arran>)


#### Additional files

[Species name corrections](<data/unassigned species list.xlsx>)
to be read through before submitting names through the WoRMS database


### Sources for Survey and Grab Data

2012 

[NatureScot Survey 2012 report](https://www.nature.scot/sites/default/files/2025-06/naturescot-commissioned-report-539.pdf)*"species_data_2012.xlsx"* 


2013 

[NatureScot Survey 2013 report](https://www.nature.scot/sites/default/files/2025-06/naturescot-commissioned-report-745.pdf)*"species_data_2013.xlsx"* 


2015

[NatureScot Survey 2015 July report](https://www.nature.scot/doc/naturescot-commissioned-report-945-infaunal-and-psa-analyses-benthic-samples-collected-around-isle)*"species_data_2015_1.xlsx"* 


[NatureScot Survey 2015 September report](https://www.nature.scot/doc/naturescot-commissioned-report-946-infaunal-and-psa-analyses-benthic-samples-collected-south-arran)*"species_data_2015_2.xlsx"*


2022

[2022 Survey data provided by Dr. Kelly Saunders, NatureScot](data/survey_data/species_data_2022.xlsx)*"species_data_2022.xlsx"*

#### Map details

Fishing gear restrictions

[Fishing measures relevant to South Arran](data/Fishing measures relevant to South Arran)

[Shapefiles of fishing management boundaries](Fishing measures relevant to South Arran.shp)

Shapefiles of the survey stations introduced from 2015 onwards

[Survey Stations](data/Fishing measures relevant to South Arran/Survey boxes - 2021 and 2016)

Provided by Dr. Kelly Saunders, NatureScot

#### Taxonomic Information 

Taxonomic names and classifications were validated using the World Register of Marine Species (WoRMS)

(https://www.marinespecies.org/aphia.php)


#### Trait Information

##### BIOTIC Database

[biotic.csv](data/biotic.csv) 

- A csv file that returns species taxonomic information and trait information. Accessed [16/06/2026]: https://api.mba.ac.uk/help_marlin


##### Websites for manual assignment

Databases to search species names for their traits to help with assignment of fuzzy coding.

[Encyclopedia of Life website](https://eol.org/)

[Sea Life Base](https://www.sealifebase.ca/search.php)

[World Register of Marine Species (WoRMS)](https://www.marinespecies.org/aphia.php)

### Additional Folders

#### outputs

**2 output folders**

- *abundance*: BTA calculated with abundance values

- *P_A*: BTA calculated with presence and absence data 

Both folders contain the four outputs returned from the BTA function that creates a empty matrix for trait information and matches species names with the biotic dataframe 


### Credits 

**Supervisors:**
Dr Richard Reeve, University of Glasgow

Dr Kelly Saunders, NatureScot

### Contact Information 
Sarah Stevenson 

3170686S@student.gla.ac.uk 

MSc Project, University of Glasgow 
