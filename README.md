# Welcome

Hello! This repository contains all of the lovely code for my dissertation. Data from this project come from the Couples Coping Study (CCS). Data was collected and cleaned in the early days of my R journey, so stuff may be a bit messy. This README file reminds future me of where important stuff comes from. 

# Quick reference guide

- data = datafiles used in this project; NOTE: b/c I don't have permissions to openly share the data from this study, this folder is ignored from the public repo but is available for request. For future me, the data folder is in the local repo.
- documentation = any files that contain information about the data or were used for cleaning that are NOT part of dissertation analyses 
- output = all output is automatically embedded in the markdown documents for analyses, this folder refers to any output used for table creation/output reporting in manuscripts; generally, this workflow is to create a dataframe with the relevant information and save it as a .csv so all numbers/data is available for manual table creation 
- fits = b/c this project relies on Bayesian models that can take some time to fit, these are saved in this folder and the code is written such that models are only re-run if there are changes to them 

# Note about data cleaning

The "CCS_data_cleaning_syntax_2021.01.10.Rmd" file is the final data cleaning file that was used to generate a (mostly) completely clean dataset. There were a couple of things in there like IPV and TLEQ scoring that were left for the future b/c of the multitude of ways these scores can be calculated. Thankfully, raw data is still available in the cleaned dataset so you can backtrace anything that is needed. This RMarkdown file is in the "documentation" folder. 

Naturally, this file creates the "CCS_data_final_2021.01.10.rds" file that contains the final, cleaned dataframe from this syntax. *THIS IS YOUR STARTING POINT!* You put this file in the "data" subfolder of this directory so it is easily accessible. 