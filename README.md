# Code from: *If you’re rare, should I care? How imperfect detection changes relationships between biodiversity and global change drivers*

## General description

This is the code for a paper exploring the importance of accounting for rare species in biodiversity metrics and their relationships with global change drivers. We accounted for detection error using multi-species occupancy and abundance models and then used results from these models for "true" abundance or presence to derive diversity metrics. We then examined the relationships with global change drivers and biodiversity using beta regressions.

We performed this modeling process on three case studies available from published sources:

| Dataset                                          | Years     | Number of sites | Number of species | Source                                                                                                              |
|--------------------------------------------------|-----------|-----------------|-------------------|---------------------------------------------------------------------------------------------------------------------|
| Konza Prairie LTER passerine birds               | 1981-2009 | 11              | 78                | [Boyle 2023](https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-knz.26.12)                           |
| Sevilleta LTER grasshoppers                      | 1992-2019 | 60              | 46                | [Lightfoot 2021](https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-sev.106.214969)                  |
| Petrified Forest National Park understory plants | 2007-2022 | 70              | 84                | [Swan and Ploughe 2023](https://irma.nps.gov/DataStore/Reference/Profile/2300890)                                   |

## Folder structure

### examples

The examples folder contains subfolders for each process in our data analysis for each example dataset. 

#### code

This folder is broken up by coding outcome, including custom functions we created for this project, data prep, model code and wrapper scripts for JAGS, and code to recreate visualizations. Each subfolder contains processes distinct to each dataset

#### 00_functions

This folder contains several scripts that are sourced in other scripts in order to perform custom functions we created for this project.

#### 01_data_prep

This folder contains scripts for compiling data for each case study and making it into list format for each JAGS model. This folder is subdivided by case study and then each case study contains two to three data prep scripts: one for preparing data for the model accounting for imperfect detection; one to compute diversity metrics from posterior samples; and the last preparing data for the beta regression model

#### 02_models

This folder contains models and their wrapper scripts for running each of the two models (MSAM/MSOM and beta regression) for each case study. It also contains template models for both of these modeling steps. This folder is first divided into the first and second models (MSAM/MSOM & beta regression) and then each of these has subfolders for each case study. These subfolders also contain templates that can be applied to new datasets. Further, each case study also includes a wrapper script that uses the `jagsUI` package to run JAGS in R

#### 03_visualizations

This folder contains code to recreate the figures in the main text and many of the figures in the Supporting Information. Which figures are created by each script are indicated in the name of that script.

### data_raw

This folder contains the raw data (if it was small enough to upload to GitHub) for each of the case studies. When raw data were too large, we provide a link to the data repositories where we accessed data.This folder is subdivided by case study. It also contains raw data used to generate Figure 1, which includes more datasets than we analyzed in this study.

### data_output

This folder contains all data generated by our data prep scripts and which we used for running models and generating summaries for each model. This folder is subdivided by case study and then by which model data prep process generated the data (MSAM/MSOM or beta regression). Further, we have subdivided those folders into whether the data generated is the data list used to run the JAGS model or whether it is additional data (e.g., a tidy dataset of all data for that model) used later in summarizing and linking site and species data back to model outputs.This folder also contains all the resulting model summaries generated by converged JAGS models that we used in results and figures.

