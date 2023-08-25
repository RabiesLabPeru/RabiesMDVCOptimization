# Optimizing the Location of Vaccination Sites to Stop a Zoonotic Epidemic

This repository contains scripts used in the analysis presented in the manuscript titled _Optimizing the Location of Vaccination Sites to Stop a Zoonotic Epidemic._ A description of the repository's contents are provided below. 

## Code files

All scripts can be found in the *code* folder. A description of each script file can be found below:
- **optimize_pcenter_top10max.R** - runs the p-center algorithm
- **optimize_pmedian.R** - runs the p-median algorithm
- **optimize_pprobability.R** - runs the p-probability algorithm
- **regressions_standard.R** - constructs fixed-effects regression models
- **regressions_mixedeffects.R** - constructs mixed-effects regression models
- **regressions_crossvalidation.R** - performs 5-fold cross validation to calculate RMSE and prediction error for the fixed- and mixed-effects regression models
- **figure1_potentialtents.R** - creates figure 1 (map of potential tent locations)
- **figure2_regressions.R** - creates figure 2 (plot of survey data and regression curves)
- **figure3_catchmentsandhistograms.R** - creates figure 3 (maps of vaccination tent catchments, distribution of walking distance, and distribution of workloads)
- **figure4_predictedparticipation.R** - creates figure 4 and supp. fig 1 (maps of houses colored by probability of MDVC participation for different sets of selected vaccination sites)

## Data files

Available data files can be found in the *data* folder. Note that the locations of houses in the study area (demandpoints.csv) and the distance matrix containing walking distances between houses and potential vaccination sites (distancematrixwalking.csv) have been removed due to privacy concerns. 
- **supplypoints.csv** - locations of potential vaccination sites 
- **MR_ALTO_SELVA_ALEGRE.kml** - polygon for Alto Selva Alegre district


