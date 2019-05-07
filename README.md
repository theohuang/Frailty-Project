# Frailty-Project

## Description

This respository contains code to apply a frailty model to improve Mendelian model risk predictions. The Simulations folder contains simulation results for improving BRCAPRO breast cancer risk prediction under 6 settings: 2 types of BRCA1/2 mutation prevalences (low risk, high risk) and 3 types of data-generating frailty distributions (discrete uniform; bivariate normal with mean (0, 0), correlation 0, and variances 0.3; bivariate normal with mean (0, 0), correlation 0, and variances 2). All code is written by Theodore Huang except where specified.

## Running the code
The way to reproduce the simulation results is as follows:

1. Generate the data using one of the FrailtySimDatGen().R files (choice of prevalence and data-generating frailty distribution).
2. Obtain the family-specific frailty distributions and risk predictions:
    1. First run the corresponding FrailtySim_RP_().R file
    2. Then gather the data files together using code in Frailty Simulation Analysis.R.
    
Note that running the code requires the BayesMendel R package, which can be found at https://projects.iq.harvard.edu/bayesmendel/bayesmendel-r-package.
    
## The files are all described below:

### General files

* Estimating Functions Discrete.R: functions used for the frailty method
* Parameters Discrete.R: getting the baseline penetrances, survival functions, and hazard functions
* Peeling Functions folder: functions for peeling-paring

### Simulations Folder

#### Subscripts

* "Dis" means discrete uniform distribution
* "V03" means bivariate normal distribution with mean 0, correlation 0, and variances of 0.3
* "V2" means bivariate normal distribution with mean 0, correlation 0, and variances of 2
* "hr" means high-risk families (using higher prevalences for AJ and non-AJ families)

#### Files

* Frailty Simulation Functions.R: functions to generate data
* FrailtySimDatGen().R: generating data
* FrailtySim_RP_().R: obtains family-specific frailty distributions and risk predictions for simulations
* Frailty Simulation Analysis.R: code to obtain final results for simulation
* simdat_().RData: simulated data
* simdat_bc_().RData: subset of families with female probands free of breast cancer in simulated data

#### Simulations subfolders:

* Family Frailty Distribution folder: has all the family-level frailty distributions
* Risk Predictions folder: has all the risk predictions
* Cluster Files folder: files to run on Harvard FAS cluster
* Generating Families Functions folder: has the code to generate families
