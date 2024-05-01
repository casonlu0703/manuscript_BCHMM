# Overview

This GitHub repo provides R/Stan codes and example simulated data developed for the manuscript - A Bayesian Circadian Hidden Markov Model to Infer Rest-Activity Rhythms Using 24-hour Actigraphy Data.

# Directory

## Data

All example data are stored in the data folder.

1. sample_data.rds: example NHANES dataset. The variables are,
- SEQN: subject ID.
- day: day of the wearable recordings.
- fivemin: indicator variable of 5-min intervals. One 24-hour of data contains 288 data points.
- activity_avg: 5-min averaged actigraphy data (in unit of MIMS/min).

2. true_param.rds: true parameters for the simulation.

## Codes

01_model_fit.R: Fit BCHMM using R Stan

02_model_diagnosis.R: Posterior diagnosis for model fits, including calculating diagnostic statistics and plotting traceplots.

03_extract_RAR_parameter.R: Extract posterior estimates such as time-varying transition probabilities and calculate RAR parameters such as Rhythmic Index and Rest Amount.

function_stan_for_bchmm.STAN: STAN file containing codes to fit BCHMM.

function_calculate_circadian_parameters.R: R function to calculate RAR parameters.

function_figure_day_profile.R: Plot 24-hour rest-activity profiles using time-varying transition probabilities fitted by BCHMM.

simulation_fit_simulated_date.R: Generate simulated data and perform BCHMM.

function_simulate data.R: R function to generate simulated data under three scenarios described in the simulation study.
