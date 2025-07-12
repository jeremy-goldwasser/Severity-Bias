# Challenges in Estimating Time-Varying Epidemic Severity Rates from Aggregate Data

This repository contains all code, data, and writing for the paper on bias in estimating time-varying severity rates, such as the case fatality rate (CFR).

## Repository Structure and Usage

All data and experimental results can be fully reproduced using the provided `R` scripts. Below is a brief summary of the contents of each folder to help orient users:

- `Code/`
  - `Preprocessing/` – Scripts for downloading and preprocessing real-world data prior to analysis.  
  - `Analysis/` – Scripts starting with `data_for_*` estimate severity rates, whereas the others generate the figures presented in the manuscript.

- `Data/`
  - `Real_data/` – Collected real-world datasets (e.g., hospitalizations, deaths, etc.).  
  - `HFR_estimates/` – Data frames containing estimated hospital fatality rates (HFRs) and related quantities, used in figure generation.

- `Figures/`
  - `Real/` – Figures derived from real-world COVID-19 data.  
  - `Simulated/` – Figures generated from simulation experiments.

## Software Environment

This analysis was conducted using **R version 4.4.1** with the following package versions:

- **lubridate**: 1.9.3  
- **epidatr**: 1.2.0  
- **extraDistr**: 1.10.0  
- **tidyverse**: 2.0.0  
- **patchwork**: 1.3.0
