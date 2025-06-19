# Challenges in Estimating Time-Varying Epidemic Severity Rates from Aggregate Data

This repository contains all code, data, and writing for the paper on bias in estimating time-varying severity rates, such as the case fatality rate (CFR).

## Usage

All data and experimental results can be fully reproduced using the provided `R` scripts.

1. **Preprocess real-world data**  
   Download and preprocess datasets using scripts in `Code/Preprocessing`.

2. **Compute severity rates**  
   Use scripts in `Code/Analysis` starting with `data_for_*` to estimate severity rates, saving results in `Data/HFR_estimates`.

3. **Generate figures**  
   Create figures using the other scripts in `Code/Analysis`. These scripts run quickly.

## Software Environment

This analysis was conducted using **R version 4.4.1** with the following package versions:

- **lubridate**: 1.9.3  
- **epidatr**: 1.2.0  
- **extraDistr**: 1.10.0  
- **tidyverse**: 2.0.0  
- **patchwork**: 1.3.0
