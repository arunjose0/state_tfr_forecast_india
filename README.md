# State-level Total Fertility Rate Forecasting in India

This repository contains the code, models, and processed data for the manuscript projecting state-level Total Fertility Rates (TFR) in India using a Bayesian Ornstein-Uhlenbeck model (via Stan).

## Project Structure

- `R/`: Analysis scripts, numbered by execution order (e.g., `01_prepare_stan_data.R` to `07_extract_F_values.R`).
- `stan/`: The Stan model file(s) for the Bayesian model.
- `data/`: Cleaned, derived, and aggregated data files ready for modeling.
- `results/`: Derived outputs (such as summary tables and metrics) and model preparation objects.
- `figures/`: Generated plots and figures.

## Data Sources

The raw data underlying the provided analytical datasets comes from:
- **Sample Registration System (SRS)**: Reports from the Office of the Registrar General, India.
- **National Family Health Survey (NFHS)**: Reports from IIPS and ICF.
- **2011 Census of India**: Office of the Registrar General & Census Commissioner, India.

Due to repository limits, the raw input files are not included. The `data/` folder contains only the finalized, cleaned outputs required for modeling. Please see `data/README.md` for more details.

## How to Reproduce

1. **Environment setup**: The project uses `renv` to manage package dependencies. Upon opening the `state_tfr_forecast_india.Rproj` project, run `renv::restore()` to install the required packages.
2. **Analysis Pipeline**: The scripts in `R/` are designed to be run sequentially:
   - `01_prepare_stan_data.R`: Prepares data structures for Stan modeling.
   - `02_model_fitting_and_diagnostics.R`: Fits the primary Bayesian model.
   - `03_sensitivity_analysis.R`: Runs sensitivity checks on prior configurations.
   - `04_holdout_validation_2018_2023.R`: Fits and validates the model using holdout data.
   - `05_census_srs_comparison.R`: Compares census data to SRS data.
   - `06_national_tfr_trajectory.R`: Estimates national trajectories.
   - `07_extract_F_values.R`: Extracts final output values.

> **Note**: Fitting the Stan models (scripts `02`, `03`, `04`) can be computationally intensive and takes considerable time. The resulting fitted model objects are large (~1-2GB each) and are thus not included in this repository due to GitHub's file size limits. Code lines that save these model fits to disk have been commented out to prevent accidental massive file creation, but you may uncomment them if you wish to save the posterior samples locally.

## License

This project is licensed under the PolyForm Noncommercial License 1.0.0. Free for academic, research, and other noncommercial use. Commercial use requires a separate license — contact [corresponding author email].

