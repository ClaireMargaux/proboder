# PROBODER

## Author: Claire Descombes

State-space model for joint prediction of transmission rate and reproduction number from SEIR or SEIRD model and public health data.

This R project is part of my master's thesis, which can be read on [Overleaf](https://www.overleaf.com/read/zvscscfpxfys#69aab2).

The algorithm I developed in R is heavily inspired by the article by J. Schmidt, P. Hennig et al. _A Probabilistic State Space Model for Joint Inference from Differential Equations and Data_, which can be found [here](https://proceedings.neurips.cc/paper/2021/hash/6734fa703f6633ab896eecbdfad8953a-Abstract.html).

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Data](#data)
- [Algorithm](#algorithm)
- [Parameters](#parameters)
- [Contributing](#contributing)
- [License](#license)
- [Credits](#credits)
- [Contact](#contact)

## Installation

Install [R](https://cran.r-project.org/), [RStudio](https://posit.co/products/open-source/rstudio/).

The inference is done using the Rmarkdown workflows:
- workflow_sim_data.Rmd
- workflow_real_data.Rmd

For the workflows to work, you need to download the Functions folder.

To store the data, you need to create a Results folder in the working directory.

## Usage

Infer the transmission rate and reproduction number from public health data (cases) using a SEIR or SEIRD model. Inference through a state-space model using Gauss-Markov priors and (extended) Kalman filtering for the prediction and update steps.

## Data

Data from the [Covid-19 Switzerland Dashboard](https://www.covid19.admin.ch/de/epidemiologic/case) of the Federal Office of Public Health (FOPH):
- COVID19Cases_geoRegion_w.csv: Weekly record timelines by geographical region for cases from February 24, 2020, to January 1, 2023.
- COVID19Cases_geoRegion.csv: Daily record timelines by geographical region for cases from February 24, 2020, to January 1, 2023.

The raw data as well as the R script used to import and format the data are in the Data folder.
 
## Algorithm

This repository contains an implementation of an (extended) Kalman filter for inference on epidemiological models, specifically SEIRD and SEIR. The algorithm integrates state-space modelling with ordinary differential equations (ODEs) to estimate and predict the dynamics of disease spread based on observed data. The SSM is based on three distinct steps:

1. **Prediction step**
   - Objective: Forecast the future state of the system.
   - How: Uses the linear stochastic differential equations (SDEs) to model the system's dynamics and predict how the state will evolve over time. This step treats the system as Gauss-Markov processes.

3. **Update on observations step**
   - Objective: Refine the predicted state using actual measurement data.
   - How: Applies the observation model to compare the predicted state with actual observations. The Kalman filter adjusts the state estimate and uncertainty based on the difference between predicted and observed data.

4. **Update on ODE step**
   - Objective: Incorporate the SEIR or SEIRD equations in the system.
   - How: Utilizes a system of ordinary differential equations (ODEs) to update the state estimate. The ODE solver, integrated with an extended Kalman filter, helps refine the state and uncertainty estimates by enforcing the satisfaction of the ODEs.

## Parameters

- Initial mean and covariance
- Latency rate (λ)
- Recovery rate (γ)
- Mortality rate (η)
- Drift matrices
- Dispersion matrices
- Length scale (ℓ)
- Observation Noise (R)
- Noise of the Wiener processes
- Overall, OBS and ODE time grids
- Jitter value on innovation covariance (ODE update)
  
## Contributing

Don't hesitate to file issues if you find code bugs.

## License

You can use, copy, modify, and distribute this code for personal, educational, and non-commercial purposes.

**Restrictions**
- Publications: You may not use this code or any derivative works for publications, including but not limited to academic papers, articles, or books, without obtaining explicit permission from the original author.
- Attribution: Any use of this code must include proper attribution to the original author.

**Disclaimer**
This code is provided "as is" without warranty of any kind, express or implied. The author is not responsible for any damages arising from using this code.


## Credits

Many thanks to @Ginsbourger and @jriou.

## Contact

Project maintainer(s): Claire Descombes (@ClaireMargaux)
