# Bayesian Regression Robustness Assessment Test (BRRAT)

<!-- badges: start -->

[![R-CMD-check](https://github.com/csiro/hydro_BRRAT_Package/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/csiro/hydro_BRRAT_Package/actions/workflows/R-CMD-check.yaml)[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.17957953-blue)](https://doi.org/10.5281/zenodo.17957953)

<!-- badges: end -->

## Overview

The Robustness Assessment Test (RAT) developed by [Nicolle et al. (2021)](https://hess.copernicus.org/articles/25/5013/2021/) tests for a relationship between the relative bias in annual streamflow and climatic variables, 
e.g. rainfall, temperature or humidity. If a statistically significant relationship is identified it is suggested that the rainfall runoff model that generated the annual streamflow series cannot be safely used for climate change impact studies. 

This package implements an extension that compares the slope of the linear regression line between relative model error and a rainfall index, 
with a slope closer to zero representative of a more robust model. This extension allows for a improvements in model structure or calibration method to be tested more easily compared to the original RAT pass/fail test. 
The regression parameters are fit using Bayesian regression (making use of the rstan package), with the modified test referred to as the Bayesian Regression Robustness Assessment Test (BRRAT). 

A hierarchical linear regression approach is used if multiple catchments are considered simultaneously.

## Getting Started

Package has been submitted to CRAN. Once on CRAN, install with install.packages("BRRAT")

To install from github source, the easiest approach is: devtools::install_github("csiro/hydro_BRRAT_Package")

## Usage 

```r
# Generate 2 models of annual streamflow
samples <- 50

# Generate some annual rainfall totals
Clim <- sort(pmax(0, rnorm(samples, 600, 100)))

# Calculate annual runoff from tanh curve
observed <- pmax(0, (Clim - 150) - 450 * tanh((Clim - 150)/450))
Obs <- data.frame(P = Clim, Q = observed, year = 1:length(Clim))

# Non-robust model, underestimates the wetter years
Sim1 <- observed * seq(1, 0.7, length.out = length(observed)) * rnorm(samples, 1, 0.05)
Sim1 <- data.frame(Q = Sim1, year = 1:length(Clim))

# Robust model, only random error
Sim2 <- observed * rnorm(samples, 1, 0.05)
Sim2 <- data.frame(Q = Sim2, year = 1:length(Clim))

# Apply BRRAT test to both (chains and iter reduced to speed up example)
fit1 <- BRRAT(Sim1, Obs, chains = 1, iter = 500)
fit2 <- BRRAT(Sim2, Obs, chains = 1, iter = 500)

# Compare posterior distributions of slope — fit2 centred on 0, fit1 has negative slope
dat <- data.frame(
  b  = c(fit1$beta_draws, fit2$beta_draws),
  id = c(rep("fit1", length(fit1$beta_draws)), rep("fit2", length(fit2$beta_draws)))
)
ggplot(dat) + geom_density(aes(b, fill = id))

# Or use the built-in plotting function
plot_BRRAT_population(fit1)
plot_BRRAT_population(fit2)

# Summarise posterior statistics
summarise_BRRAT(fit1)
summarise_BRRAT(fit2)
```

## Citation

Gibbs MS, Wang B, Crosbie R, Montazeri M, Vaze J and Yang A (2026) Improving distributed hydrological model robustness through multi-variable calibration. Hydrological Processes, In Press.

## Contribute

Feel free to add pull requests, raise issues here on github or email.



