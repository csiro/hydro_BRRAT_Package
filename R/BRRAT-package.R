#' The 'BRRAT' package.
#'
#' @description Function to apply the Bayesian Regression Robustness Assessment Test,
#' to evaluate if model simulations are independent of the climate input used,
#' suggesting a robust model.
#'
#' This package makes use of Stan to fit the linear regression model,
#' and the precompiled model created for the package with rstantools.
#'
#' @docType package
#' @name BRRAT-package
#' @aliases BRRAT-package
#' @useDynLib BRRAT, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import MASS
#' @import ggplot2
#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom stats lm qchisq
#'
#' @references
#' Stan Development Team (NA). RStan: the R interface to Stan. R package version 2.32.6. https://mc-stan.org
#'
#' Gibbs et al. (2025)
#'
"_PACKAGE"
NULL
