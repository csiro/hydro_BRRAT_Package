#' Bayesian Regression Robustness Assessment Test
#'
#' Function to apply the BRRAT. The input `Sim` and `Obs` time series are first
#' transformed to near normal distributions using a boxcox transformation, with
#' the lambda exponent fit to the `obs` data.
#'
#' A linear regression is fit between the error between the transformed `Sim`
#' and `Obs` data and the corresponding `Clim` data, to test for a dependent
#' relationship between the error and climate, indicating a model that is not
#' robust to the input climate data. Stan is used to fit the linear regression,
#' and hence a distribution of slope values are returned.
#'
#' @param Sim Time series of simulated values (typically annual streamflow volume)
#' @param Obs Time series of observed values
#' @param Clim Time series of independent climate values
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#'
#' @return a list with members:
#' \describe{
#' \item{beta}{samples of slope values from the posterior distribution}
#' \item{fit}{rstan model  object for testing of convergence etc}
#' \item{lambda}{value of exponent used for BC transformation}
#'}
#'
#' @export
#'
#' @examples
#'
#' #generate 2 models of annual streamflow
#' samples <- 50
#' #generate some annual rainfall totals
#' Clim <- sort(pmax(0,rnorm(samples,600,100)))
#' #calculate annual runoff from tanh curve
#' observed <- pmax(0,(Clim - 150) - 450 * tanh((Clim - 150)/450))
#' #non-robust model, underestimate the wetter years
#' Sim1 <- observed*seq(1,0.7,length.out=length(observed))*rnorm(samples,1,0.05)
#' #robust model, only random error
#' Sim2 <- observed*rnorm(samples,1,0.05)
#'
#' #apply BRRAT test to both
#' slopes1 <- BRRAT(Sim1,observed,Clim)
#' slopes2 <- BRRAT(Sim2,observed,Clim)
#'
#' #compare slopes of regression between box-cox transformed errors and climate
#' dat <- data.frame(beta = c(slopes1$beta, slopes2$beta),
#'                   model = c(rep("Sim1",length(slopes1$beta)),
#'                   rep("Sim2",length(slopes2$beta))))
#' ggplot2::ggplot(dat)+
#' ggplot2::geom_density(ggplot2::aes(beta,fill=model), alpha=0.5)+
#' ggplot2::geom_vline(xintercept=0, linetype="dashed")
#'

BRRAT<-function(Sim, Obs, Clim, ...){

  if(length(Sim)!=length(Obs) | length(Sim)!=length(Clim)){
    stop("length of input vectors must be the same")
  }

  if(any(is.na(Sim)) | any(is.na(Obs)) | any(is.na(Clim))) stop("remove NAs")

  if(min(Obs)<0) stop("Obs must be positive")

  #find lambda for box-cox transformation
  lambda <- 1
  dframe <- data.frame(x = Obs)
  lm_m <- stats::lm(x~1, data=dframe)
  out <- MASS::boxcox(lm_m, data=dframe, plotit=FALSE)
  range <- range(out$x[out$y > max(out$y)-stats::qchisq(0.95,1)/2])

  #if 1 is outside the 95% confidence interval then transform.
  #Assumes already near normal if 1 is in the range
  if(min(range)>1 | max(range<1)){

    lambda <- out$x[which.max(out$y)]
    Obs <- (Obs^lambda - 1) / lambda
    Sim <- (Sim^lambda - 1) / lambda
  }

  Qerr <- Sim/Obs - 1
  Clim_Index <- Clim/mean(Clim) - 1

  stan_data <- list(
      N = length(Obs),
      x = Clim_Index,
      y = Qerr
    )

  fit <- rstan::sampling(stanmodels$lm, data = stan_data, ...)
  posterior_samples <- rstan::extract(fit)

  return_list <- list(beta = posterior_samples$beta,
                     fit = fit,
                     lambda = lambda)

  return(return_list)
}
