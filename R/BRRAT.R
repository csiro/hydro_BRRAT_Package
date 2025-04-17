#https://cran.r-project.org/web/packages/rstantools/vignettes/minimal-rstan-package.html

#' Bayesian Regression Robustness Assessment Test
#'
#' @param Sim Annual time series of simulated values (typically streamflow volume)
#' @param Obs Annual time series of observed values
#' @param Clim Annual time series of independent climate values
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#'
#' @return list with slopes for Sim (beta1), Sim2 if provided (beta2), and rstan model fit object for testing of convergence etc, lambda exponent used for BC transformation
#' @export
#'
#' @examples
#'
#' #generate 2 models of annual streamflow
#' samples <- 50
#' Clim <- sort(pmax(0,rnorm(samples,600,100))) #generate some annual rainfall totals
#' Obs <- pmax(0,(Clim - 150) - 450 * tanh((Clim - 150)/450)) #calculate annual runoff from tanh curve
#' Sim1 <- Obs*seq(1,0.7,length.out=length(Obs))*rnorm(samples,1,0.05) #non-robust model, underestimate the wetter years
#' Sim2 <- Obs*rnorm(samples,1,0.05) #robust model, only random error
#'
#' #apply BRRAT test to both
#' slopes1 <- BRRAT(Sim1,Obs,Clim)
#' slopes2 <- BRRAT(Sim2,Obs,Clim)
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

  if(length(Sim)!=length(Obs) | length(Sim)!=length(Clim)) stop("length of input vectors must be the same")

  if(any(is.na(Sim)) | any(is.na(Obs)) | any(is.na(Clim))) stop("remove NAs")

  if(min(Obs)<0) stop("Obs must be positive")

  #find lambda for box-cox transformation
  lambda <- 1
  out <- MASS::boxcox(lm(Obs~1), plotit=FALSE)
  range <- range(out$x[out$y > max(out$y)-qchisq(0.95,1)/2])

  #if 1 is outside the 95% confidence interval of the best lambda, then transform.
  #Assumes already near normal if 1 is in the range
  if(min(range)>1 | max(range<1)){

    lambda <- out$x[which.max(out$y)]
    Obs <- (Obs^lambda - 1) / lambda
    Sim <- (Sim^lambda - 1) / lambda
    if(!is.null(Sim2)) Sim2 <- (Sim2^lambda - 1) / lambda
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
