#' Bayesian Regression Robustness Assessment Test
#'
#' Function to apply the BRRAT. A linear regression is fit between the error between the transformed `Sim`
#' and `Obs` data and the corresponding `Clim` data, to test for a dependent
#' relationship between the error and climate, indicating a model that is not
#' robust to the input climate data. Stan is used to fit the linear regression,
#' and hence a distribution of slope values are returned.
#'
#' Multiple catchments are supported through random factors in the hierarchical
#' linear regression model,indicated by different values in an `ID` column.
#'
#' Multiple replicates of the Simulated time series are supported, matched based
#' on the `year` column.
#'
#' @param Sim Data frame type containing time series of simulated values (typically annual streamflow).
#' @param Obs Data frame type containing time series of observed values (typically annual streamflow and rainfall).
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains, cores).
#'
#' `Sim` requires column names of `year` and `Q`.
#' `year` is used to match the simulated an observed values, and
#' `Q` are the simulated values to test, typically streamflow volume/depth
#'#'
#' `Obs` requires column names of `year`, `P` and `Q`.
#' `year` is the year of observation, matched to `year` in `Sim` to calculate errors
#' `Q` is the observed values to compare to the `Q` column from `Sim`, using a log ratio
#' `P` is the independent climate variable, to test if there is a relationship in the errors
#'
#' An optional `ID` column can be included to indicate multiple sites to be used
#' as random factors in the test. If an `ID` column is used, it is required in both `Sim` and `Obs`.
#'
#' @return a list with members:
#' \describe{
#' \item{fit}{rstan model  object for testing of convergence etc}
#' \item{data}{list of data in stan format used to fit the model, after transformation}
#' \item{id_levels}{order of random factors if used}
#' \item{beta_draws}{samples of slope values from the posterior distribution}
#' \item{statistics}{statistics indicating the probability that the slope is: greater than 0 `prob_beta_gt0`, less than 0 `prob_beta_lt0`, and non-zero `p_nonzero`}
#'}
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#' #generate 2 models of annual streamflow
#' samples <- 50
#' #generate some annual rainfall totals
#' Clim <- sort(pmax(0,rnorm(samples,600,100)))
#' #calculate annual runoff from tanh curve
#' observed <- pmax(0,(Clim - 150) - 450 * tanh((Clim - 150)/450))
#' Obs <- data.frame(P=Clim,Q=observed,year=1:length(Clim))
#'
#' #non-robust model, underestimate the wetter years
#' Sim1 <- observed*seq(1,0.7,length.out=length(observed))*rnorm(samples,1,0.05)
#' Sim1 <- data.frame(Q=Sim1,year=1:length(Clim))
#'
#' #robust model, only random error
#' Sim2 <- observed*rnorm(samples,1,0.05)
#' Sim2 <- data.frame(Q=Sim2,year=1:length(Clim))
#'
#' #apply BRRAT test to both, reduce chains and iters to speed up example.
#' fit1 <- BRRAT(Sim1,Obs,chains=1,iter=500)
#' fit2 <- BRRAT(Sim2,Obs,chains=1,iter=500)
#'
#' dat <- data.frame(b = c(fit1$beta_draws,fit2$beta_draws),
#' id = c(rep("fit1",length(fit1$beta_draws)),rep("fit2",length(fit2$beta_draws)))
#' ggplot(dat)+geom_density(aes(b,fill=id)) #fit2 is centered on 0, fit1 has a negative slope
#' #can also be viewed using the plotting function
#' plot_BRRAT_population(fit1)
#' plot_BRRAT_population(fit2)
#'
#' #Summarise statistics of the mean can be viewed using
#' summarise_BRRAT(fit1)
#' summarise_BRRAT(fit2)
#'}
#'

BRRAT<-function(Sim, Obs, ...){

  #check inputs are OK
  if(!("data.frame" %in% class(Sim))){
    stop("Sim must be of type that inherits data.frame (data.frame, tibble, etc)")
  }

  if(!("data.frame" %in% class(Obs))){
    stop("Obs must be of type that inherits data.frame (data.frame, tibble, etc)")
  }

  expected_cols <- c("year", "Q", "P")
  if(!all(expected_cols %in% names(Obs))){
    stop(paste("Obs input requires columns of:\nyear (index of observations to match to Sim)\nP (independent climate variable)\nQ (observed flow response to evaluate against Sim"))
  }

  if(any(is.na(Sim)) | any(is.na(Obs))) stop("remove NAs")

  if(min(nrow(Sim), nrow(Obs))<3) stop("data frame too small, not enough rows")

  if(min(Obs$Q)<=0) stop("Obs$Q must be positive")

  if(!("ID" %in% names(Obs))){
    warning("No ID column in Obs. Assuming 1 site, adding dummy ID column")
    Obs <- Obs |> dplyr::mutate(ID = "ID")
    Sim <- Sim |> dplyr::mutate(ID = "ID")
  }

  expected_cols <- c("year", "Q", "ID")
  if(!all(expected_cols %in% names(Sim))){
    stop(paste("Sim input requires columns of:\nyear (index of observations to match to Obs)\nQ (Simulated flow response to evaluate against Obs)"))
  }

  #Standardise Rain
  meanP <- Obs |> dplyr::group_by(.data$ID) |> dplyr::summarise(Pm = mean(.data$P),
                                              SD = stats::sd(.data$P))

  Obs <- Obs |> dplyr::left_join(meanP,by="ID") |>
    dplyr::mutate(P=(.data$P-.data$Pm)/.data$SD) |>
    dplyr::select(dplyr::all_of(c("year","P","Obs"="Q","ID")))

  #join data, calculate error as log ratio

  df <- Sim |>
    dplyr::left_join(Obs,by=c("ID","year")) |>
    dplyr::filter(!is.na(.data$Obs)) |> #only years with matching obs
    dplyr::mutate(Qerr = log(.data$Q/.data$Obs),#log ratio
           ID = factor(.data$ID)) |>
    dplyr::select(dplyr::all_of(c("P","Qerr","ID")))

  #setup stan list
  x <- df$P
  y <- df$Qerr
  group <- as.integer(df$ID)
  id_levels <- levels(df$ID)
  J <- length(id_levels)
  N <- nrow(df)

  stan_data <- list(N = N, x = as.vector(x), y = as.vector(y), J = J, group = group)

  if(J==1){
    fit <- rstan::sampling(stanmodels$lm, data = stan_data, ...)
  }else{
    fit <- rstan::sampling(stanmodels$lm_random_id, data = stan_data, ...)
  }


  # Posterior summaries
  summ <- rstan::summary(fit)$summary
  summ_df <- as.data.frame(summ)
  summ_df$param <- rownames(summ)

  # Slope evidence
  beta_draws <- rstan::extract(fit, pars = 'beta')[[1]]
  prob_beta_gt0 <- mean(beta_draws > 0)
  prob_beta_lt0 <- mean(beta_draws < 0)
  p_nonzero <- 2 * min(prob_beta_gt0, prob_beta_lt0)

  out <- list(
    fit = fit, data = stan_data,
    id_levels = id_levels,
    summary = summ_df,
    beta_draws = beta_draws,
    prob_beta_gt0 = prob_beta_gt0, prob_beta_lt0 = prob_beta_lt0, p_nonzero = p_nonzero
  )
  class(out) <- 'BRRAT'

  return(out)
}

