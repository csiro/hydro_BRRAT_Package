
#' Summarise global slope and per-site coefficients into tidy tables
#'
#' @title Summarise population and site-level posterior summaries
#' @description
#' Computes tidy posterior summaries for the **population slope** (`beta`) and
#' **per-site intercepts and slopes** (from `intercept_id` and `slope_id` generated
#' quantities) returned by the Stan fit. The function reports posterior means,
#' standard deviations, equal-tail credible intervals at a chosen level, and a
#' two-sided tail probability for non-zero slope at the population level.
#'
#' @details
#' The two-sided tail probability is computed as:
#' \deqn{p_{\mathrm{nonzero}} = 2\min\{ \Pr(\beta > 0 \mid \mathrm{data}), \Pr(\beta < 0 \mid \mathrm{data}) \}}{
#' p_nonzero = 2 * min(Pr(beta > 0 | data), Pr(beta < 0 | data))}
#' which is analogous to a two-sided p-value in spirit, but is a **Bayesian posterior tail
#' probability**. It shrinks toward 0 as the posterior mass concentrates away from zero.
#'
#' @param fit_obj A fitted object returned by [BRRAT()], of class `BRRAT`.
#' @param level Credible interval level (default `0.95`). The lower and upper
#' intervals are computed as `(1 - level) / 2` and `1 - (1 - level) / 2`.
#'
#' @return A list with one or two data frames:
#' \itemize{
#'   \item `global`: one row summarising the population slope `beta` (mean, sd,
#'   lower/upper credible interval, `prob_gt0`, `prob_lt0`, and `p_nonzero`).
#'   \item `per_id`: if applicable (more than one site `ID` in the data) one
#'   row per site `ID` with posterior summaries for
#'   `intercept_id` and `slope_id` (mean, sd, lower/upper credible intervals).
#' }
#'
#' @section Notes:
#' Assumes the Stan model emits `intercept_id` and `slope_id` in `generated quantities`,
#' i.e., per-site fully-composed coefficients:
#' \deqn{\mathrm{intercept\_id}_j = \alpha + b_{0j},\quad \mathrm{slope\_id}_j = \beta + b_{1j}}{
#' intercept_id[j] = alpha + b0[j]; slope_id[j] = beta + b1[j]}
#'
#' @examples
#' \dontrun{
#' res <- BRRAAT(Sim,Obs)
#' tabs <- summarise_BRRAT(res, level = 0.95)
#' tabs$global
#' head(tabs$per_id)
#' }
#'
#' @export
summarise_BRRAT <- function(fit_obj, level = 0.95) {
  stopifnot(inherits(fit_obj, 'BRRAT'))
  sf   <- fit_obj$fit
  q_lo <- (1 - level) / 2
  q_hi <- 1 - q_lo

  # Population slope beta
  beta_draws <- rstan::extract(sf, pars = 'beta')[[1]]
  global <- data.frame(
    param     = 'beta',
    mean      = mean(beta_draws),
    sd        = stats::sd(beta_draws),
    q_lo      = stats::quantile(beta_draws, probs = q_lo),
    q_hi      = stats::quantile(beta_draws, probs = q_hi),
    prob_gt0  = mean(beta_draws > 0),
    prob_lt0  = mean(beta_draws < 0),
    p_nonzero = 2 * min(mean(beta_draws > 0), mean(beta_draws < 0))
  )

  if(fit_obj$data$J>1){
    # Per-ID intercepts and slopes
    intercept_draws <- rstan::extract(sf, pars = 'intercept_id')[[1]]  # D x J
    slope_draws     <- rstan::extract(sf, pars = 'slope_id')[[1]]      # D x J
    id_levels       <- fit_obj$id_levels

    per_id <- data.frame(
      ID               = id_levels,
      intercept_mean   = apply(intercept_draws, 2, mean),
      intercept_sd     = apply(intercept_draws, 2, stats::sd),
      intercept_q_lo   = apply(intercept_draws, 2, stats::quantile, probs = q_lo),
      intercept_q_hi   = apply(intercept_draws, 2, stats::quantile, probs = q_hi),
      slope_mean       = apply(slope_draws, 2, mean),
      slope_sd         = apply(slope_draws, 2, stats::sd),
      slope_q_lo       = apply(slope_draws, 2, stats::quantile, probs = q_lo),
      slope_q_hi       = apply(slope_draws, 2, stats::quantile, probs = q_hi)
    )
  }else{
    per_id <- NULL
  }

  list(global = global, per_id = per_id)
}

#' Export summary tables to CSV files
#'
#' @title Export global and per-site summaries to CSV
#' @description
#' Writes the output of [summarise_BRRAT()] to two CSV files: one for the
#' population slope (`file_global`) and one for per-site coefficients (`file_per_id`).
#'
#' @param fit_obj A fitted object returned by [BRRAT()], of class `BRRAT`.
#' @param file_global Path/filename for the global table CSV (default `'global.csv'`).
#' @param file_per_id If applicable to the model, Path/filename
#' for the per-site table CSV (default `'per_id.csv'`).
#' @param level Credible interval level passed to [summarise_BRRAT()] (default `0.95`).
#'
#' @return (Invisibly) a list with the two data frames written:
#' `list(global = <data.frame>, per_id = <data.frame>)`.
#'
#' @examples
#' \dontrun{
#' res <- BRRAT(Sim,Obs)
#' export_BRRAT_csv(res,
#'   file_global = "qerr_global.csv",
#'   file_per_id = "qerr_per_id.csv",
#'   level = 0.95
#' )
#' }
#'
#' @export
export_BRRAT_csv <- function(fit_obj,
                                file_global = 'lobal.csv',
                                file_per_id = 'per_id.csv',
                                level = 0.95) {
  tabs <- summarise_BRRAT(fit_obj, level = level)
  utils::write.csv(tabs$global,  file_global,  row.names = FALSE)
  if(!is.null(tabs$per_id)){
    utils::write.csv(tabs$per_id,  file_per_id,  row.names = FALSE)
  }
  invisible(tabs)
}
