
#' Plot population-level regression Qerr ~ P with uncertainty band
#'
#' @title Population-level line and credible band
#' @description
#' Plots the **population-level** regression line \eqn{\alpha + \beta x} and its
#' equal-tail credible band by transforming the x-grid to the centered/scaled space
#' used during fitting. The plot overlays a (potentially downsampled) scatter of the
#' original observations for context.
#'
#' @param fit_obj A fitted object returned by [BRRAT()], of class `BRRAT`.
#' @param n_grid Integer; number of x points for the prediction grid (default `200`).
#' @param level Credible interval level for the band (default `0.95`).
#'
#' @return A `ggplot2` object showing the population regression line and credible band.
#'
#' @details
#' The function extracts posterior draws of `alpha` and `beta` from the Stan fit and
#' evaluates \eqn{\mu(x) = \alpha + \beta x} along a grid of x values on the original
#' scale, internally converting to the centered/scaled space via
#' \eqn{x_{\mathrm{std}} = (x - x_{\mathrm{center}}) / x_{\mathrm{scale}}}.
#'
#' @examples
#' \dontrun{
#' res <- BRRAT(Sim,Obs)
#' p_pop <- plot_BRRAT_population(res, n_grid = 200, level = 0.95)
#' print(p_pop)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_ribbon geom_line labs theme_minimal
#' @export
plot_BRRAT_population <- function(fit_obj, n_grid = 200, level = 0.95) {
  stopifnot(inherits(fit_obj, 'BRRAT'))

  draws <- rstan::extract(fit_obj$fit, pars = c('alpha','beta'))
  alpha <- draws$alpha
  beta  <- draws$beta

  x_min  <- min(fit_obj$data$x, na.rm = TRUE)
  x_max  <- max(fit_obj$data$x, na.rm = TRUE)
  x_std <- seq(x_min, x_max, length.out = n_grid)

  mu_draws <- sapply(x_std, function(xi) alpha + beta * xi)
  mu_mean  <- apply(mu_draws, 2, mean)
  q_lo     <- (1 - level) / 2
  q_hi     <- 1 - q_lo
  mu_lo    <- apply(mu_draws, 2, stats::quantile, probs = q_lo)
  mu_hi    <- apply(mu_draws, 2, stats::quantile, probs = q_hi)

  df_plot <- data.frame(P = x_std, mu_mean = mu_mean, mu_lo = mu_lo, mu_hi = mu_hi)

  data_sc <- data.frame(x=fit_obj$data$x,y=fit_obj$data$y)
  if (nrow(data_sc) > 20000) {
    set.seed(1)
    data_sc <- data_sc[sample(seq_len(nrow(data_sc)), 20000), ]
  }

  if (!requireNamespace('ggplot2', quietly = TRUE)) stop('ggplot2 is required.')

  ggplot2::ggplot() +
    ggplot2::geom_point(data = data_sc,
                        ggplot2::aes(x = .data$x, y = .data$y), alpha = 0.15, size = 0.6) +
    ggplot2::geom_ribbon(data = df_plot,
                         ggplot2::aes(x = .data$P, ymin = .data$mu_lo, ymax = .data$mu_hi),
                         fill = '#1f77b4', alpha = 0.15) +
    ggplot2::geom_line(data = df_plot, ggplot2::aes(x = .data$P, y = .data$mu_mean),
                       colour = '#1f77b4', linewidth = 1) +
    ggplot2::labs(x = 'P (normalized rainfall)', y = 'Qerr (transformed)',
                  title = sprintf('Population regression: Qerr ~ P (%.0f%% credible band)', level * 100)) +
    ggplot2::theme_minimal()
}


#' Plot per-site lines (random intercepts & slopes) with credible intervals
#'
#' @title Per-site partial-pooling lines and credible bands
#' @description
#' Plots **site-specific** regression lines using draws of `intercept_id[j]` and
#' `slope_id[j]` for each site `ID`. Requires more than 1 ID in the initial
#' regression by [BRRAT()]. The output overlays lines and credible bands
#' across sites; optionally facetted into separate panels for readability.
#'
#' @param fit_obj A fitted object returned by [BRRAT()], of class `BRRAT`.
#' @param n_grid Integer; number of x points per site for the prediction grids
#' (default `150`).
#' @param level Credible interval level (default `0.80`) used for site-specific bands.
#' @param facet Logical; if `TRUE`, facet the plot by `ID`; otherwise overlay all sites
#' with colour and fill aesthetics (default `TRUE`).
#'
#' @return A `ggplot2` object showing per-site lines and credible intervals.
#'
#' @details
#' For site \eqn{j}, the line is computed as:
#' \deqn{\mu_j(x) = \mathrm{intercept\_id}_j + \mathrm{slope\_id}_j \, x_{\mathrm{std}}}{
#' mu_j(x) = intercept_id[j] + slope_id[j] * x_std}
#' evaluated on a grid of x values transformed to the centered/scaled space used
#' during fitting. The function reconstructs equal-tail credible bands by propagating
#' posterior draws of `intercept_id[j]` and `slope_id[j]`.
#'
#' @examples
#' \dontrun{
#' res <- BRRAT(Sim,Obs)
#' p_sites <- plot_BRRAT_by_id(res, n_grid = 150, level = 0.80, facet = TRUE)
#' print(p_sites)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_ribbon geom_line labs theme_minimal facet_wrap guides
#' @export
plot_BRRAT_by_id <- function(fit_obj, n_grid = 150, level = 0.80, facet = TRUE) {
  stopifnot(inherits(fit_obj, 'BRRAT'))

  id_levels <- fit_obj$id_levels
  J <- length(id_levels)

  # Extract per-ID generated quantities (iterations x J)
  intercept_draws <- rstan::extract(fit_obj$fit, pars = 'intercept_id')[[1]]
  slope_draws     <- rstan::extract(fit_obj$fit, pars = 'slope_id')[[1]]

  x_min  <- min(fit_obj$data$x, na.rm = TRUE)
  x_max  <- max(fit_obj$data$x, na.rm = TRUE)
  x_std <- seq(x_min, x_max, length.out = n_grid)

  q_lo <- (1 - level) / 2
  q_hi <- 1 - q_lo

  out_list <- vector('list', J)
  for (j in seq_len(J)) {
    mu_draws_j <- sapply(x_std, function(xi) intercept_draws[, j] + slope_draws[, j] * xi)
    mu_mean <- apply(mu_draws_j, 2, mean)
    mu_lo   <- apply(mu_draws_j, 2, stats::quantile, probs = q_lo)
    mu_hi   <- apply(mu_draws_j, 2, stats::quantile, probs = q_hi)
    out_list[[j]] <- data.frame(ID = id_levels[j], P = x_std,
                                mu_mean = mu_mean, mu_lo = mu_lo, mu_hi = mu_hi)
  }
  df_lines <- do.call(rbind, out_list)

  # Downsample points per ID if very large
  data_sc <- data.frame(x=fit_obj$data$x,
                        y=fit_obj$data$y,
                        ID = fit_obj$id_levels[fit_obj$data$group])
  if (nrow(data_sc) > 30000) {
    set.seed(1)
    data_sc <- do.call(rbind, lapply(split(data_sc, data_sc$ID), function(dd) {
      keep <- min(5000, nrow(dd))
      dd[sample(seq_len(nrow(dd)), keep), ]
    }))
  }

  if (!requireNamespace('ggplot2', quietly = TRUE)) stop('ggplot2 is required.')

  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = data_sc, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$ID),
                        alpha = 0.12, size = 0.6) +
    ggplot2::geom_ribbon(data = df_lines, ggplot2::aes(x = .data$P, ymin = .data$mu_lo, ymax = .data$mu_hi, fill = .data$ID),
                         alpha = 0.15) +
    ggplot2::geom_line(data = df_lines, ggplot2::aes(x = .data$P, y = .data$mu_mean, colour = .data$ID),
                       linewidth = 1) +
    ggplot2::labs(x = 'P (normalized)', y = 'Qerr (transformed)',
                  title = sprintf('Per-site partial-pooling lines (%.0f%% credible bands)', level * 100)) +
    ggplot2::theme_minimal()

  if (facet) {
    p <- p + ggplot2::facet_wrap(~ID, scales = 'free_y', ncol = 3) +
      ggplot2::guides(colour = 'none', fill = 'none')
  }
  return(p)
}
