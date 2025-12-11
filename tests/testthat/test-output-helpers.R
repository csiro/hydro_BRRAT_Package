
test_that("Population plot returns a ggplot object", {
  set.seed(30)
  samples <- 50
  #' #generate some annual rainfall totals
  Clim <- sort(pmax(0,rnorm(samples,600,100)))
  #calculate annual runoff from tanh curve
  observed <- pmax(0,(Clim - 150) - 450 * tanh((Clim - 150)/450))
  Obs <- data.frame(P=Clim,Q=observed,year=1:length(Clim),ID=1)
  #'
  #' #non-robust model, underestimate the wetter years
  Sim <- observed*seq(1,0.7,length.out=length(observed))*rnorm(samples,1,0.05)
  Sim <- data.frame(Q=Sim,year=1:length(Clim),ID=1)

  expect_warning(
  fit <- BRRAT(Sim, Obs, iter = iters, warmup = warmup, chains = chains),
  regexp = "Samples Size")

  p   <- plot_BRRAT_population(fit, n_grid = 50, level = 0.8)
  expect_s3_class(p, "ggplot")
})

test_that("By-ID plot and statistics are produced when model has per-site parameters", {

  set.seed(30)
  samples <- 50
  #' #generate some annual rainfall totals
  Clim <- sort(pmax(0,rnorm(samples,600,100)))
  #calculate annual runoff from tanh curve
  observed <- pmax(0,(Clim - 150) - 450 * tanh((Clim - 150)/450))
  Obs <- data.frame(P=c(Clim,Clim),
                    Q=c(observed,observed*0.8),
                    year=c(1:length(Clim),1:length(Clim)),
                    ID=c(rep(1,length(Clim)),rep(2,length(Clim))))

  #' #non-robust model, underestimate the wetter years
  Sim <- Obs[,c("year","Q","ID")]
  Sim$Q <- Sim$Q*rnorm(samples*2,1,0.05) #add some noise, but no trend

  expect_warning(
    fit <- BRRAT(Sim, Obs, iter = iters, warmup = warmup, chains = chains),
    regexp = "Samples Size")

  p   <- plot_BRRAT_population(fit, n_grid = 50, level = 0.8)
  expect_s3_class(p, "ggplot")

  p <- plot_BRRAT_by_id(fit, n_grid = 50, level = 0.8, facet = TRUE)
  expect_s3_class(p, "ggplot")

  out <- summarise_BRRAT(fit, level = 0.95)

  expect_type(out, "list")
  expect_true(all(c("global", "per_id") %in% names(out)))

  # Global table should have slope estimates (median and interval)
  expect_true(all(c("param", "mean", "sd",  "q_lo", "q_hi", "prob_gt0", "prob_lt0","p_nonzero") %in% names(out$global)),
              info = "Global summary must include slope and probability statistics")

  # test per_id outputs exist
  expect_true(all(c("ID", "intercept_mean", "intercept_sd", "intercept_q_lo",
                      "intercept_q_hi", "slope_mean", "slope_sd", "slope_q_lo",
                      "slope_q_hi" ) %in% names(out$per_id)),info = "per_id summary should include intercept and slope information")
  expect_true(nrow(out$per_id)==2, info = "there should be two rows in out$per_id for the two IDs tested")

})
