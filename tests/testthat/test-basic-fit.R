
set.seed(10)
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

test_that("BRRAT returns expected object", {
  expect_warning(
  res <- BRRAT(Sim, Obs, iter = iters, warmup = warmup, chains = chains),
  regexp = "Samples Size")

  # Class and core members
  expect_s3_class(res, "BRRAT")
  expect_true(!is.null(res$fit), info = "Stan fit object should be present")

  # Posterior slope draws should be accessible either via $beta or rstan extract
  b <- NULL
  if (!is.null(res$beta)) {
    b <- res$beta
    expect_true(is.numeric(b))
  } else {
    # fall back to extracting from Stan fit
    pars <- try(rstan::extract(res$fit, pars = "beta"), silent = TRUE)
    expect_false(inherits(pars, "try-error"))
    b <- pars$beta
    expect_true(is.numeric(b))
  }

  expect_true(length(b) > 50, info = "Should have many posterior draws")
  expect_true(all(is.finite(b)))
})
