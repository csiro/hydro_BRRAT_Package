
set.seed(123)
options(mc.cores = 1L)

# Let tests override Stan sampling args if your BRRAT() forwards "..." to rstan::sampling()
options(BRRAT.iter   = 400L)
options(BRRAT.warmup = 200L)
options(BRRAT.chains = 1L)

iters   <- getOption("BRRAT.iter",   400L)
warmup  <- getOption("BRRAT.warmup", 200L)
chains  <- getOption("BRRAT.chains", 1L)
