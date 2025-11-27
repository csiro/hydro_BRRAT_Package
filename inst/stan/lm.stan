//inst/stan/lm.stan
data {
  int<lower=1> N;
  vector[N] x;
  vector[N] y;
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}
model {
  alpha ~ normal(0, 5);
  beta  ~ normal(0, 2);
  sigma ~ exponential(0.1);

  for (n in 1:N) {
    real mu_n = alpha + beta * x[n];
    y[n] ~ normal(mu_n, sigma);
  }
}
