
// hierarchical linear regression with random intercepts and slopes by ID
// Likelihood: y ~ Normal(mu, sigma)
// mu = alpha + beta * x + b0[group] + b1[group] * x
// Non-centered parameterization with correlated random effects

data {
  int<lower=1> N;                        // number of observations
  vector[N] x;                           // predictor (P)
  vector[N] y;                           // response (Qerr)
  int<lower=1> J;                        // number of groups (sites)
  int<lower=1, upper=J> group[N];        // group index for each observation
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;

  vector<lower=0>[2] tau;
  cholesky_factor_corr[2] L_Omega;
  matrix[2, J] z;
}
transformed parameters {
  matrix[2, J] b;
  b = diag_pre_multiply(tau, L_Omega) * z;
}
model {
  alpha ~ normal(0, 5);
  beta  ~ normal(0, 2);
  sigma ~ exponential(0.1);

  tau ~ exponential(1);
  L_Omega ~ lkj_corr_cholesky(2);
  to_vector(z) ~ normal(0, 1);

  for (n in 1:N) {
    real mu_n = alpha + beta * x[n] + b[1, group[n]] + b[2, group[n]] * x[n];
    y[n] ~ normal(mu_n, sigma);
  }
}
generated quantities {
  corr_matrix[2] Omega;
  vector[J] intercept_id;
  vector[J] slope_id;
  Omega = multiply_lower_tri_self_transpose(L_Omega);
  for (j in 1:J) {
    intercept_id[j] = alpha + b[1, j];
    slope_id[j]     = beta  + b[2, j];
  }
}
