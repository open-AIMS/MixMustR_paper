// Code built by MixMustR on 2025-06-21 19:31:36.178196
// User options:
  // sample_tracer = FALSE
  // fix_unsampled = TRUE
  // hierarchical = FALSE
data {
  int<lower=1> N; // number of observations
  int<lower=1> J; // number of sampled (a.k.a considered) sources
  int<lower=1> K; // number of tracers
  matrix[N, K] Y; // data stream 1
  matrix[N, J + 1] ln_rho; // Log mixing proportions, data stream 2
  matrix[N, J + 1] sigma_ln_rho; // Confidence around `ln_rho`
  matrix[J, K] x; // sample means for each source and tracer
}
parameters {
  vector[J + 1] zeta[N]; // zeta from softmax
  real<lower=0> sigma[K]; // residual SDs for each tracer
}
transformed parameters {
  simplex[J + 1] p[N]; // probability vector for each site
  for (n in 1:N) {
    p[n] = softmax(zeta[n]);
  }
}
model {
  // Mixture model
  for (n in 1:N) {
    for (k in 1:K) {
      real phi[J + 1];
      // Contribution from sampled sources
      for (j in 1:J) {
        phi[j] = log(p[n][j]) + normal_lpdf(Y[n, k] | x[j, k], sigma[k]);
      }
      // Contribution from unsampled sources
      real mean_x_k = mean(x[, k]);
      phi[J + 1] = log(p[n][J + 1]) + normal_lpdf(Y[n, k] | mean_x_k, sigma[k]);
      // Log-sum-exp for numerical stability
      target += log_sum_exp(phi);
    }
  }
  // Independent evaluation of mixing proportion (zeta scale)
  for (n in 1:N) {
    zeta[n] ~ normal(ln_rho[n, ], sigma_ln_rho[n, ]);
  }
}
