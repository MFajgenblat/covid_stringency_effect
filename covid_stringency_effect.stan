functions {
  
  // One-dimensional Gaussian process
  vector GP_1D(real[] x, real rho, real alpha, vector eta) {
    return(cholesky_decompose(add_diag(gp_exp_quad_cov(x, alpha, rho), 1e-9)) * eta);
  }

}

data {
  
  int<lower=0> N;
  int<lower=0> N_countries;
  int<lower=0> N_cov;
  int<lower=0> N_bf;
  int<lower=0> N_lags;
  int country[N];
  matrix[N,N_lags] X[N_cov];
  matrix[N,N_lags] bf[N_bf];
  real lag_range[N_lags];
  vector[N] y;
  int<lower=0> N_pred;
  matrix[N_pred,N_bf] bf_pred;
  int countrystart[N_countries];
  int countryend[N_countries];
  int N_pred_country;
  matrix[N_pred_country,N_cov] X_pred_country;
  matrix[N_pred_country,N_bf] bf_pred_lowres;
  int<lower=0,upper=1> t;
  int<lower=2> nu_t;
  real prior_intercept[2];
  real prior_tau[2];
  real prior_Omega;
  real prior_nu_mvt[2];
  real prior_alpha_smooth[2];
  real prior_rho_lags[2];
  real prior_alpha_lags[2];
  real prior_sigma[2];
  
}

parameters {
  
  // Model intercept
  real intercept;
  
  // Multivariate t random country effect parameters
  vector[N_cov+1] beta_country[N_countries];
  corr_matrix[N_cov+1] Omega;
  vector<lower=0>[N_cov+1] tau;
  real<lower=3> nu_mvt;
  
  // Tensor-product spline smooth parameters
  real<lower=0> alpha_smooth;
  vector[N_bf] eta_smooth;
  
  // Delay distribution GP parameters
  real<lower=0> rho_lags;
  real<lower=0> alpha_lags;
  vector[N_lags] eta_lags;
  
  // Residual variation scale parameter
  real<lower=0> sigma[N_countries];
  
}


transformed parameters {
 
  // Multivariate t random country effect covariance matrix
  matrix[N_cov+1,N_cov+1] beta_sigma = quad_form_diag(Omega, tau);
  // Delay distribution based on a softmax-transformed Gaussian process
  vector[N_lags] delay_distribution = softmax(GP_1D(lag_range, rho_lags, alpha_lags, eta_lags));
 
}

model {
  
  // Declaring the temporally weighted predictor matrices
  matrix[N,N_cov+1] X_2D;
  matrix[N,N_bf] bf_2D;
  
  // Prior on the model intercept
  target += normal_lpdf(intercept | prior_intercept[1], prior_intercept[2]);
  
  // Priors on the multivariate t random country effect parameters
  target += normal_lpdf(tau | prior_tau[1], prior_tau[2]);
  target += lkj_corr_lpdf(Omega | prior_Omega);
  target += gamma_lpdf(nu_mvt | prior_nu_mvt[1], prior_nu_mvt[2]);
  
  // Priors on the tensor-product spline smooth parameters
  target += std_normal_lpdf(eta_smooth);
  target += normal_lpdf(alpha_smooth | prior_alpha_smooth[1], prior_alpha_smooth[2]);
  
  // Priors on the delay distribution GP parameters
  target += inv_gamma_lpdf(rho_lags | prior_rho_lags[1], prior_rho_lags[2]);
  target += normal_lpdf(alpha_lags | prior_alpha_lags[1], prior_alpha_lags[2]);
  target += std_normal_lpdf(eta_lags);
  
  // Prior on the residual variation scale parameter
  for (i in 1:N_countries) {target += normal_lpdf(sigma[i] | prior_sigma[1], prior_sigma[2]);}
  
  // Temporally weigh the predictors
  X_2D[,1] = rep_vector(1.0, N);
  for (i in 1:N_cov) {
    X_2D[,i+1] = X[i] * delay_distribution;
  }
  for (i in 1:N_bf) {
    bf_2D[,i] = bf[i] * delay_distribution;
  }
  
  for (i in 1:N_countries) {
    // Multivariate t random country effects
    target += multi_student_t_lpdf(beta_country[i] | nu_mvt, rep_vector(0, N_cov+1), beta_sigma);
    // Model likelihood
    if (t == 0) {
      // Normal likelihood
      target += normal_id_glm_lpdf(y[countrystart[i]:countryend[i]] | append_col(X_2D[countrystart[i]:countryend[i],], bf_2D[countrystart[i]:countryend[i],]), intercept, append_row(beta_country[i], eta_smooth * alpha_smooth), sigma[i]);
    } else if (t == 1) {
      // Student's t likelihood
      target += student_t_lpdf(y[countrystart[i]:countryend[i]] | nu_t, intercept + append_col(X_2D[countrystart[i]:countryend[i],], bf_2D[countrystart[i]:countryend[i],]) * append_row(beta_country[i], eta_smooth * alpha_smooth), sigma[i]);
    }
  }
  
 }

generated quantities {
  
  vector[N_pred] preds;
  vector[N_pred_country] preds_country[N_countries];
  matrix[N,N_cov+1] X_2D;
  matrix[N,N_bf] bf_2D;
  real y_rep[N];
  
  // Average predicted change for a grid of predictor values
  preds = intercept + bf_pred * (eta_smooth * alpha_smooth);
  
  // Temporally weigh the predictors
  X_2D[,1] = rep_vector(1.0, N);
  for (i in 1:N_cov) {
    X_2D[,i+1] = X[i] * delay_distribution;
  }
  for (i in 1:N_bf) {
    bf_2D[,i] = bf[i] * delay_distribution;
  }
  
  for (i in 1:N_countries) {
    // Posterior predictive distribution
    if (t == 0) {
      // Normal likelihood
      y_rep[countrystart[i]:countryend[i]] = normal_rng(intercept + append_col(X_2D[countrystart[i]:countryend[i],], bf_2D[countrystart[i]:countryend[i],]) * append_row(beta_country[i], eta_smooth * alpha_smooth), sigma[i]);
      // Student's t likelihood
      y_rep[countrystart[i]:countryend[i]] = student_t_rng(nu_t, intercept + append_col(X_2D[countrystart[i]:countryend[i],], bf_2D[countrystart[i]:countryend[i],]) * append_row(beta_country[i], eta_smooth * alpha_smooth), sigma[i]);
    }
    // Country-specific predicted change for a grid of predictor values
    preds_country[i] = intercept + append_col(rep_vector(1, N_pred_country), X_pred_country) * beta_country[i] + bf_pred_lowres * (eta_smooth * alpha_smooth);
  }
  
}
