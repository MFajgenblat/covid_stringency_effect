functions {

  // One-dimensional Gaussian process with exponentiated quadratic kernel
  vector GP_1D_exp_quad(real[] x, real rho, real alpha, vector eta) {
    return(cholesky_decompose(add_diag(gp_exp_quad_cov(x, alpha, rho), 1e-9)) * eta);
  }

}

data {
  
  int<lower=0> N;
  int<lower=0> N_countries;
  int<lower=0> N_cov;
  int<lower=0> N_bf;
  int<lower=0> N_lags;
  int<lower=0> N_days;
  int<lower=0> N_bf_days;
  int country[N];
  matrix[N,N_lags] X[N_cov];
  matrix[N,N_lags] bf[N_bf];
  real lag_range[N_lags];
  vector[N] y;
  int<lower=0,upper=1> t;
  int<lower=2> nu_t;
  int<lower=0,upper=1> temporal;
  matrix[N_days,N_bf_days] bf_days;
  real day_range[N_bf_days];
  int<lower=1,upper=N_days> day[N];
  real prior_intercept[2];
  real prior_tau[2];
  real prior_Omega;
  real prior_nu_mvt[2];
  real prior_alpha_smooth[2];
  real prior_rho_lags[2];
  real prior_alpha_lags[2];
  real prior_rho_temporal[2];
  real prior_alpha_temporal[2];
  real prior_sigma[2];
  int<lower=0> N_pred;
  matrix[N_pred,N_bf] bf_pred;
  int countrystart[N_countries];
  int countryend[N_countries];
  int N_pred_country;
  matrix[N_pred_country,N_cov] X_pred_country;
  matrix[N_pred_country,N_bf] bf_pred_lowres;

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
  
  // Temporal (spline-projected) Gaussian process parameters
  real<lower=0> rho_temporal[temporal,N_countries];
  real<lower=0> alpha_temporal[temporal,N_countries];
  vector[N_bf_days] eta_temporal[temporal,N_countries];
  
  // Residual variation scale parameter
  real<lower=0> sigma[N_countries];
  
}

transformed parameters {
 
  // Multivariate t random country effect covariance matrix
  matrix[N_cov+1,N_cov+1] beta_sigma = quad_form_diag(Omega, tau);
  
  // Delay distribution based on a softmax-transformed Gaussian process
  vector[N_lags] delay_distribution = softmax(GP_1D_exp_quad(lag_range, rho_lags, alpha_lags, eta_lags));
  
  // Temporal (spline-projected) Gaussian process
  matrix[N_bf_days,N_countries] GP_temporal_effect[temporal];
  if (temporal == 1) {
    for (i in 1:N_countries) {
      GP_temporal_effect[temporal,,i] = GP_1D_exp_quad(day_range, rho_temporal[temporal,i], alpha_temporal[temporal,i], eta_temporal[temporal,i]);
    }
  }
  
}

model {
  
  // Declaring the temporally weighted predictor matrices
  matrix[N_days,N_countries] temporal_effect;
  
  // Temporally weigh the predictors
  matrix[N,N_cov+1] X_2D;
  matrix[N,N_bf] bf_2D;
  
  // Computing the country-specific temporal effects
  for (i in 1:N_countries) {
    if (temporal == 0) {
      temporal_effect[,i] = rep_vector(0.0, N_days);
    } else if (temporal == 1) {
      temporal_effect[,i] = bf_days * GP_temporal_effect[temporal,,i];
    }
  }
  
  // Temporally weigh the predictors
  X_2D[,1] = rep_vector(1.0, N);
  for (i in 1:N_cov) {
    X_2D[,i+1] = X[i] * delay_distribution;
  }
  for (i in 1:N_bf) {
    bf_2D[,i] = bf[i] * delay_distribution;
  }
  
  // Prior on the model intercept
  target += normal_lpdf(intercept | prior_intercept[1], prior_intercept[2]);
  
  // Priors on the multivariate t random country effect parameters
  target += normal_lpdf(tau | prior_tau[1], prior_tau[2]);
  target += lkj_corr_lpdf(Omega | prior_Omega);
  target += gamma_lpdf(nu_mvt | prior_nu_mvt[1], prior_nu_mvt[2]);
  
  // Priors on the tensor-product spline smooth parameters
  target += std_normal_lpdf(eta_smooth);
  target += normal_lpdf(alpha_smooth | prior_alpha_smooth[1], prior_alpha_smooth[2]);
  
  // Priors on the delay distribution hyperparameters
  target += inv_gamma_lpdf(rho_lags | prior_rho_lags[1], prior_rho_lags[2]);
  target += normal_lpdf(alpha_lags | prior_alpha_lags[1], prior_alpha_lags[2]);
  target += std_normal_lpdf(eta_lags);
  
  // Priors on the temporal (spline-projected) Gaussian process hyperparameters
  if (temporal == 1) {
    for (i in 1:N_countries) {
      target += inv_gamma_lpdf(rho_temporal[temporal,i] | prior_rho_temporal[1], prior_rho_temporal[2]);
      target += normal_lpdf(alpha_temporal[temporal,i] | prior_alpha_temporal[1], prior_rho_temporal[2]);
      target += std_normal_lpdf(eta_temporal[temporal,i]);
    }
  }
  
  // Prior on the residual variation scale parameter
  for (i in 1:N_countries) {target += normal_lpdf(sigma[i] | prior_sigma[1], prior_sigma[2]);}
  
  for (i in 1:N_countries) {
    // Multivariate t random country effects
    target += multi_student_t_lpdf(beta_country[i] | nu_mvt, rep_vector(0, N_cov+1), beta_sigma);
    // Model likelihood
    if (t == 0) {
      // Normal likelihood
      target += normal_lpdf(y[countrystart[i]:countryend[i]] | intercept + temporal_effect[day[countrystart[i]:countryend[i]],i] + append_col(X_2D[countrystart[i]:countryend[i],], bf_2D[countrystart[i]:countryend[i],]) * append_row(beta_country[i], eta_smooth * alpha_smooth), sigma[i]);
    } else if (t == 1) {
      // Student's t likelihood
      target += student_t_lpdf(y[countrystart[i]:countryend[i]] | nu_t, intercept + temporal_effect[day[countrystart[i]:countryend[i]],i] + append_col(X_2D[countrystart[i]:countryend[i],], bf_2D[countrystart[i]:countryend[i],]) * append_row(beta_country[i], eta_smooth * alpha_smooth), sigma[i]);
    }
  }
  
 }
