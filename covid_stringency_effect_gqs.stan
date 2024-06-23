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

generated quantities {
  
  # 2D predictions at high resolution
  vector[N_pred] preds;
  # 2D country-specific predictions at lower resolution
  vector[N_pred_country] preds_country[N_countries];
  # Posterior predictive distribution
  real y_rep[N];
  
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
  
  // Average predicted change for a grid of predictor values
  preds = intercept + bf_pred * (eta_smooth * alpha_smooth);
  
  for (i in 1:N_countries) {
    // Posterior predictive distribution
    if (t == 0) {
      // Normal likelihood
      y_rep[countrystart[i]:countryend[i]] = normal_rng(intercept + temporal_effect[day[countrystart[i]:countryend[i]],i] + append_col(X_2D[countrystart[i]:countryend[i],], bf_2D[countrystart[i]:countryend[i],]) * append_row(beta_country[i], eta_smooth * alpha_smooth), sigma[i]);
    } else if (t == 1) {
      // Student's t likelihood
      y_rep[countrystart[i]:countryend[i]] = student_t_rng(nu_t, intercept + temporal_effect[day[countrystart[i]:countryend[i]],i] + append_col(X_2D[countrystart[i]:countryend[i],], bf_2D[countrystart[i]:countryend[i],]) * append_row(beta_country[i], eta_smooth * alpha_smooth), sigma[i]);
    }
    // Country-specific predicted change for a grid of predictor values
    preds_country[i] = intercept + append_col(rep_vector(1, N_pred_country), X_pred_country) * beta_country[i] + bf_pred_lowres * (eta_smooth * alpha_smooth);
  }
  
}
