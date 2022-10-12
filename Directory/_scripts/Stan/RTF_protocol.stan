data {
  
  int<lower=0> nSubjects;
  int<lower=0> protocol[nSubjects, 3, 2];
  int<lower=0> dummy1[nSubjects, 3, 2];
  int<lower=0> dummy2[nSubjects, 3, 2];
  real<lower=0> reps[nSubjects, 3, 2];
  
  vector[2] prior_a;
  vector[2] prior_b;
  vector[2] prior_c;
  vector[2] prior_d;
  vector[2] prior_e;
  vector[2] prior_f;
  
  
  vector<lower=0>[2] prior_sigma;
  vector<lower=0>[2] prior_tau1;
  vector<lower=0>[2] prior_tau2;
  vector<lower=0>[2] prior_tau3;
  vector<lower=0>[2] prior_tau4;
  vector<lower=0>[2] prior_tau5;
  vector<lower=0>[2] prior_tau6;
  
  real<lower=0> prior_eta;
  
}


parameters {
  
  // group-level parameters
  real<lower=0> sigma;
  
  real mu_a;
  real mu_b;
  real mu_c;
  real mu_d;
  real mu_e;
  real mu_f;
  
  // parameters for Cholesky factorization
  matrix[6, nSubjects] z_u;
  cholesky_factor_corr[6] L_u;
  vector<lower=0>[6] tau_u;
  
}


transformed parameters {
  // matrix of correlated parameters (Cholesky factorization)
  matrix[nSubjects, 6] u = (diag_pre_multiply(tau_u, L_u) * z_u)';
  
  // calculating subject level parameters; 
  vector[nSubjects] RE_a = mu_a + u[,1];
  vector[nSubjects] RE_b = mu_b + u[,2];
  vector[nSubjects] RE_c = mu_c + u[,3];
  vector[nSubjects] RE_d = mu_d + u[,4];
  vector[nSubjects] RE_e = mu_e + u[,5];
  vector[nSubjects] RE_f = mu_f + u[,6];
  
} 


model {
  // group-level prior
  target += cauchy_lpdf(sigma | prior_sigma[1], prior_sigma[2]);

  target += cauchy_lpdf(mu_a | prior_a[1], prior_a[2]); 
  target += cauchy_lpdf(mu_b | prior_b[1], prior_b[2]);
  target += cauchy_lpdf(mu_c | prior_c[1], prior_c[2]); 
  target += cauchy_lpdf(mu_d | prior_d[1], prior_d[2]);
  target += cauchy_lpdf(mu_e | prior_e[1], prior_e[2]);
  target += cauchy_lpdf(mu_f | prior_f[1], prior_f[2]);

  // generate uncorrelated vectors for Cholesky facorization
  target += std_normal_lpdf(to_vector(z_u));
  
  // prior for Cholesky facorization
  target += lkj_corr_cholesky_lpdf(L_u | prior_eta);
  
  // priors for scaling of correlated parameters
  target += cauchy_lpdf(tau_u[1] | prior_tau1[1], prior_tau1[2]); 
  target += cauchy_lpdf(tau_u[2] | prior_tau2[1], prior_tau2[2]); 
  target += cauchy_lpdf(tau_u[3] | prior_tau3[1], prior_tau3[2]); 
  target += cauchy_lpdf(tau_u[4] | prior_tau4[1], prior_tau4[2]); 
  target += cauchy_lpdf(tau_u[5] | prior_tau5[1], prior_tau5[2]);
  target += cauchy_lpdf(tau_u[6] | prior_tau6[1], prior_tau6[2]);



  // likelihood
  for (s in 1:nSubjects) {
    
    for (t in 1:3){
      
      for (p in 1:2){
        
        target += normal_lpdf(reps[s, t, p] | RE_a[s] + RE_b[s] * dummy1[s, t, p] + RE_c[s] * dummy2[s, t, p] + protocol[s, t, p] * (RE_d[s] + RE_e[s] * dummy1[s, t, p] + RE_f[s] * dummy2[s, t, p]), sigma);
        
      }
      
    }
    
  }
  
}


generated quantities {
  
  vector[nSubjects] log_lik;
  
  real reps_pred[nSubjects,3, 2];
  real Residuals[nSubjects,3, 2];
  
  matrix[6,6] tau_m;
  cov_matrix[6] Sigma_cov;
  vector[6] sim;
  
  real mu_SV_90 = mu_a;
  real mu_SV_80 = mu_a + mu_b;
  real mu_SV_70 = mu_a + mu_c;
  
  real mu_delta_90 = mu_d;
  real mu_delta_80 = mu_d + mu_e;
  real mu_delta_70 = mu_d + mu_f;
  
  vector[nSubjects] SV_90 = RE_a;
  vector[nSubjects] SV_80 = RE_a + RE_b;
  vector[nSubjects] SV_70 = RE_a + RE_c;
  
  vector[nSubjects] delta_90 = RE_d;
  vector[nSubjects] delta_80 = RE_d + RE_e;
  vector[nSubjects] delta_70 = RE_d + RE_f;
  
  corr_matrix[6] Sigma_corr = multiply_lower_tri_self_transpose(L_u);
  
  tau_m = diag_matrix(tau_u);
  Sigma_cov = tau_m * L_u * L_u' * tau_m;
  sim = multi_normal_rng([mu_a, mu_b, mu_c, mu_d, mu_e, mu_f], Sigma_cov);
  
  for (s in 1:nSubjects) {
    
    log_lik[s] = 0;
    
    for (t in 1:3) {
      
      for (p in 1:2){
        
        log_lik[s] += normal_lpdf(reps[s, t, p] | RE_a[s] + RE_b[s] * dummy1[s, t, p] + RE_c[s] * dummy2[s, t, p] + protocol[s, t, p] * (RE_d[s] + RE_e[s] * dummy1[s, t, p] + RE_f[s] * dummy2[s, t, p]), sigma);
        reps_pred[s, t, p] = RE_a[s] + RE_b[s] * dummy1[s, t, p] + RE_c[s] * dummy2[s, t, p] + protocol[s, t, p] * (RE_d[s] + RE_e[s] * dummy1[s, t, p] + RE_f[s] * dummy2[s, t, p]);
        Residuals[s, t, p] = reps[s, t, p] - reps_pred[s, t, p];
        
      }
      
    }
    
  }
  
}
