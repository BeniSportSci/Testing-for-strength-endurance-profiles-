data {
  int<lower=1> nSubjects;
  matrix<lower=0>[nSubjects, 4] load;
  matrix<lower=0>[nSubjects, 4] reps;
  
  vector[2] prior_a;
  vector[2] prior_b;
  
  vector<lower=0>[2] prior_sigma;
  vector<lower=0>[2] prior_tau1;
  vector<lower=0>[2] prior_tau2;
  
  real<lower=0> prior_eta;
}

transformed data {
  matrix<lower=0>[nSubjects, 4] load_reciprocal = 1 ./ load;
  matrix[nSubjects, 4] load_std = (load_reciprocal-mean(load_reciprocal))/sd(load_reciprocal);
  matrix[nSubjects, 4] reps_std = (reps-mean(reps))/sd(reps);
}

parameters {
  // group-level parameters
  real<lower=0> sigma;
  
  real alpha;
  real beta;
  
  // cholesky factorization
  matrix[2, nSubjects] z_u;
  cholesky_factor_corr[2] L_u;
  vector<lower=0>[2] tau_u;
  
}


transformed parameters {
  // matrix of correlated parameters (Cholesky factorization)
  matrix[nSubjects, 2] u = (diag_pre_multiply(tau_u, L_u) * z_u)';
  
  // calculating absolute subject level parameters; 
  vector[nSubjects] a = alpha + u[,1];
  vector[nSubjects] b = beta + u[,2];
  
} 


model {
  // group-level prior
  target += normal_lpdf(sigma | prior_sigma[1], prior_sigma[2]);

  target += normal_lpdf(alpha | prior_a[1], prior_a[2]); 
  target += normal_lpdf(beta | prior_b[1], prior_b[2]);
  
  
  // generate uncorrelated vectors for Cholesky factorization
  target += std_normal_lpdf(to_vector(z_u));
  
  // prior for Cholesky factorization
  target += lkj_corr_cholesky_lpdf(L_u | prior_eta);
  
  // priors for scaling of correlated parameters
  target += normal_lpdf(tau_u[1] | prior_tau1[1], prior_tau1[2]); 
  target += normal_lpdf(tau_u[2] | prior_tau2[1], prior_tau2[2]);
  
  
  // likelihood
  for (s in 1:nSubjects) {
    
    for (t in 1:4) {
        
        target += normal_lpdf(load_std[s, t] | a[s] + b[s] * reps_std[s, t], sigma);
        
    }
    
  }
  
}


generated quantities {
  
  real log_lik[nSubjects];
  
  real alpha_rec;
  real beta_rec;
  real sigma_rec;
  
  real a_rec[nSubjects];
  real b_rec[nSubjects];
  
  matrix[nSubjects,4] load_pred;
  real Residuals[nSubjects,4];
  real<lower=0> SEE;
  
  corr_matrix[2] Sigma_corr;
  Sigma_corr = multiply_lower_tri_self_transpose(L_u);
  
  alpha_rec = alpha * sd(load_reciprocal) + mean(load_reciprocal) - beta * sd(load_reciprocal) * mean(reps) / sd(reps);
  beta_rec = beta * sd(load_reciprocal) / sd(reps);
  sigma_rec = sigma * sd(load_reciprocal);
  
  for (s in 1:nSubjects) {
    
    log_lik[s] = 0;
    
    a_rec[s] = a[s] * sd(load_reciprocal) + mean(load_reciprocal) - b[s] * sd(load_reciprocal) * mean(reps) / sd(reps);
    b_rec[s] = b[s] * sd(load_reciprocal) / sd(reps);
    
    for (t in 1:4) {
        
        // log likelihood
        log_lik[s] += normal_lpdf(load_std[s, t] | a[s] + b[s] * reps_std[s, t], sigma);
        
        //predicted load
        load_pred[s, t] = 1 / (a_rec[s] + b_rec[s] * reps[s, t]);
        Residuals[s, t] = load[s, t] - load_pred[s, t];
    }
    
  }
  
  SEE = sqrt(sum(square(load_pred - load)) / (nSubjects * 4 - 2));
  
}
