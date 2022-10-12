## ---------------------------
##
## Project name: Data Collection for Strength-endurance Profiles: Can Assessments Be Completed within a Single Session?
##
## Purpose of script: simulation-based analysis of the SV vs. MV protocols; evaluation of RTF._protocol.stan
##
## Author: Benedikt Mitter MSc
##
## Date Created: 2022-Oct-12
##
## Copyright (c) Benedikt Mitter, 2022
## Email: benedikt.mitter@univie.ac.at
##
## ---------------------------


#### Libraries ####
library(clusterGeneration)
library(ggplot2)
library(ggpubr)
library(mvtnorm)
library(truncnorm)


#### Load functions ####
source("./_scripts/functions.R")
set.seed(1234)

#### Stan setup ####
rstan_options(auto_write = TRUE)
options(mc.cores = 4) # number of cores for sampling

#### Define population ####

# define group-level parameters (fixed effects)
mu_a <- 4
mu_b <- 4
mu_c <- 10
mu_d <- 0
mu_e <- 1
mu_f <- 2

sigma <- 1

# define group-level parameter scales; express them in a diagonal matrix
tau = c(1, 1, 2, 0.4, 0.7, 1)
m_diag_tau <- diag(6)
for (i in 1:ncol(m_diag_tau)){
  m_diag_tau[i, i] <- tau[i]
}

# define custom correlation matrix
rho <- c(0.8, 0.8, 0.8, 0, 0, 0, 0.1, 0.1, 0.1, 0, 0.2, 0.2, 0.2, 0, 0.8)
m_R <- matrix(0, 6, 6)
k <- 1
for (i in 2:ncol(m_R)){
  for (j in 1:(i-1)){
    m_R[j, i] <- rho[k]
    k <- k + 1
  }
}
m_R <- m_R + t(m_R) + diag(6)

# alternatively: define random correlation matrix with weak correlations (remove "#" in next line)
# m_R <- rcorrmatrix(6, 100)

# covariance matrix
m_Sigma <- m_diag_tau %*% m_R %*% m_diag_tau


#### Build data structure ####
nSubjects <- 14

# sample random effects
rand_eff <- mvrnorm(nSubjects, rep(0, 6), m_Sigma)

# protocol structure
protocol <- array(c(matrix(0, nrow = nSubjects, ncol = 3),
                    matrix(1, nrow = nSubjects, ncol = 3)),
                  dim = c(nSubjects, 3, 2)
)

# dummy variables for trial B (80% 1-RM) and trial C (70% 1-RM)
trialB <- t(matrix(c(0, 1, 0), 6, nSubjects))
trialB <- array(trialB, dim = c(nSubjects, 3, 2))

trialC <- t(matrix(c(0, 0, 1), 6, nSubjects))
trialC <- array(trialC, dim = c(nSubjects, 3, 2))

# simulate reps from model
reps <- array(0, dim = c(nSubjects, 3, 2))

for (s in 1:nSubjects){
  for (t in 1:3){
    for (p in 1:2){
      
      mu <- (mu_a + rand_eff[s, 1]) + 
        (mu_b + rand_eff[s, 2]) * trialB[s, t, p] + 
           (mu_c + rand_eff[s, 3]) * trialC[s, t, p] + 
        protocol[s, t, p] * (
          (mu_d + rand_eff[s, 4]) + 
            (mu_e + rand_eff[s, 5]) * trialB[s, t, p] + 
            (mu_f + rand_eff[s, 6]) * trialC[s, t, p]
          )
      reps[s, t, p] <- rtruncnorm(1, 0, Inf, mu, sigma)
      
    }
  }
}


# data list for Stan

scale1 <- 10
scale2 <- 10

dataList <- list(
  nSubjects = nSubjects,
  protocol = protocol,
  reps = floor(reps),
  dummy1 = trialB,
  dummy2 = trialC,
  
  prior_a = c(0,scale1),
  prior_b = c(0,scale1),
  prior_c = c(0,scale1),
  prior_d = c(0,scale1),
  prior_e = c(0,scale1),
  prior_f = c(0,scale1),
  
  prior_sigma = c(0,scale2),
  prior_tau1 = c(0,scale2),
  prior_tau2 = c(0,scale2),
  prior_tau3 = c(0,scale2),
  prior_tau4 = c(0,scale2),
  prior_tau5 = c(0,scale2),
  prior_tau6 = c(0,scale2),
  
  prior_eta = 1
)

#### Evaluate simulation ####

# run model
fit <- run.model("_scripts/Stan/RTF_protocol.stan", dataList, adapt_delta = 0.95)

# model diagnostics
cat("Diagnostics for fitted model: \n", diagnose.model(fit))
pars <- c("mu_a", "mu_b", "mu_c", "mu_d", "mu_e", "mu_f")
rstan::traceplot(fit, par = pars, nrow = ceiling(length(pars)/2), inc_warmup = TRUE)

# evaluate posteriors
plot_list <- list()
for (i in 1:length(pars)){
  
  mu_tmp <- eval(parse(text = pars[i]))
  plot_tmp <- stan_plot(fit, pars = pars[i], show_density=T, fill_color = 'skyblue', ci_level = 0.95, outer_level = 0.99) + 
    geom_vline(xintercept = mu_tmp, color ="red", size = 1.5) +
    theme_bw() +
    theme(
          panel.grid = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 12)
          )
  assign(paste0("plot_", pars[i]), plot_tmp)
  plot_list[[i]] <- eval(parse(text = paste0("plot_", pars[i])))
  
}

plots_merged <- ggarrange(plotlist = plot_list, 
                          labels = pars,
                          label.x = 0,
                          label.y = 0.95,
                          ncol = 2, nrow = 3)

plots_merged <- annotate_figure(plots_merged, 
                top = text_grob(paste0("Population parameters used for simulations are shown as red lines \n", 
                "Blue areas mark the 2.5th and 97.5th percentile"), 
                face = "bold", size = 14))

print(plots_merged)
