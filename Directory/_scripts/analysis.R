## ---------------------------
##
## Project name: Data Collection for Strength-endurance Profiles: Can Assessments Be Completed within a Single Session?
##
## Purpose of script: analysis of SV vs. MV protocols
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
library(bayestestR)
library(fitdistrplus)
library(greekLetters)
library(loo)
library(matrixStats)
library(posterior)
library(readxl)
library(rstan)
library(writexl)


#### Load functions ####
source("_scripts/functions.R")


#### Stan Setup ####
rstan_options(auto_write = TRUE)
options(mc.cores = 4) # number of cores for sampling


#### Prepare dataList ####
data_directory <- ('./_data/rawdata.xlsx')
rawdata <- as.data.frame(read_xlsx(data_directory, col_names = T))

nSubjects <- nrow(rawdata)

protocol <- array(c(matrix(0, nrow = nSubjects, ncol = 3),
                    matrix(1, nrow = nSubjects, ncol = 3)),
                  dim = c(nSubjects, 3, 2)
                  )

reps <- rawdata[, (ncol(rawdata)-5):ncol(rawdata)]
colnames(reps) <- NULL
reps <- as.matrix(reps)
reps <- array(reps, dim = c(nSubjects, 3, 2))

trialA <- t(matrix(c(1, 0, 0), 6, nSubjects))
trialA <- array(trialA, dim = c(nSubjects, 3, 2))

trialB <- t(matrix(c(0, 1, 0), 6, nSubjects))
trialB <- array(trialB, dim = c(nSubjects, 3, 2))

trialC <- t(matrix(c(0, 0, 1), 6, nSubjects))
trialC <- array(trialC, dim = c(nSubjects, 3, 2))

# data list for Stan
scale <- 5 # choose weakly informative scale from the prior sensitivity analysis

dataList <- list(
  nSubjects = nSubjects,
  protocol = protocol,
  reps = reps,
  dummy0 = trialA,
  dummy1 = trialB,
  dummy2 = trialC,
  
  prior_a = c(0,scale),
  prior_b = c(0,scale),
  prior_c = c(0,scale),
  prior_d = c(0,scale),
  prior_e = c(0,scale),
  prior_f = c(0,scale),
  
  prior_sigma = c(0,scale),
  prior_tau1 = c(0,scale),
  prior_tau2 = c(0,scale),
  prior_tau3 = c(0,scale),
  prior_tau4 = c(0,scale),
  prior_tau5 = c(0,scale),
  prior_tau6 = c(0,scale),
  
  prior_eta = 1
)


#### Modeling ####
set.seed(NULL)
fit <- run.model("_scripts/Stan/RTF_protocol.stan", datalist = dataList, nIter = 8000, adapt_delta = 0.95)

# safe model for backup
save("fit", file ="fit_final")
load(file ="fit_final")

# model diagnostics
cat("Diagnostics for fitted model: \n", diagnose.model(fit))
pars <- c("mu_SV_90",
          "mu_SV_80",
          "mu_SV_70",
          "mu_delta_90",
          "mu_delta_80",
          "mu_delta_70",
          "sigma")
rstan::traceplot(fit, par = pars, nrow = ceiling(length(pars)/2), inc_warmup = FALSE)


#### Summarize posteriors #####
draws <- rstan::extract(fit, inc_warmup =FALSE)
digits <- 1
ROPE <- 1

# group-level effects

# RTF in SV protocol
posterior.summary(draws$mu_SV_90, digits = digits)
posterior.summary(draws$mu_SV_80, digits = digits)
posterior.summary(draws$mu_SV_70, digits = digits)

# delta effects
effects_protocol <- data.frame(delta90 = draws$mu_delta_90,
                               delta80 = draws$mu_delta_80,
                               delta70 = draws$mu_delta_70)

for (i in 1: length(effects_protocol)){
  
  lower <- rope(effects_protocol[, i], range = c(-Inf, -ROPE), ci = 1)[[4]] * 100
  trivial <- rope(effects_protocol[, i], range = c(-ROPE, ROPE), ci = 1)[[4]] * 100
  higher <- rope(effects_protocol[, i], range = c(ROPE, Inf), ci = 1)[[4]] * 100
  
  cat("\n", colnames(effects_protocol)[i], "... \n", 
      posterior.summary(effects_protocol[, i], point = "mean", digits = digits), "\n",
      "p(\u03B8 < ROPE | data)", format(round(lower, digits), nsmall = digits)[1], "\n", 
      "p(\u03B8 \u2208 ROPE | data)", format(round(trivial, digits), nsmall = digits)[1], "\n",
      "p(\u03B8 > ROPE | data)", format(round(higher, digits), nsmall = digits)[1], "\n"
  )
  
}

# SD of delta effects
posterior.summary(draws$tau_u[, 4], digits = digits)
posterior.summary(draws$tau_u[, 5], digits = digits)
posterior.summary(draws$tau_u[, 6], digits = digits)

# posterior predictive distributions (PPDs)
df_PPD <- data.frame(delta = c("90%", "80%", "70%"), mean = rep(0, 3), sd = rep (0, 3))
delta_list = c(90, 80, 70)
for (e in 1:3){
  
  PPD <- rep(0, (nIter-nWarmup)*4)
  for (i in 1:((nIter-nWarmup)*4)){
    PPD[i] <- draws$sim[i, 3+e]
  }
  MM <- mmedist(PPD, distr = "norm")$estimate
  df_PPD[e, 2:3] <- round(MM, digits = digits)
}
print(df_PPD)

# subject-level effects
delta90_s <- draws$delta_90
delta80_s <- draws$delta_80
delta70_s <- draws$delta_70
stan_plot(fit, pars = c("delta_90", "delta_80", "delta_70"), show_density=T, fill_color = 'skyblue')

# R²
digits <- 2
draws <- rstan::extract(fit, inc_warmup =FALSE)
reps_pred <- t(apply(draws$reps_pred, 1L, c))
Residuals <- t(apply(draws$Residuals, 1L, c))
R2 <- bayesian.R2(reps_pred, Residuals)
adjR2 <- 1 - ((1-R2^2)*(nSubjects-1))/(nSubjects-2-1)
posterior.summary(R2, digits = digits)
posterior.summary(adjR2, digits = digits)
