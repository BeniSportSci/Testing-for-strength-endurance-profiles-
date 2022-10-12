## ---------------------------
##
## Project name: Data Collection for Strength-endurance Profiles: Can Assessments Be Completed within a Single Session?
##
## Purpose of script: prior sensitivity analysis for the paper
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
library(greekLetters)
library(loo)
library(matrixStats)
library(posterior)
library(readxl)
library(rstan)
library(writexl)

#### Functions ####
source( "./_scripts/functions.R")


#### Seed ####
set.seed(1234)


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

trialB <- t(matrix(c(0, 1, 0), 6, nSubjects))
trialB <- array(trialB, dim = c(nSubjects, 3, 2))

trialC <- t(matrix(c(0, 0, 1), 6, nSubjects))
trialC <- array(trialC, dim = c(nSubjects, 3, 2))

# define different scales
scale <- c(1, 5, 10, 20, 30, 40, 50)
model_list <- rep(0, length(scale))
for (i in 1:length(model_list)){
  model_list[i] <- paste0("fit_", scale[i])
}

for (i in 1:length(scale)){
  
  def_prior <- c(0, scale[i])
  
  dataList <- list(
    nSubjects = nSubjects,
    protocol = protocol,
    reps = floor(reps),
    dummy1 = trialB,
    dummy2 = trialC,
    
    prior_a = def_prior,
    prior_b = def_prior,
    prior_c = def_prior,
    prior_d = def_prior,
    prior_e = def_prior,
    prior_f = def_prior,
    
    prior_sigma = def_prior,
    prior_tau1 = def_prior,
    prior_tau2 = def_prior,
    prior_tau3 = def_prior,
    prior_tau4 = def_prior,
    prior_tau5 = def_prior,
    prior_tau6 = def_prior,
    
    prior_eta = 1
  )
  
  model_name = paste0("fit_", scale[i])
  
  fit_tmp <- run.model("_scripts/Stan/RTF_protocol.stan", dataList, adapt_delta = 0.95, save_as = NA)
  
  assign(model_name, fit_tmp)
  
}

# model diagnostics
for (i in 1:length(model_list)){
  cat(model_list[i], ": \n")
  diagnose.model(eval(parse(text = model_list[i])))
}

# sensitivity analysis w/ plots
pars <- c("mu_SV_90",
          "mu_SV_80",
          "mu_SV_70",
          "mu_delta_90",
          "mu_delta_80",
          "mu_delta_70",
          "sigma")

for (i in 2: length(model_list)){
  
  model1 <- model_list[i-1]
  model2 <- model_list[i]
  
  plot.prior.sens(model1, model2, parameter = pars, plottype = "dens", smooth = 1)
  
}

# sensitivity analysis (full parameter spectrum)
pars <- names(eval(parse(text = model_list[1])))
pars <- pars[pars != "lp__"]
df_PSA <- as.data.frame(matrix(0, length(pars), length(model_list)))
df_PSA[, 1] <- pars
colnames(df_PSA)[1] <- "parameter"

for (i in 2: length(model_list)){
  
  model1 <- model_list[i-1]
  model2 <- model_list[i]
  colnames(df_PSA)[i] <- paste0(model1, "_vs_", model2)
  
  for (j in 1: length(pars)){
    
    df_PSA[j, i] <- prior.sens(model1, model2, parameter = pars[j], smooth = 1)
    
  }
  
  cat("Smallest overlap in ", colnames(df_PSA)[i], ": \n", 
      "parameter: ", df_PSA[which.min(df_PSA[, i]), 1], "\n",
      "overlap: ", round(min(df_PSA[, i]), 3)
  )
  
}
