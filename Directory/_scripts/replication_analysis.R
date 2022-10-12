## ---------------------------
##
## Project name: Data Collection for Strength-endurance Profiles: Can Assessments Be Completed within a Single Session?
##
## Purpose of script: analysis of model functions (Lin, Ex2, Rec, Ex3, Crt)
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
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(loo)
library(matrixStats)
library(posterior)
library(readxl)
library(reshape2)
library(rstan)
library(utils)
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

loads <- matrix(0, nrow = nSubjects, ncol = 4)
loads[, 1] <- rawdata$`1RM_true(kg)`
loads[, 2] <- rawdata$`load90(n)`
loads[, 3] <- rawdata$`load80(n)`
loads[, 4] <- rawdata$`load70(n)`

reps <- rawdata[, (ncol(rawdata)-2):ncol(rawdata)]
colnames(reps) <- NULL
reps <- as.matrix(reps)
reps <- cbind(rep(1, nSubjects), reps)


#### Prior sensitivity analysis ####

# modeling

scale <- c(5, 10, 20, 30, 40, 50)
model_types <- c("Lin", "Ex2", "Rec", "Ex3", "Crt")
model_list <- rep(0, length(scale))

t <- 1 # select index from model_types to choose a Stan script

diagnostics <- data.frame(scale = scale, divergences = rep(0, length(scale)), low_bfmi = rep(0, length(scale)))

if (t == 4){
  init = function(){list(alpha = 1, beta = -0.5, gamma = -1, tau_u = c(1, 1, 1))}
  adapt_delta = 0.95
} else if (t == 5){
  init = function(){list(alpha = 20, beta = -5, gamma = -5, tau_u = c(1, 1, 1))}
  adapt_delta = 0.95
} else{
  init = "random"
  adapt_delta = 0.8
}

for (i in 1:length(scale)){
  
  model_name = paste0("fit_", model_types[t], "_", scale[i])
  model_list[i] <- model_name
  
  dataList <- list(
    nSubjects = nSubjects,
    load = loads,
    reps = reps,
    
    prior_a = c(0, scale[i]),
    prior_b = c(0, scale[i]),
    prior_c = c(0, scale[i]),
    
    prior_sigma = c(0, scale[i]),
    prior_tau1 = c(0, scale[i]),
    prior_tau2 = c(0, scale[i]),
    prior_tau3 = c(0, scale[i]),
    
    prior_eta = 1
  )
  
  fit_tmp <- run.model(paste0("_scripts/Stan/", model_types[t], ".stan"), dataList, adapt_delta = adapt_delta, seed = 2468, init = init)
  
  assign(model_name, fit_tmp)
  
  save("fit_tmp", file = paste0("_files/", model_name))
  
  diagnostics[i, 2] <- get_num_divergent(fit_tmp)
  diagnostics[i, 3] <- if (length(get_low_bfmi_chains(fit_tmp)) == 0) {0} else {get_low_bfmi_chains(fit_tmp)}
  
}


# sensitivity analysis (full parameter spectrum)

pars <- names(eval(parse(text = model_list[1])))
pars <- pars[(pars != "lp__")]
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
                "overlap: ", round(min(df_PSA[, i]*100), 1), "\n"
  )
  
}

print(diagnostics) # minimalist model diagnostics


#### Modeling ####

# data list using appropriate prior scales
scale <- 40 # select appropriate scale from prior sensitivity analysis

dataList <- list(
  nSubjects = nSubjects,
  load = loads,
  reps = reps,
  
  prior_a = c(0, scale),
  prior_b = c(0, scale),
  prior_c = c(0, scale),
  
  prior_sigma = c(0, scale),
  prior_tau1 = c(0, scale),
  prior_tau2 = c(0, scale),
  prior_tau3 = c(0, scale),
  
  prior_eta = 1
)

model_types <- c("Lin", "Ex2", "Rec", "Ex3", "Crt")

pars <- c("alpha_rec", "beta_rec", "sigma_rec", "gamma_rec")

for (t in 1:(length(model_types))){
  
  model_name = paste0("fit_", model_types[t])
  
  if (t == 3){
    init <- "random"
    adapt_delta <- 0.8
    pars.tmp <- c(pars[1:2], "SEE")
  } else if (t == 4){
    init <- function(){list(alpha = 1, beta = -0.5, gamma = -1, tau_u = c(1, 1, 1))}
    adapt_delta <- 0.95
    pars.tmp <- pars
  } else if (t == 5){
    init <- function(){list(alpha = 20, beta = -5, gamma = -5, tau_u = c(1, 1, 1))}
    adapt_delta <- 0.95
    pars.tmp <- pars
  } else{
    init <- "random"
    adapt_delta <- 0.8
    pars.tmp <- pars[1:3]
  }
  
  fit_tmp <- run.model(paste0("_scripts/Stan/", model_types[t], ".stan"), dataList, nIter = 8000, adapt_delta = adapt_delta, init = init)
  assign(model_name, fit_tmp)
  
  cat("\n Diagnostic reults for ", model_types[t], "\n")
  diagnose.model(fit_tmp)
  rstan::traceplot(fit_tmp, par = pars.tmp, nrow = ceiling(length(pars)/2), inc_warmup = FALSE)
  
}

# safe models
save(list=c(
  "fit_Lin", 
  "fit_Ex2",
  "fit_Rec",
  "fit_Ex3",
  "fit_Crt"
  ), file ="./_files/final_models")

# load models
load(file ="./_files/final_models")

# model diagnostics
diagnose.model(fit)
pars <- c("alpha_rec",
          "beta_rec",
          "gamma_rec",
          "sigma_rec")
rstan::traceplot(fit_Ex3, par = pars, nrow = ceiling(length(pars)/2), inc_warmup = FALSE)


#### Model evaluation ####

model_types <- c("Lin", "Ex2", "Rec", "Ex3", "Crt")

for (t in 1:length(model_types)){
  
  model <- eval(parse(text = paste0("fit_", model_types[t])))
  LL <- extract_log_lik(model)
  assign(paste0("loo", t), loo(LL, is_method = "sis", r_eff = NA))
  
}

ranking <- substr(rownames(loo_compare(loo1, loo2, loo3, loo4, loo5)), 6, 6)

# looic
choose_top <- 3
names <- ranking
loo_top <- eval(parse(text = paste0("loo", names[1])))
loo_a <- eval(parse(text = paste0("loo", names[2])))
loo_b <- eval(parse(text = paste0("loo", names[3])))
loo_c <- eval(parse(text = paste0("loo", names[4])))
loo_d <- eval(parse(text = paste0("loo", names[5])))

looic_matrix <- data.frame(model = names, 
                           LOOIC_diff = rep(0, 5), 
                           LOOIC_diff_se = rep(0, 5), 
                           relative = rep(0, 5))

looic_matrix[2,2] <- -2 * loo_compare(loo_top, loo_a)[2,1] * (if (as.numeric(substr(row.names(loo_compare(loo_top, loo_a))[1], 6, 6)) == 2){-1} else {1})
looic_matrix[3,2] <- -2 * loo_compare(loo_top, loo_b)[2,1] * (if (as.numeric(substr(row.names(loo_compare(loo_top, loo_b))[1], 6, 6)) == 2){-1} else {1})
looic_matrix[4,2] <- -2 * loo_compare(loo_top, loo_c)[2,1] * (if (as.numeric(substr(row.names(loo_compare(loo_top, loo_c))[1], 6, 6)) == 2){-1} else {1})
looic_matrix[5,2] <- -2 * loo_compare(loo_top, loo_d)[2,1] * (if (as.numeric(substr(row.names(loo_compare(loo_top, loo_d))[1], 6, 6)) == 2){-1} else {1})

looic_matrix[2,3] <- 2 * loo_compare(loo_top, loo_a)[2,2]
looic_matrix[3,3] <- 2 * loo_compare(loo_top, loo_b)[2,2]
looic_matrix[4,3] <- 2 * loo_compare(loo_top, loo_c)[2,2]
looic_matrix[5,3] <- 2 * loo_compare(loo_top, loo_d)[2,2]

for (i in 1:5){
  looic_matrix[i,4] <- looic_matrix[i,2]/looic_matrix[i,3]
}

looic_matrix[, 2:4] <- round(looic_matrix[, 2:4], digits = 1)

looic_matrix

# R²
model_types <- c("Lin", "Ex2", "Rec", "Ex3", "Crt")
digits <- 2
CrInt <- 0.95
df_R2 <- data.frame(model = model_types, R2 = rep(0, length(model_types)), adjR2 = rep(0, length(model_types)))

for (t in 1: length(model_types)){
  
  if (t > 3){
    k <- 1 + 3 * (2 + nSubjects)
  } else{
    k <- 1 + 2 * (2 + nSubjects)
  }
  
  model <- eval(parse(text = paste0("fit_", model_types[t])))
  draws <- rstan::extract(model, inc_warmup =FALSE)
  load_pred <- t(apply(draws$load_pred, 1L, c))
  Residuals <- t(apply(draws$Residuals, 1L, c))
  R2 <- bayesian.R2(load_pred, Residuals, adjR2 = FALSE)
  adjR2 <- bayesian.R2(load_pred, Residuals, adjR2 = TRUE, N = length(loads), k = k)
  
  R2 <- posterior.summary(R2, ci = CrInt, digits = digits)
  
  adjR2 <- posterior.summary(adjR2, ci = CrInt, digits = digits)
  
  df_R2[t, 2] <- R2
  df_R2[t, 3] <- adjR2
  
}

print(df_R2[ranking, ])

# compare R² (posterior overlap)
comparisons <- as.data.frame(t(combn(model_types, 2)))
comparisons$R2_p_OL <- 0

for (c in 1:nrow(comparisons)){
  
  for (m in 1:2){
    
    model <- eval(parse(text = paste0("fit_", comparisons[c, m])))
    draws <- rstan::extract(model, inc_warmup =FALSE)
    load_pred <- t(apply(draws$load_pred, 1L, c))
    Residuals <- t(apply(draws$Residuals, 1L, c))
    if ((comparisons[c, m] == "Ex3") | (comparisons[c, m] == "Crt")){
      k <- 1 + 3 * (2 + nSubjects)
    } else{
      k <- 1 + 2 * (2 + nSubjects)
    }
    assign(paste0("R2_", m), bayesian.R2(load_pred, Residuals, adjR2 = FALSE, N = length(loads), k = k))
    
  }
  
  comparisons[c, 3] <- compare_densities(R2_1, R2_2)
  
}

print(comparisons)
