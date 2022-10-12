## ---------------------------
##
## Project name: Data Collection for Strength-endurance Profiles: Can Assessments Be Completed within a Single Session?
##
## Purpose of script: gather all functions used in scripts
##
## Author: Benedikt Mitter MSc
##
## Date Created: 2022-Oct-12
##
## Copyright (c) Benedikt Mitter, 2022
## Email: benedikt.mitter@univie.ac.at
##
## ---------------------------


#### Run Stan model ####

run.model <- function(model, datalist, nChains = 4, nIter = 4000, nWarmup = floor(nIter/2), nThin = 1, adapt_delta = 0.8, init = "random", seed = NA, save_as = NA){
  
  require(rstan)
  
  # Stan setup
  rstan_options(auto_write = TRUE)
  if (parallel::detectCores() > nChains){
    options(mc.cores = nChains)
  } else {
    options(mc.cores = (parallel::detectCores()-1))
  }
  
  # clean object names
  if (exists("fit.tmp")){
    rm(fit)
  }
  
  # Stan import
  cat("Estimating", model, "model... \n")
  startTime = Sys.time()
  cat("Calling", nChains, "simulations in Stan... \n")
  
  fit.tmp = rstan::stan(model, 
             data    = datalist, 
             chains  = nChains,
             iter    = nIter,
             warmup  = nWarmup,
             thin    = nThin,
             control = list(adapt_delta = adapt_delta),
             init    = init,
             seed    = seed
  )
  
  cat("Finishing", model, "model simulation ... \n")
  endTime = Sys.time()
  cat("It took",as.character.Date(endTime - startTime), "\n")
  
  return(fit.tmp)
  
  # safe model
  if (!is.na(save_as)){
    save(fit.tmp, file = paste0(save_as))
  }

}


#### Model diagnostics ####

diagnose.model <- function(model, save_as = NA){
  
  require(rstan)
  require(posterior)
  require(writexl)
  model <- model
  chain_diagnostics <- posterior::summarise_draws(model, "rhat", "ess_bulk", "ess_tail")
  Rhat_max <- chain_diagnostics[which.max(chain_diagnostics$rhat),]
  ess_bulk_min <- chain_diagnostics[which.min(chain_diagnostics$ess_bulk),]
  ess_tail_min <- chain_diagnostics[which.min(chain_diagnostics$ess_tail),]
  diagnostic_table <- data.frame(diagnostic = c("max Rhat", 
                                                "min bluk ess", 
                                                "min tail ess"), 
                                 parameter = c(unlist(Rhat_max[1]),
                                               unlist(ess_bulk_min[1]),
                                               unlist(ess_tail_min[1])),
                                 value = c(unlist(Rhat_max[2]),
                                           unlist(ess_bulk_min[3]),
                                           unlist(ess_tail_min[4]))
                                 )
  
  
  if (!is.na(save_as)){
    chain_diagnostics <- na.omit(chain_diagnostics)
    chain_diagnostics <- chain_diagnostics[order(-chain_diagnostics$ess_bulk),]
    writexl::write_xlsx(chain_diagnostics, path = save_as, col_names = TRUE)
  }
  
  divergences <- get_num_divergent(model)
  bfmi <- if (length(get_low_bfmi_chains(model)) == 0) {0} else {get_low_bfmi_chains(model)}
  cat("divergences: ", divergences, "\n", "BFMI: ", bfmi, "\n")
  print(diagnostic_table)
}


#### Prior Sensitivity Analysis ####

plot.prior.sens <- function(model1, model2, parameter, plottype = "dens", smooth = 1, low = NA, high = NA){
  
  require(rstan)
  require(sfsmisc)
  require(ggplot2)
  require(ggpubr)
  
  output <- tryCatch(
    {
      
      plot_list <- list()
      
      for (i in 1:length(parameter)){
        
        # prepare data frame
        posterior1 <- unlist(rstan::extract(object = eval(parse(text = model1)), pars = parameter[i], inc_warmup = FALSE))
        posterior2 <- unlist(rstan::extract(object = eval(parse(text = model2)), pars = parameter[i], inc_warmup = FALSE))
        PSA.tbl <- data.frame(model = c(rep(model1, length(posterior1)), rep(model2, length(posterior2))), 
                              value = c(posterior1, posterior2)
        )
        PSA.tbl$model <- as.factor(PSA.tbl$model)
        
        # prepare x axis
        if (is.na(low)){
          xlow <- min(c(posterior1, posterior2))
        }
        if (is.na(high)){
          xhigh <- max(c(posterior1, posterior2))
        }
        
        # prepare MAPs and overlap
        lower <- min(c(posterior1, posterior2)) - 1 
        upper <- max(c(posterior1, posterior2)) + 1
        
        dens1 <- density(posterior1, adjust = smooth, kernel = "g", from = lower, to = upper)
        dens2 <- density(posterior2, adjust = smooth, kernel = "g", from = lower, to = upper)
        MAP1 <- dens1$x[which.max(dens1$y)]
        MAP2 <- dens2$x[which.max(dens2$y)]
        
        dens <- data.frame(x=dens1$x, a=dens1$y, b=dens2$y)
        
        dens$w <- pmin(dens$a, dens$b)
        
        total <- sfsmisc::integrate.xy(dens$x, dens$a) + sfsmisc::integrate.xy(dens$x, dens$b)
        intersection <- sfsmisc::integrate.xy(dens$x, dens$w)
        
        # compute overlap
        overlap <- 2 * intersection / total
        if (overlap < 0.001){
          OL <- "Overlap < 0.001"
        } else {
          OL <- paste0("Overlap = ", round(overlap, 3))
        }
        
        # plotting
        if (plottype == "dens"){
          plot_tmp <- ggplot(PSA.tbl, aes(x = value, color = model, fill = model, alpha = model)) +
            theme_bw() +
            ggtitle(paste0(model1, " vs. ", model2, ": ", parameter[i], "\n\n", "(", OL, ")")) +
            theme(plot.title = element_text(hjust = 0.5,
                                            size = 12
                                            ),
                  panel.grid = element_blank(),
                  axis.title.y =  element_blank(),
                  axis.ticks.y =  element_blank(),
                  axis.text.y =  element_blank()
                  ) +
            xlim(xlow, xhigh) +
            geom_density(adjust = smooth, kernel = "g", n = 512) +
            scale_alpha_manual(values = c(0.8, 0.5)) +
            scale_color_manual(values = c("blue", "red")) +
            scale_fill_manual(values = c("blue", "red")) +
            geom_vline(xintercept = MAP1, color = "black", linetype = "32") +
            geom_vline(xintercept = MAP2, color = "black", linetype = "32") 
          
        } else if (plottype == "hist"){
          plot_tmp <- hist(posterior1, col = rgb(0, 0, 1, 0.8), breaks = 50) +
            hist(posterior2, col = rgb(1, 0, 0, 0.5), breaks = 50, add = T)
        } else {
          print("Error! Plot type unknown!")
        }
        
        assign(paste0("plot_", i), plot_tmp)
        
        plot_list[[i]] <- eval(parse(text = paste0("plot_", i)))
        
      }
      
      plot_PSA <- ggarrange(plotlist = plot_list)
      
      print(plot_PSA)
      
    }, 
    error = function(e){
      message("Syntax error!")
      message(e)
      message()
      return(NULL)
    },
    warning = function(w){
      message("Syntax warning!")
      message(w)
      message()
      return(output)
    }
  )
  
  return(output)
  
}


prior.sens <- function(model1, model2, parameter, smooth = 1){
  
  require(rstan)
  require(sfsmisc)
  
  output <- rep(0, length(parameter))
  
  for (i in 1:length(parameter)){
    
    # prepare data frame
    posterior1 <- unlist(rstan::extract(object = eval(parse(text = model1)), pars = parameter[i], inc_warmup = FALSE))
    posterior2 <- unlist(rstan::extract(object = eval(parse(text = model2)), pars = parameter[i], inc_warmup = FALSE))
    
    # prepare MAPs and overlap
    lower <- min(c(posterior1, posterior2)) - 1 
    upper <- max(c(posterior1, posterior2)) + 1
    
    dens1 <- density(posterior1, adjust = smooth, kernel = "g", from = lower, to = upper)
    dens2 <- density(posterior2, adjust = smooth, kernel = "g", from = lower, to = upper)
    MAP1 <- dens1$x[which.max(dens1$y)]
    MAP2 <- dens2$x[which.max(dens2$y)]
    
    dens <- data.frame(x=dens1$x, a=dens1$y, b=dens2$y)
    
    dens$w <- pmin(dens$a, dens$b)
    
    total <- sfsmisc::integrate.xy(dens$x, dens$a) + sfsmisc::integrate.xy(dens$x, dens$b)
    intersection <- sfsmisc::integrate.xy(dens$x, dens$w)
    
    # compute overlap
    output[i] <- 2 * intersection / total
    
  }
  
  return(output)
  
}


#### Summarize posterior ####
posterior.summary <- function(x, point = "mean", uncertainty = "hdi", ci = 0.95, digits = 2, correct = 0){
  
  if (!require(bayestestR)) install.packages('bayestestR')
  library(bayestestR)
  
  if (point == "mean"){
    pe <- mean(x)
  } else if (point == "median"){
    pe <- median(x)
  } else if (point == "map"){
    pe <- map_estimate(x)
  }
  
  pe <- format(round(pe, digits), nsmall = digits)[1]
  
  if (uncertainty == "hdi"){
    interval <- round(hdi(x, ci = ci), digits)
    interval <- paste0("[", format(interval[[2-correct]], nsmall = digits)[1], ", ", format(interval[[3-correct]], nsmall = digits)[1], "]")
  } else if (point == "ci"){
    interval <- round(quantile(x, probs = c((1-ci)/2, (1+ci)/2)), digits)
    interval <- paste0("[", format(interval[1], nsmall = digits)[1], ", ", format(interval[2], nsmall = digits)[1], "]")
  } 
  
  
  output <- paste0(point, " [", ci*100, "% HDI]: ", pe, " ", interval)
  
  return(output)
}


#### Bayesian R2 ####
bayesian.R2 <- function(pred, res, adjR2 = FALSE, N, k){
  require(matrixStats)
  var_pred <- matrixStats::rowVars(pred)
  var_res <- matrixStats::rowVars(res)
  R2 <- var_pred/(var_pred + var_res)
  
  if (adjR2 == TRUE){
    R2 <- 1 - ((1-R2^2)*(N-1))/(N-k-1)
  }
  
  return(R2)
}


#### Posterior overlap ####
compare_densities <- function(posterior1, posterior2, digits = 3) { 
  
  require(sfsmisc)
  
  lower <- min(c(posterior1, posterior2)) - 1 
  upper <- max(c(posterior1, posterior2)) + 1
  
  # generate kernel densities
  dmodel1 <- density(posterior1, from=lower, to=upper)
  dmodel2 <- density(posterior2, from=lower, to=upper)
  dens <- data.frame(x=dmodel1$x, a=dmodel1$y, b=dmodel2$y)
  
  # calculate intersection densities
  dens$w <- pmin(dens$a, dens$b)
  
  # integrate areas under curves
  total <- integrate.xy(dens$x, dens$a) + integrate.xy(dens$x, dens$b)
  intersection <- integrate.xy(dens$x, dens$w)
  
  # compute overlap
  overlap <- 2 * intersection / total
  
  
  
  if (overlap < 0.001){
    "< 0.001"
  }else{
    round(overlap, digits)
  }
}
