## ---------------------------
##
## Project name: Data Collection for Strength-endurance Profiles: Can Assessments Be Completed within a Single Session?
##
## Purpose of script: plots for the paper
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
library(dagitty)
library(ggdag)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(readxl)
library(rstan)
library(tidyr)


#### Load functions ####
source("_scripts/functions.R")


#### Theme presets ####
theme_pub <- theme_bw() + 
  theme(
    text = element_text(family = "Times New Roman"),
    legend.position = "None",
    axis.title.y = element_text(margin = margin(r = 9)),
    axis.title.x = element_text(margin = margin(t = 9)),
    axis.title = element_text(size = 11), 
    panel.grid = element_blank()
  )


#### Figure 2 (DAG) ####
coordinates <- tribble(
  ~name, ~x, ~y,
  "P", 1, 1,
  "RTF", 2, 1,
  "T", 1.5, 0, 
  "holder", 1.5, 1
)

DAG <- dagify(
  RTF ~ P,
  eval(parse(text = "holder ~ T")),
  eval(parse(text = "RTF ~ T")),
  coords = coordinates
)

Figure2 <- ggplot(DAG, aes(x = x, y = y, xend = xend, yend = yend)) + 
geom_dag_point(data = function(x) filter(x, (name == "P") | 
                                           (name == "T")),
               size = 12, shape = 15
               ) +
geom_dag_point(data = function(x) filter(x, (name != "holder") & 
                                           (name != "P") & 
                                           (name != "T")),
               size = 12, shape = 19
               ) +
geom_dag_edges() + 
geom_dag_text(data = function(x) filter(x, (name != "holder")),
              size = 3
              ) +
theme_dag()


#### Figure 3A (data) ####

# data frame
data_directory <- ('./_data/rawdata.xlsx')
rawdata <- as.data.frame(read_xlsx(data_directory, col_names = T))
nSubjects <- nrow(rawdata)

plot_data <- rawdata[1:nSubjects, (ncol(rawdata)-5):ncol(rawdata)]
ID <- seq(1, nSubjects, 1)
plot_data <- cbind(ID, 
                   plot_data)
plot_data <- reshape(plot_data, direction = "long", 
                varying = list(c(colnames(plot_data[2:4])), c(colnames(plot_data[5:7]))),
                timevar = "load", 
                times = c("90", "80", "70"),
                v.names = c("SV", "MV"),
                idvar = c("ID")
                )
plot_data$ID <- as.factor(plot_data$ID)
plot_data$load <- as.factor(plot_data$load)
plot_data <- gather(plot_data, protocol, reps, colnames(plot_data[, (ncol(plot_data)-1):ncol(plot_data)]), factor_key = TRUE)

# plot
ylim <- c(0, 2 * ceiling(max(plot_data$reps)/2))
loads <- levels(plot_data$load)

for (i in 1:length(loads)){
  
  plot_tmp <- ggplot(plot_data[plot_data$load == loads[i], ], aes(x = protocol, y = reps)) + 
  lemon::geom_pointpath(aes(group = ID),
                        position = position_jitter(width = 0.15, height = 0.15),
                        distance = 0, size = 2, alpha = 0.6, color = "grey60"
  ) +
  xlab("Protocol") +
  ylab("RTF") +
  theme_pub +
  scale_x_discrete(expand = c(0, 0.4)) +
  scale_y_continuous(limits = ylim, breaks = seq(0, ylim[2], 4))
  
  assign(paste0("plot_", loads[i]), plot_tmp)
  
}


#### Figure 3B (Posteriors) ####
load(file = "fit_final") # "analysis.R" has to be run before plotting

# extract posteriors
draws <- rstan::extract(fit, inc_warmup = FALSE)

SV_90 <- draws$mu_SV_90
SV_80 <- draws$mu_SV_80
SV_70 <- draws$mu_SV_70
delta_90 <- draws$mu_delta_90
delta_80 <- draws$mu_delta_80
delta_70 <- draws$mu_delta_70
sigma <- draws$sigma

adj <- 1
bw <- "SJ"
kernel <- "gaussian"
n_dens <- 1000
ci <- 0.95

posterior_list <- c("delta_90", "delta_80", "delta_70")
xlabs <- c(
  expression(mu[Delta*a]),
  expression(mu[Delta*b]),
  expression(mu[Delta*c])
)

for (i in 1: length(posterior_list)){
  
  post_tmp <- eval(parse(text = posterior_list[i]))
  dens_tmp <- density(post_tmp, bw = bw, adjust = adj, kernel = kernel, n = n_dens) 
  hdi_tmp <- unlist(hdi(post_tmp, ci = ci))
  plot_data <- data.frame(x = dens_tmp$x, y = dens_tmp$y)
  plot_hdi <- data.frame(x = dens_tmp$x[(dens_tmp$x > hdi_tmp[2]) & (dens_tmp$x < hdi_tmp[3])], 
                         y = dens_tmp$y[(dens_tmp$x > hdi_tmp[2]) & (dens_tmp$x < hdi_tmp[3])])
  
  plot_tmp <- ggplot(data = plot_data, aes(x = x, y = y)) + 
    geom_area(data = plot_hdi, aes(x=x, y=y), color = "grey60", alpha = 0.3) +
    geom_line(size = 1) +
    xlab(xlabs[i]) +
    geom_vline(xintercept = 0, size = 0.5) +
    geom_vline(xintercept = c(-1, 1), linetype = "22", size = 0.5) +
    scale_x_continuous(limits = c(-3, 5), breaks = seq(-3, 5, 1)) +
    theme_pub +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  assign(paste0("plot_", posterior_list[i]), plot_tmp)
  
}

Figure3 <- ggarrange(plot_90, plot_delta_90, plot_80, plot_delta_80, plot_70, plot_delta_70,
                     labels = c("A", "", "B", "", "C", ""),
                     label.x = 0.01,
                     label.y = 0.99,
                     ncol = 2,
                     nrow = 3)


#### Figure 4 ####

# functions
fun.Ex2 <- function(x, a, b){
  y = a * exp(b * x)
  return(y)
}

fun.Rec <- function(x, a, b){
  y = 1/(a + b * x)
  return(y)
}

# load data
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


# load parameters
load("_files/final_models") # "replication_analysis.R" has to be run before plotting
draws_Ex2 <- rstan::extract(fit_Ex2, inc_warmup = FALSE)
draws_Rec <- rstan::extract(fit_Rec, inc_warmup = FALSE)

# find Subjects with largest and lowest residuals
Residuals_Ex2 <- matrix(0, nrow = nSubjects, ncol = 5)
for (i in 1:nSubjects){
  
  for (t in 1:4){
    
    Residuals_Ex2[i, t] <- mean(draws_Ex2$Residuals[, i, t])
    
  }
  
  Residuals_Ex2[i, 5] <- sum(Residuals_Ex2[i,]^2)
  
}
min_Ex2 <- which(Residuals_Ex2[, 5] == min(Residuals_Ex2[, 5]))
max_Ex2 <- which(Residuals_Ex2[, 5] == max(Residuals_Ex2[, 5]))

Residuals_Rec <- matrix(0, nrow = nSubjects, ncol = 5)
for (i in 1:nSubjects){
  
  for (t in 1:4){
    
    Residuals_Rec[i, t] <- mean(draws_Rec$Residuals[, i, t])
    
  }
  
  Residuals_Rec[i, 5] <- sum(Residuals_Rec[i,]^2)
  
}
min_Rec <- which(Residuals_Rec[, 5] == min(Residuals_Rec[, 5]))
max_Rec <- which(Residuals_Rec[, 5] == max(Residuals_Rec[, 5]))

cat("Ex2: \n", 
    "Best fit: ", min_Ex2, "Worst fit: ", max_Ex2, "\n",
    "Rec: \n", 
    "Best fit: ", min_Rec, "Worst fit: ", max_Rec, "\n"
    )

# select sunject IDs with min and max residuals
min_select <- 10
max_select <- 7

# build data frame
plot_data_min <- data.frame(x = reps[min_select, ], y = loads[min_select, ])
plot_data_max <- data.frame(x = reps[max_select, ], y = loads[max_select, ])

# parameter means
a_min_Ex2 <- mean(draws_Ex2$a_rec[, min_select])
a_max_Ex2 <- mean(draws_Ex2$a_rec[, max_select])
b_min_Ex2 <- mean(draws_Ex2$b_rec[, min_select])
b_max_Ex2 <- mean(draws_Ex2$b_rec[, max_select])

a_min_Rec <- mean(draws_Rec$a_rec[, min_select])
a_max_Rec <- mean(draws_Rec$a_rec[, max_select])
b_min_Rec <- mean(draws_Rec$b_rec[, min_select])
b_max_Rec <- mean(draws_Rec$b_rec[, max_select])

plot_base <- ggplot(data = data.frame(x = 0, y = 0), aes(x=x, y=y)) +
  geom_point(size = 3) +
  xlim(0, 25) +
  ylim(70, 160) +
  ylab("Load (kg)") + 
  xlab("RTF (n)") +
  theme_pub +
  theme(plot.title = element_text(size = 12, hjust = 0.5))

colors <- c("grey60", "black")
alphas <- c(0.5, 0.6)
linetypes <- c("solid", "22")
functions <- c("Ex2", "Rec")
subj <- c("min", "max")
ID <- c("A", "B")

for (s in 1:length(ID)){
  
  plot.tmp <- plot_base
  
  for (m in 1:2){

    args <- list(a = eval(parse(text = paste0("a_", subj[s], "_", functions[m]))), 
                 b = eval(parse(text = paste0("b_", subj[s], "_", functions[m]))))
    
    plot.tmp <- plot.tmp +
      stat_function(fun = eval(parse(text = paste0("fun.", functions[m]))), 
                    args = args, 
                    color = colors[m], 
                    alpha = alphas[m], 
                    linetype = linetypes[m], 
                    size =1.5)

  }
  
  plot.tmp <- plot.tmp +
    geom_point(data = eval(parse(text = paste0("plot_data_", subj[s]))), 
               aes(x=x, y=y), 
               size = 3,
               stroke = 2,
               shape = 3,
               color = "black", 
               alpha = 1) +
    ggtitle(paste0("Strength-endurance profile of athlete ", ID[s]))
  
  assign(paste0("plot_", ID[s], i), plot.tmp)
  
}

Figure4 <- ggarrange(eval(parse(text = paste0("plot_A", i))), eval(parse(text = paste0("plot_B", i))),
                      ncol = 2,
                      nrow = 1)


#### Print figures ####
from <- 2
n_plots <- 4

output_width <- c(90, 190, 190) # in mm
output_height <- c(output_width[1] * 0.8, output_width[2] * 3 / 2, output_width[3] / 2)

for (i in from:n_plots){
  
  file_name <- paste0("_plots/Figure", i, ".png")
  ggsave(filename = file_name, plot = eval(parse(text = paste0("Figure", i))),  
       width = output_width[i-1], height = output_height[i-1],   
       units = "mm", device="png", dpi = 1200
       , bg ="white"
       )
  
}
