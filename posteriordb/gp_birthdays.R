cmdstanr::cmdstan_make_local(cpp_options=list(STAN_THREADS=TRUE),append=TRUE)
cmdstanr::rebuild_cmdstan()
## Load necessary libraries ##
library(rstan)
library(tidyverse)
library(readr)
library(bayesplot)
library(patchwork)
library(bridgesampling)
library(cmdstanr)
library(loo)
library(parallel)
library(foreach)
set.seed(1)
## Set up the environment ##
setwd("/scratch/work/micaleg1/pathfinder")
options(mc.cores = parallel::detectCores())
print(parallel::detectCores())
rstan_options(auto_write = TRUE)
options(pillar.neg = FALSE, pillar.subtle = FALSE, pillar.sigfig = 2)
theme_set(bayesplot::theme_default(base_family = "sans"))

## Load utility scripts ##
source("./utils/sim_pf.R")
source("./utils/lp_utils.R")

## Color palette for plots ##
set1 <- RColorBrewer::brewer.pal(7, "Set1")

## Define the function to load data, compile and run the model ##
run_birthday_model <- function(data_path, stan_file, seed) {
  ## Load data
  data <- read_csv(data_path) %>%
    mutate(date = as.Date("1968-12-31") + id,
           births_relative100 = births / mean(births) * 100)
    

  ## Data list from the original file
  data_list <- list(x=data$id,
                    y=log(data$births_relative100),
                    N=length(data$id),
                    c_f1=1.5, # factor c of basis functions for GP for f1
                    M_f1=20, # number of basis functions for GP for f1
                    J_f2=20, # number of basis functions for periodic f2
                    day_of_week=data$day_of_week,
                    monday = 0,
                    day_of_year=data$day_of_year2) # 1st March = 61 every year
  
  
  ## Compile and run the Stan model
  model <- cmdstan_model(stan_file, force_recompile = TRUE)
  pth6 <- model$pathfinder(data = data_list, init=0.1, seed = seed,
                          num_paths=10, single_path_draws=40, draws=400,
                          history_size=50, max_lbfgs_iters=100,
                          refresh=0, psis_resample=FALSE)

  fit <- model$sample(data = data_list, seed = seed,
                      chains = 4, iter_warmup = 1000, iter_sampling = 4000, thin = 1, init = pth6)

  ## Return the fit object
  return(fit)
}
fit_result <- run_birthday_model("./births_usa_1969.csv", "gpbf6.stan", seed = 1)
res <- bridge_sampler(fit_result, num_splits = 6, total_perms = 100, seed = 1, return_always = TRUE)
results <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
for (j in 1:length(res)) {
  results <- rbind(results, data.frame(logml = res[[j]]$logml, 
                                       pareto_k_numi = res[[j]]$pareto_k_numi, 
                                       pareto_k_deni = res[[j]]$pareto_k_deni,
                                       mcse_logml = res[[j]]$mcse_logml
                                       ))
}
write.csv(results, file = "pathfinderbirthdays.csv", row.names = FALSE)
print("FINISHED :)")

res <- bridge_sampler(fit_result, num_splits = 6, total_perms = 100, seed = 1, return_always = TRUE, verbose = TRUE, cores = parallel::detectCores(), pareto_smoothing_all = TRUE)
results <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
for (j in 1:length(res)) {
  results <- rbind(results, data.frame(logml = res[[j]]$logml, 
                                       pareto_k_numi = res[[j]]$pareto_k_numi, 
                                       pareto_k_deni = res[[j]]$pareto_k_deni,
                                       mcse_logml = res[[j]]$mcse_logml
                                       ))
}
write.csv(results, file = "pathfinderbirthdays_smoothed.csv", row.names = FALSE)


max_iter <- 100
for(j in 1:max_iter){
  fit_result <- run_birthday_model("./births_usa_1969.csv", "gpbf6.stan", seed = j)
  res <- bridgesampling:::bridge_sampler.CmdStanMCMC(fit_result, num_splits = 2, total_perms = 1, seed = j, return_always = TRUE)
  results <- rbind(results, data.frame(logml = res$logml, 
                                       pareto_k_numi = res$pareto_k_numi, 
                                       pareto_k_deni = res$pareto_k_deni,
                                       mcse_logml = res$mcse_logml
                                       ))
}
write.csv(results, file = "pathfinderbirthdays_bruteforce.csv", row.names = FALSE)
print("FINISHED :)")

results <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
max_iter <- 100
for(j in 1:max_iter){
  fit_result <- run_birthday_model("./births_usa_1969.csv", "gpbf6.stan", seed = j)
  res <- bridge_sampler(fit_result, num_splits = 2, total_perms = 1, seed = j, return_always = TRUE, verbose = TRUE, cores = parallel::detectCores(), pareto_smoothing_all = TRUE)
  results <- rbind(results, data.frame(logml = res$logml, 
                                       pareto_k_numi = res$pareto_k_numi, 
                                       pareto_k_deni = res$pareto_k_deni,
                                       mcse_logml = res$mcse_logml
                                       ))
}
write.csv(results, file = "pathfinderbirthdays_bruteforce_smoothed.csv", row.names = FALSE)
