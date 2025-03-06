##Note: The bridgesampling version of CmdstanR must be installed, comment the line below if already installed
remotes::install_github("stan-dev/cmdstanr@bridge_sampler-method")
Sys.setenv(GITHUB_PAT = "YOUR_TOKEN")
library(dplyr)
library(readr)
cmdstanr::cmdstan_make_local(cpp_options=list(STAN_THREADS=TRUE),append=TRUE)
cmdstanr::rebuild_cmdstan()
setwd("posteriordb")
library(rstan)
library(parallel)
library(foreach)
options(mc.cores = parallel::detectCores())
library(bridgesampling)
rstan_options(auto_write = TRUE)
library(cmdstanr)
library(bayesplot)
library(posteriordb)
library(posterior)
source("./utils/sim_pf.R")
source("./utils/lp_utils.R")
set.seed(123)

run_birthday_model <- function(data_path, stan_file, seed) {
  data <- read_csv(data_path) %>%
    mutate(date = as.Date("1968-12-31") + id,
           births_relative100 = births / mean(births) * 100)
    

  data_list <- list(x=data$id,
                    y=log(data$births_relative100),
                    N=length(data$id),
                    c_f1=1.5, # factor c of basis functions for GP for f1
                    M_f1=20, # number of basis functions for GP for f1
                    J_f2=20, # number of basis functions for periodic f2
                    day_of_week=data$day_of_week,
                    monday = 0,
                    day_of_year=data$day_of_year2) # 1st March = 61 every year
  
  
  model <- cmdstan_model(stan_file, force_recompile = TRUE)
  pth6 <- model$pathfinder(data = data_list, init=0.1, seed = seed,
                          num_paths=10, single_path_draws=40, draws=400,
                          history_size=50, max_lbfgs_iters=100,
                          refresh=0, psis_resample=FALSE)

  fit <- model$sample(data = data_list, seed = seed,
                      chains = 4, iter_warmup = 1000, iter_sampling = 4000, thin = 1, init = pth6)

  return(fit)
}

fit_result <- run_birthday_model("./births_usa_1969.csv", "gpbf6.stan", seed = 1)
res <- bridge_sampler(fit_result, 
                      num_splits = 6, 
                      total_perms = 100, 
                      seed = 1, 
                      return_always = TRUE, 
                      verbose = TRUE, 
                      cores = parallel::detectCores(), 
                      pareto_smoothing_all = FALSE,
                      calculate_covariance = TRUE)[[1]]
numi_split <- I(list(lapply(res, function(x) x$numi)))
deni_split <- I(list(lapply(res, function(x) x$deni)))
pareto_k_numi_10 <- extract_khat(numi_split, n_draws = 10)[[1]]
pareto_k_deni_10 <- extract_khat(deni_split, n_draws = 10)[[1]]
pareto_k_numi_20 <- extract_khat(numi_split, n_draws = 20)[[1]]
pareto_k_deni_20 <- extract_khat(deni_split, n_draws = 20)[[1]]
print(length(pareto_k_numi_10))
print(length(pareto_k_deni_10))
print(length(pareto_k_numi_20))
print(length(pareto_k_deni_20))
print(length(res))
results <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric(), pareto_k_numi10 = numeric(), pareto_k_deni10 = numeric(), pareto_k_numi20 = numeric(), pareto_k_deni20 = numeric())
for (j in 1:length(res)) {
  results <- rbind(results, data.frame(logml = res[[j]]$logml, 
                                       pareto_k_numi = res[[j]]$pareto_k_numi, 
                                       pareto_k_deni = res[[j]]$pareto_k_deni,
                                       mcse_logml = res[[j]]$mcse_logml,
                                        pareto_k_numi10 = pareto_k_numi_10[[j]],
                                        pareto_k_deni10 = pareto_k_deni_10[[j]],
                                        pareto_k_numi20 = pareto_k_numi_20[[j]],
                                        pareto_k_deni20 = pareto_k_deni_20[[j]]
                                       ))
}
write.csv(results, file = "pathfinderbirthdays.csv", row.names = FALSE)



res <- bridge_sampler(fit_result, 
                      num_splits = 6, 
                      total_perms = 100, 
                      seed = 1, 
                      return_always = TRUE, 
                      verbose = TRUE, 
                      cores = parallel::detectCores(), 
                      pareto_smoothing_all = FALSE,
                      calculate_covariance = FALSE)[[1]]
numi_split <- I(list(lapply(res, function(x) x$numi)))
deni_split <- I(list(lapply(res, function(x) x$deni)))
pareto_k_numi_10 <- extract_khat(numi_split, n_draws = 10)[[1]]
pareto_k_deni_10 <- extract_khat(deni_split, n_draws = 10)[[1]]
pareto_k_numi_20 <- extract_khat(numi_split, n_draws = 20)[[1]]
pareto_k_deni_20 <- extract_khat(deni_split, n_draws = 20)[[1]]
print(length(pareto_k_numi_10))
print(length(pareto_k_deni_10))
print(length(pareto_k_numi_20))
print(length(pareto_k_deni_20))
print(length(res))
results <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric(), pareto_k_numi10 = numeric(), pareto_k_deni10 = numeric(), pareto_k_numi20 = numeric(), pareto_k_deni20 = numeric())
for (j in 1:length(res)) {
  results <- rbind(results, data.frame(logml = res[[j]]$logml, 
                                       pareto_k_numi = res[[j]]$pareto_k_numi, 
                                       pareto_k_deni = res[[j]]$pareto_k_deni,
                                       mcse_logml = res[[j]]$mcse_logml,
                                        pareto_k_numi10 = pareto_k_numi_10[[j]],
                                        pareto_k_deni10 = pareto_k_deni_10[[j]],
                                        pareto_k_numi20 = pareto_k_numi_20[[j]],
                                        pareto_k_deni20 = pareto_k_deni_20[[j]]
                                       ))
}
write.csv(results, file = "pathfinderbirthdays_no_cov.csv", row.names = FALSE)



res <- bridge_sampler(fit_result, 
                      num_splits = 6, 
                      total_perms = 100, 
                      seed = 1, 
                      return_always = TRUE, 
                      verbose = TRUE, 
                      cores = parallel::detectCores(), 
                      pareto_smoothing_all = TRUE,
                      calculate_covariance = TRUE)[[1]]
numi_split <- I(list(lapply(res, function(x) x$numi)))
deni_split <- I(list(lapply(res, function(x) x$deni)))
pareto_k_numi_10 <- extract_khat(numi_split, n_draws = 10)[[1]]
pareto_k_deni_10 <- extract_khat(deni_split, n_draws = 10)[[1]]
pareto_k_numi_20 <- extract_khat(numi_split, n_draws = 20)[[1]]
pareto_k_deni_20 <- extract_khat(deni_split, n_draws = 20)[[1]]
print(length(pareto_k_numi_10))
print(length(pareto_k_deni_10))
print(length(pareto_k_numi_20))
print(length(pareto_k_deni_20))
print(length(res))
results <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric(), pareto_k_numi10 = numeric(), pareto_k_deni10 = numeric(), pareto_k_numi20 = numeric(), pareto_k_deni20 = numeric())
for (j in 1:length(res)) {
  results <- rbind(results, data.frame(logml = res[[j]]$logml, 
                                       pareto_k_numi = res[[j]]$pareto_k_numi, 
                                       pareto_k_deni = res[[j]]$pareto_k_deni,
                                       mcse_logml = res[[j]]$mcse_logml,
                                        pareto_k_numi10 = pareto_k_numi_10[[j]],
                                        pareto_k_deni10 = pareto_k_deni_10[[j]],
                                        pareto_k_numi20 = pareto_k_numi_20[[j]],
                                        pareto_k_deni20 = pareto_k_deni_20[[j]]
                                       ))
}
write.csv(results, file = "pathfinderbirthdays_smoothed.csv", row.names = FALSE)



res <- bridge_sampler(fit_result, 
                      num_splits = 6, 
                      total_perms = 100, 
                      seed = 1, 
                      return_always = TRUE, 
                      verbose = TRUE, 
                      cores = parallel::detectCores(), 
                      pareto_smoothing_all = TRUE,
                      calculate_covariance = FALSE)[[1]]
numi_split <- I(list(lapply(res, function(x) x$numi)))
deni_split <- I(list(lapply(res, function(x) x$deni)))
pareto_k_numi_10 <- extract_khat(numi_split, n_draws = 10)[[1]]
pareto_k_deni_10 <- extract_khat(deni_split, n_draws = 10)[[1]]
pareto_k_numi_20 <- extract_khat(numi_split, n_draws = 20)[[1]]
pareto_k_deni_20 <- extract_khat(deni_split, n_draws = 20)[[1]]
print(length(pareto_k_numi_10))
print(length(pareto_k_deni_10))
print(length(pareto_k_numi_20))
print(length(pareto_k_deni_20))
print(length(res))
results <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric(), pareto_k_numi10 = numeric(), pareto_k_deni10 = numeric(), pareto_k_numi20 = numeric(), pareto_k_deni20 = numeric())
for (j in 1:length(res)) {
  results <- rbind(results, data.frame(logml = res[[j]]$logml, 
                                       pareto_k_numi = res[[j]]$pareto_k_numi, 
                                       pareto_k_deni = res[[j]]$pareto_k_deni,
                                       mcse_logml = res[[j]]$mcse_logml,
                                        pareto_k_numi10 = pareto_k_numi_10[[j]],
                                        pareto_k_deni10 = pareto_k_deni_10[[j]],
                                        pareto_k_numi20 = pareto_k_numi_20[[j]],
                                        pareto_k_deni20 = pareto_k_deni_20[[j]]
                                       ))
}
write.csv(results, file = "pathfinderbirthdays_smoothed_no_cov.csv", row.names = FALSE)



results_bruteforce <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
numi <- numeric(100)
deni <- numeric(100)
max_iter <- 100
for(j in 1:max_iter){
  fit_result <- run_birthday_model("./births_usa_1969.csv", "gpbf6.stan", seed = j)
  res <- bridge_sampler(fit_result, 
                        num_splits = 2, 
                        total_perms = 1, 
                        seed = j, 
                        return_always = TRUE, 
                        verbose = TRUE, 
                        cores = parallel::detectCores(), 
                        pareto_smoothing_all = FALSE,
                        calculate_covariance = TRUE)
  numi[[j]] <- res$numi
  deni[[j]] <- res$deni
  results_bruteforce <- rbind(results_bruteforce, data.frame(logml = res$logml, 
                                       pareto_k_numi = res$pareto_k_numi, 
                                       pareto_k_deni = res$pareto_k_deni,
                                       mcse_logml = res$mcse_logml
                                       ))
}
numi <- I(list(numi))
deni <- I(list(deni))
pareto_k_numi_10_brute <- extract_khat(numi, n_draws = 10)[[1]]
pareto_k_deni_10_brute <- extract_khat(deni, n_draws = 10)[[1]]
pareto_k_numi_20_brute <- extract_khat(numi, n_draws = 20)[[1]]
pareto_k_deni_20_brute <- extract_khat(deni, n_draws = 20)[[1]]
results_bruteforce <- cbind(results_bruteforce, pareto_k_numi_10_brute, pareto_k_deni_10_brute, pareto_k_numi_20_brute, pareto_k_deni_20_brute)
write.csv(results_bruteforce, file = "pathfinderbirthdays_bruteforce.csv", row.names = FALSE)



results_bruteforce <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
numi <- numeric(100)
deni <- numeric(100)
max_iter <- 100
for(j in 1:max_iter){
  fit_result <- run_birthday_model("./births_usa_1969.csv", "gpbf6.stan", seed = j)
  res <- bridge_sampler(fit_result, 
                        num_splits = 2, 
                        total_perms = 1, 
                        seed = j, 
                        return_always = TRUE, 
                        verbose = TRUE, 
                        cores = parallel::detectCores(), 
                        pareto_smoothing_all = FALSE,
                        calculate_covariance = FALSE)
  numi[[j]] <- res$numi
  deni[[j]] <- res$deni
  results_bruteforce <- rbind(results_bruteforce, data.frame(logml = res$logml, 
                                       pareto_k_numi = res$pareto_k_numi, 
                                       pareto_k_deni = res$pareto_k_deni,
                                       mcse_logml = res$mcse_logml
                                       ))
}
numi <- I(list(numi))
deni <- I(list(deni))
pareto_k_numi_10_brute <- extract_khat(numi, n_draws = 10)[[1]]
pareto_k_deni_10_brute <- extract_khat(deni, n_draws = 10)[[1]]
pareto_k_numi_20_brute <- extract_khat(numi, n_draws = 20)[[1]]
pareto_k_deni_20_brute <- extract_khat(deni, n_draws = 20)[[1]]
results_bruteforce <- cbind(results_bruteforce, pareto_k_numi_10_brute, pareto_k_deni_10_brute, pareto_k_numi_20_brute, pareto_k_deni_20_brute)
write.csv(results_bruteforce, file = "pathfinderbirthdays_bruteforce_no_cov.csv", row.names = FALSE)



results_bruteforce <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
numi <- numeric(100)
deni <- numeric(100)
max_iter <- 100
for(j in 1:max_iter){
  fit_result <- run_birthday_model("./births_usa_1969.csv", "gpbf6.stan", seed = j)
  res <- bridge_sampler(fit_result, 
                        num_splits = 2, 
                        total_perms = 1, 
                        seed = j, 
                        return_always = TRUE, 
                        verbose = TRUE, 
                        cores = parallel::detectCores(), 
                        pareto_smoothing_all = TRUE,
                        calculate_covariance = TRUE)
  numi[[j]] <- res$numi
  deni[[j]] <- res$deni
  results_bruteforce <- rbind(results_bruteforce, data.frame(logml = res$logml, 
                                       pareto_k_numi = res$pareto_k_numi, 
                                       pareto_k_deni = res$pareto_k_deni,
                                       mcse_logml = res$mcse_logml
                                       ))
}
numi <- I(list(numi))
deni <- I(list(deni))
pareto_k_numi_10_brute <- extract_khat(numi, n_draws = 10)[[1]]
pareto_k_deni_10_brute <- extract_khat(deni, n_draws = 10)[[1]]
pareto_k_numi_20_brute <- extract_khat(numi, n_draws = 20)[[1]]
pareto_k_deni_20_brute <- extract_khat(deni, n_draws = 20)[[1]]
results_bruteforce <- cbind(results_bruteforce, pareto_k_numi_10_brute, pareto_k_deni_10_brute, pareto_k_numi_20_brute, pareto_k_deni_20_brute)
write.csv(results_bruteforce, file = "pathfinderbirthdays_bruteforce_smoothed.csv", row.names = FALSE)



results_bruteforce <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
numi <- numeric(100)
deni <- numeric(100)
max_iter <- 100
for(j in 1:max_iter){
  fit_result <- run_birthday_model("./births_usa_1969.csv", "gpbf6.stan", seed = j)
  res <- bridge_sampler(fit_result, 
                        num_splits = 2, 
                        total_perms = 1, 
                        seed = j, 
                        return_always = TRUE, 
                        verbose = TRUE, 
                        cores = parallel::detectCores(), 
                        pareto_smoothing_all = TRUE,
                        calculate_covariance = FALSE)
  numi[[j]] <- res$numi
  deni[[j]] <- res$deni
  results_bruteforce <- rbind(results_bruteforce, data.frame(logml = res$logml, 
                                       pareto_k_numi = res$pareto_k_numi, 
                                       pareto_k_deni = res$pareto_k_deni,
                                       mcse_logml = res$mcse_logml
                                       ))
}
numi <- I(list(numi))
deni <- I(list(deni))
pareto_k_numi_10_brute <- extract_khat(numi, n_draws = 10)[[1]]
pareto_k_deni_10_brute <- extract_khat(deni, n_draws = 10)[[1]]
pareto_k_numi_20_brute <- extract_khat(numi, n_draws = 20)[[1]]
pareto_k_deni_20_brute <- extract_khat(deni, n_draws = 20)[[1]]
results_bruteforce <- cbind(results_bruteforce, pareto_k_numi_10_brute, pareto_k_deni_10_brute, pareto_k_numi_20_brute, pareto_k_deni_20_brute)
write.csv(results_bruteforce, file = "pathfinderbirthdays_bruteforce_smoothed_no_cov.csv", row.names = FALSE)
