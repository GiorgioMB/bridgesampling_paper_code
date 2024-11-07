##Note: The bridgesampling version of CmdstanR must be installed, comment the line below if already installed
remotes::install_github("stan-dev/cmdstanr@bridge_sampler-method")
cmdstanr::cmdstan_make_local(cpp_options=list(STAN_THREADS=TRUE),append=TRUE)
cmdstanr::rebuild_cmdstan()
Sys.setenv(GITHUB_PAT = "YOUR_TOKEN")
setwd("../posteriordb/")
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
pd <- pdb_github()
po <- posterior("low_dim_gauss_mix-low_dim_gauss_mix", pdb = pd)
sc <- stan_code(po)
data <- get_data(po)
model <- stan_model(model_code = sc)
write_stan_file(sc, dir = getwd(), basename = "low_dim_gauss_mix.stan")
model_cmdstanr <- cmdstan_model("low_dim_gauss_mix.stan", force_recompile = TRUE)
init_val <- model_cmdstanr$pathfinder(data = data, 
                                      num_paths = 10, 
                                      single_path_draws = 40, 
                                      draws = 400, 
                                      history_size = 50, 
                                      max_lbfgs_iters = 100, 
                                      psis_resample = FALSE)
fit_stan <- model_cmdstanr$sample(data = data,
                                  chains = 4, 
                                  iter_warmup = 1000, 
                                  iter_sampling = 4000, 
                                  thin = 1, 
                                  init = init_val,
                                  seed = 1)
print("Finished fitting the model")
res <- bridge_sampler(fit_stan, 
                      num_splits = 6, total_perms = 100, seed = 1, return_always = TRUE, verbose = TRUE, cores = parallel::detectCores())
split_1 <- res[[1]]
split_2 <- res[[2]]
results <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
for (j in 1:length(split_1)) {
  results <- rbind(results, data.frame(logml = split_1[[j]]$logml, 
                                       pareto_k_numi = split_1[[j]]$pareto_k_numi, 
                                       pareto_k_deni = split_1[[j]]$pareto_k_deni,
                                       mcse_logml = split_1[[j]]$mcse_logml
                                       ))
}

write.csv(results, file = "low_dim_gauss_mix_pathfinder.csv", row.names = FALSE)
write.csv(split_2, file = "low_dim_gauss_mix_pathfinder_split2.csv", row.names = FALSE)

res <- bridge_sampler(fit_stan,
                      num_splits = 6, total_perms = 300, seed = 1, return_always = TRUE, verbose = TRUE, cores = parallel::detectCores(), pareto_smoothing_all = TRUE)
split_1 <- res[[1]]
split_2 <- res[[2]]
results <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
for (j in 1:length(split_1)) {
  results <- rbind(results, data.frame(logml = split_1[[j]]$logml, 
                                       pareto_k_numi = split_1[[j]]$pareto_k_numi, 
                                       pareto_k_deni = split_1[[j]]$pareto_k_deni,
                                       mcse_logml = split_1[[j]]$mcse_logml
                                       ))
}

write.csv(results, file = "low_dim_gauss_mix_pathfinder_smoothed.csv", row.names = FALSE)
write.csv(split_2, file = "low_dim_gauss_mix_pathfinder_split2_smoothed.csv", row.names = FALSE)


res <- bridge_sampler(fit_stan, 
                      method = "student_t", 
                      num_splits = 6, total_perms = 100, seed = 1, return_always = TRUE, verbose = TRUE, cores = parallel::detectCores())
split_1 <- res[[1]]
split_2 <- res[[2]]
results <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
for (j in 1:length(split_1)) {
  results <- rbind(results, data.frame(logml = split_1[[j]]$logml, 
                                       pareto_k_numi = split_1[[j]]$pareto_k_numi, 
                                       pareto_k_deni = split_1[[j]]$pareto_k_deni,
                                       mcse_logml = split_1[[j]]$mcse_logml
                                       ))
}

write.csv(results, file = "low_dim_gauss_mix_t_dist_pathfinder.csv", row.names = FALSE)
write.csv(split_2, file = "low_dim_gauss_mix_t_dist_pathfinder_split2.csv", row.names = FALSE)

res <- bridge_sampler(fit_stan,
                      method = "student_t", 
                      num_splits = 6, total_perms = 300, seed = 1, return_always = TRUE, verbose = TRUE, cores = parallel::detectCores(), pareto_smoothing_all = TRUE)
split_1 <- res[[1]]
split_2 <- res[[2]]
results <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
for (j in 1:length(split_1)) {
  results <- rbind(results, data.frame(logml = split_1[[j]]$logml, 
                                       pareto_k_numi = split_1[[j]]$pareto_k_numi, 
                                       pareto_k_deni = split_1[[j]]$pareto_k_deni,
                                       mcse_logml = split_1[[j]]$mcse_logml
                                       ))
}

write.csv(results, file = "low_dim_gauss_mix_t_dist_pathfinder_smoothed.csv", row.names = FALSE)
write.csv(split_2, file = "low_dim_gauss_mix_t_dist_pathfinder_split2_smoothed.csv", row.names = FALSE)





results_bruteforce <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
for (i in 1:100) {
  iteration_successful <- FALSE
  while(!iteration_successful){
    try({
      print(paste("Iteration", i))
      init_val <- model_cmdstanr$pathfinder(data = data, 
                                      num_paths = 10, 
                                      single_path_draws = 40, 
                                      draws = 400, 
                                      history_size = 50, 
                                      max_lbfgs_iters = 100, 
                                      psis_resample = FALSE)
      fit_stan <- model_cmdstanr$sample(data = data,
                                          chains = 4, 
                                          iter_warmup = 1000, 
                                          iter_sampling = 4000, 
                                          thin = 1, 
                                          init = init_val)
      res <- bridge_sampler(fit_stan, 
                            #method = "student_t", 
                            return_always = TRUE, verbose = TRUE, cores = parallel::detectCores())
      results_bruteforce <- rbind(results_bruteforce, data.frame(logml = res$logml, 
                                          pareto_k_numi = res$pareto_k_numi, 
                                          pareto_k_deni = res$pareto_k_deni,
                                          mcse_logml = res$mcse_logml
                                          ))
      iteration_successful <- TRUE
    })
  }
}
write.csv(results_bruteforce, file = "low_dim_gauss_mix_pathfinder_bruteforce.csv", row.names = FALSE)

results_bruteforce <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
for (i in 1:100) {
  iteration_successful <- FALSE
  while(!iteration_successful){
    try({
      print(paste("Iteration", i))
      init_val <- model_cmdstanr$pathfinder(data = data, 
                                      num_paths = 10, 
                                      single_path_draws = 40, 
                                      draws = 400, 
                                      history_size = 50, 
                                      max_lbfgs_iters = 100, 
                                      psis_resample = FALSE)
      fit_stan <- model_cmdstanr$sample(data = data,
                                          chains = 4, 
                                          iter_warmup = 1000, 
                                          iter_sampling = 4000, 
                                          thin = 1, 
                                          init = init_val)
      res <- bridge_sampler(fit_stan, 
                            #method = "student_t", 
                            return_always = TRUE, verbose = TRUE, cores = parallel::detectCores(), pareto_smoothing_all = TRUE)
      results_bruteforce <- rbind(results_bruteforce, data.frame(logml = res$logml, 
                                          pareto_k_numi = res$pareto_k_numi, 
                                          pareto_k_deni = res$pareto_k_deni,
                                          mcse_logml = res$mcse_logml
                                          ))
      iteration_successful <- TRUE
    })
  }
}
write.csv(results_bruteforce, file = "low_dim_gauss_mix_pathfinder_bruteforce_smoothed.csv", row.names = FALSE)


results_bruteforce <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
for (i in 1:100) {
  iteration_successful <- FALSE
  while(!iteration_successful){
    try({
      print(paste("Iteration", i))
      init_val <- model_cmdstanr$pathfinder(data = data, 
                                      num_paths = 10, 
                                      single_path_draws = 40, 
                                      draws = 400, 
                                      history_size = 50, 
                                      max_lbfgs_iters = 100, 
                                      psis_resample = FALSE)
      fit_stan <- model_cmdstanr$sample(data = data,
                                          chains = 4, 
                                          iter_warmup = 1000, 
                                          iter_sampling = 4000, 
                                          thin = 1, 
                                          init = init_val)
      res <- bridge_sampler(fit_stan, 
                            method = "student_t", 
                            return_always = TRUE, verbose = TRUE, cores = parallel::detectCores())
      results_bruteforce <- rbind(results_bruteforce, data.frame(logml = res$logml, 
                                          pareto_k_numi = res$pareto_k_numi, 
                                          pareto_k_deni = res$pareto_k_deni,
                                          mcse_logml = res$mcse_logml
                                          ))
      iteration_successful <- TRUE
    })
  }
}
write.csv(results_bruteforce, file = "low_dim_gauss_mix_t_dist_pathfinder_bruteforce.csv", row.names = FALSE)

results_bruteforce <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
for (i in 1:100) {
  iteration_successful <- FALSE
  while(!iteration_successful){
    try({
      print(paste("Iteration", i))
      init_val <- model_cmdstanr$pathfinder(data = data, 
                                      num_paths = 10, 
                                      single_path_draws = 40, 
                                      draws = 400, 
                                      history_size = 50, 
                                      max_lbfgs_iters = 100, 
                                      psis_resample = FALSE)
      fit_stan <- model_cmdstanr$sample(data = data,
                                          chains = 4, 
                                          iter_warmup = 1000, 
                                          iter_sampling = 4000, 
                                          thin = 1, 
                                          init = init_val)
      res <- bridge_sampler(fit_stan, 
                            method = "student_t", 
                            return_always = TRUE, verbose = TRUE, cores = parallel::detectCores(), pareto_smoothing_all = TRUE)
      results_bruteforce <- rbind(results_bruteforce, data.frame(logml = res$logml, 
                                          pareto_k_numi = res$pareto_k_numi, 
                                          pareto_k_deni = res$pareto_k_deni,
                                          mcse_logml = res$mcse_logml
                                          ))
      iteration_successful <- TRUE
    })
  }
}
write.csv(results_bruteforce, file = "low_dim_gauss_mix_t_dist_pathfinder_bruteforce_smoothed.csv", row.names = FALSE)
