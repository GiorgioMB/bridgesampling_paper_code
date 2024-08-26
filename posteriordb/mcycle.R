##Note: The bridgesampling version of CmdstanR must be installed

cmdstanr::cmdstan_make_local(cpp_options=list(STAN_THREADS=TRUE),append=TRUE)
cmdstanr::rebuild_cmdstan()
setwd("../posteriordb/")
library(rstan)
library(parallel)
library(foreach)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(cmdstanr)
library(bayesplot)
library(bridgesampling)
## install the beta release version of R package posterior
# install.packages("posterior", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(posterior)
library(posteriordb)
source("./utils/sim_pf.R")
source("./utils/lp_utils.R")

set.seed(123)
pd <- pdb_github()
po <- posterior("mcycle_gp-accel_gp", pdb = pd)
data <- get_data(po)
model <- cmdstan_model("mcycle_gp_accel_gp.stan", force_recompile = TRUE)
init_gp <- model$pathfinder(data = data, 
                            num_paths = 10, 
                            single_path_draws = 40, 
                            draws = 400, 
                            history_size = 50, 
                            max_lbfgs_iters = 100, 
                            psis_resample = FALSE)
fit_stan <- model$sample(data = data, 
                         chains = 4, 
                         adapt_delta = 0.9,
                         init = init_gp,
                         iter_warmup = 1000, 
                         iter_sampling = 4000, 
                         thin = 1, 
                         seed = 1)
res <- bridge_sampler(fit_stan, num_splits = 6, total_perms = 100, seed = 1, return_always = TRUE, verbose = TRUE)
results <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
for (j in 1:length(res)) {
  results <- rbind(results, data.frame(logml = res[[j]]$logml, 
                                       pareto_k_numi = res[[j]]$pareto_k_numi, 
                                       pareto_k_deni = res[[j]]$pareto_k_deni,
                                       mcse_logml = res[[j]]$mcse_logml
                                       ))
}
write.csv(results, file = "mcycle_gp_accel_gp_pathfinder.csv", row.names = FALSE)
res <- bridge_sampler(fit_stan, num_splits = 6, total_perms = 100, seed = 1, return_always = TRUE, verbose = TRUE, pareto_smoothing_all = TRUE, cores = parallel::detectCores())
results <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
for (j in 1:length(res)) {
  results <- rbind(results, data.frame(logml = res[[j]]$logml, 
                                       pareto_k_numi = res[[j]]$pareto_k_numi, 
                                       pareto_k_deni = res[[j]]$pareto_k_deni,
                                       mcse_logml = res[[j]]$mcse_logml
                                       ))
}
write.csv(results, file = "mcycle_gp_accel_gp_pathfinder_smoothed.csv", row.names = FALSE)
results_bruteforce <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
for (i in 1:100) {
  print(paste("Iteration", i))
  init_gp <- model$pathfinder(data = data, 
                              num_paths = 10, 
                              single_path_draws = 40, 
                              draws = 400, 
                              history_size = 50, 
                              max_lbfgs_iters = 100, 
                              psis_resample = FALSE)
  fit_stan <- model$sample(data = data, 
                           chains = 4, 
                           adapt_delta = 0.9,
                           init = init_gp,
                           iter_warmup = 1000, 
                           iter_sampling = 4000, 
                           thin = 1, 
                           seed = i)
  res <- bridge_sampler(fit_stan, seed = i, return_always = TRUE)
  results_bruteforce <- rbind(results_bruteforce, data.frame(logml = res$logml, 
                                       pareto_k_numi = res$pareto_k_numi, 
                                       pareto_k_deni = res$pareto_k_deni,
                                       mcse_logml = res$mcse_logml
                                       ))
}
write.csv(results_bruteforce, file = "mcycle_gp_accel_gp_pathfinder_bruteforce.csv", row.names = FALSE)
results_bruteforce <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
for (i in 1:100) {
  iteration_successful <- FALSE
  while(!iteration_successful){
    tryCatch({
      print(paste("Iteration", i))
      model <- cmdstan_model("mcycle_gp_accel_gp.stan", force_recompile = TRUE)
      init_gp <- model$pathfinder(data = data, 
                                  num_paths = 10, 
                                  single_path_draws = 40, 
                                  draws = 400, 
                                  history_size = 50, 
                                  max_lbfgs_iters = 100, 
                                  psis_resample = FALSE)
      fit_stan <- model$sample(data = data, 
                              chains = 4, 
                              adapt_delta = 0.9,
                              init = init_gp,
                              iter_warmup = 1000, 
                              iter_sampling = 4000, 
                              thin = 1, 
                              seed = i)
      res <- bridge_sampler(fit_stan, seed = i, return_always = TRUE, verbose = TRUE, cores = parallel::detectCores(), pareto_smoothing_all = TRUE)
      if(any(is.na(res$pareto_k_deni[[1]]$khat)) | any(is.infinite(res$pareto_k_deni[[1]]$khat))){
        iteration_successful <- FALSE
      } else {
        if(any(is.na(res$pareto_k_numi[[1]]$khat)) | any(is.infinite(res$pareto_k_numi[[1]]$khat))){
          iteration_successful <- FALSE
        } else {
          iteration_successful <- TRUE
          results_bruteforce <- rbind(results_bruteforce, data.frame(logml = res$logml, 
                                          pareto_k_numi = res$pareto_k_numi, 
                                          pareto_k_deni = res$pareto_k_deni,
                                          mcse_logml = res$mcse_logml
                                          ))
        }
      }
    })
  }
}
write.csv(results_bruteforce, file = "mcycle_gp_accel_gp_pathfinder_bruteforce_smoothed.csv", row.names = FALSE)
