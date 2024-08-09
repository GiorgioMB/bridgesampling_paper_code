cmdstanr::cmdstan_make_local(cpp_options=list(STAN_THREADS=TRUE),append=TRUE)
cmdstanr::rebuild_cmdstan()
setwd("/scratch/work/micaleg1/pathfinder")
library(rstan)
library(parallel)
library(foreach)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(cmdstanr)
library(bayesplot)
library(bridgesampling)
Sys.setenv(GITHUB_PAT = "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")

## install the beta release version of R package posterior
# install.packages("posterior", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(posterior)
library(posteriordb)
source("./utils/sim_pf.R")
source("./utils/lp_utils.R")

set.seed(123)

### generate reference posterior samples ###

# compile the model

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
print("Error is not here")
fit_stan <- model$sample(data = data, 
                         chains = 4, 
                         adapt_delta = 0.9,
                         init = init_gp,
                         iter_warmup = 1000, 
                         iter_sampling = 4000, 
                         thin = 1, 
                         seed = 1)
print("Erorr is not here either")
res <- bridgesampling:::bridge_sampler.CmdStanMCMC(fit_stan, num_splits = 6, total_perms = 100, seed = 1, return_always = TRUE, verbose = TRUE)
results <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
for (j in 1:length(res)) {
  results <- rbind(results, data.frame(logml = res[[j]]$logml, 
                                       pareto_k_numi = res[[j]]$pareto_k_numi, 
                                       pareto_k_deni = res[[j]]$pareto_k_deni,
                                       mcse_logml = res[[j]]$mcse_logml
                                       ))
}
write.csv(results, file = "mcycle_gp_accel_gp_pathfinder.csv", row.names = FALSE)

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
  res <- bridgesampling:::bridge_sampler.CmdStanMCMC(fit_stan, seed = i, return_always = TRUE)
  results_bruteforce <- rbind(results_bruteforce, data.frame(logml = res$logml, 
                                       pareto_k_numi = res$pareto_k_numi, 
                                       pareto_k_deni = res$pareto_k_deni,
                                       mcse_logml = res$mcse_logml
                                       ))
}
write.csv(results_bruteforce, file = "mcycle_gp_accel_gp_pathfinder_bruteforce.csv", row.names = FALSE)
