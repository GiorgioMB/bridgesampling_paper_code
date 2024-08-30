cmdstanr::cmdstan_make_local(cpp_options=list(STAN_THREADS=TRUE),append=TRUE)
cmdstanr::rebuild_cmdstan()
setwd("../posteriordb/")
library(rstan)
library(parallel)
library(foreach)
options(mc.cores = parallel::detectCores())
library(bridgesampling)
rstan_options(auto_write = TRUE)
library(cmdstanr)
library(bayesplot)
## install the beta release version of R package posterior
# install.packages("posterior", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(posteriordb)
library(posterior)
source("./utils/sim_pf.R")
source("./utils/lp_utils.R")
set.seed(123)
pd <- pdb_github()
po <- posterior("nes2000-nes", pdb = pd)
sc <- stan_code(po)
data <- get_data(po)
model <- stan_model(model_code = sc)
write_stan_file(sc, dir = getwd(), basename = "nes2000.stan")
model_cmdstanr <- cmdstan_model("nes2000.stan", force_recompile = TRUE)
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
res <- bridge_sampler(fit_stan, num_splits = 6, total_perms = 100, seed = 1, return_always = TRUE, verbose = TRUE, cores = parallel::detectCores())
results <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
for (j in 1:length(res)) {
  results <- rbind(results, data.frame(logml = res[[j]]$logml, 
                                       pareto_k_numi = res[[j]]$pareto_k_numi, 
                                       pareto_k_deni = res[[j]]$pareto_k_deni,
                                       mcse_logml = res[[j]]$mcse_logml
                                       ))
}

write.csv(results, file = "nes2000_pathfinder.csv", row.names = FALSE)
results_bruteforce <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
for (i in 1:100) {
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
                                    init = init_val,
                                    seed = i)
  res <- bridge_sampler(fit_stan, seed = i, return_always = TRUE, verbose = TRUE, cores = parallel::detectCores())
  results_bruteforce <- rbind(results_bruteforce, data.frame(logml = res$logml, 
                                       pareto_k_numi = res$pareto_k_numi, 
                                       pareto_k_deni = res$pareto_k_deni,
                                       mcse_logml = res$mcse_logml
                                       ))
}
write.csv(results_bruteforce, file = "nes2000_pathfinder_bruteforce.csv", row.names = FALSE)
