##
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
po <- posterior("bball_drive_event_1-hmm_drive_1", pdb = pd)
sc <- stan_code(po)
data <- get_data(po)
model <- stan_model(model_code = sc)
write_stan_file(sc, dir = getwd(), basename = "bball.stan")
model_cmdstanr <- cmdstan_model("bball.stan", force_recompile = TRUE)
init_fun <- function() {
  list(
    theta1 = c(0.5, 0.5),
    theta2 = c(0.5, 0.5),
    phi = c(0, 3),
    lambda = c(0, 3)
  )
}
iteration_successful <- FALSE
while(!iteration_successful){
  tryCatch({
    fit_stan <- model_cmdstanr$sample(data = data,
                                      chains = 4, 
                                      iter_warmup = 1000, 
                                      iter_sampling = 4000, 
                                      thin = 1, 
                                      init = init_fun,
                                      sig_figs=9
                                      )
    print("Finished fitting the model")
    res <- bridge_sampler(fit_stan, num_splits = 6, total_perms = 100, seed = 1, return_always = TRUE, verbose = TRUE, cores = parallel::detectCores())[[1]]
    iteration_successful <- TRUE
  }, error = function(e) {
    print("Iteration failed, retrying")
  })
}
print("Iteration successful")
results <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
for (j in 1:length(res)) {
  results <- rbind(results, data.frame(logml = res[[j]]$logml, 
                                       pareto_k_numi = res[[j]]$pareto_k_numi, 
                                       pareto_k_deni = res[[j]]$pareto_k_deni,
                                       mcse_logml = res[[j]]$mcse_logml
                                       ))
}
write.csv(results, file = "bball_pathfinder.csv", row.names = FALSE)
res <- bridge_sampler(fit_stan, num_splits = 6, total_perms = 100, seed = 1, return_always = TRUE, verbose = TRUE, cores = parallel::detectCores(), pareto_smoothing_all = TRUE)[[1]]
results <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
for (j in 1:length(res)) {
  results <- rbind(results, data.frame(logml = res[[j]]$logml, 
                                       pareto_k_numi = res[[j]]$pareto_k_numi, 
                                       pareto_k_deni = res[[j]]$pareto_k_deni,
                                       mcse_logml = res[[j]]$mcse_logml
                                       ))
}
write.csv(results, file = "bball_pathfinder_smoothed.csv", row.names = FALSE)
results_bruteforce <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
for (i in 1:100) {
  print(paste("Iteration", i))
  iteration_successful <- FALSE
  while(!iteration_successful){
    tryCatch({
      fit_stan <- model_cmdstanr$sample(data = data,
                                        chains = 4, 
                                        iter_warmup = 1000, 
                                        iter_sampling = 4000, 
                                        thin = 1, 
                                        sig_figs=9
                                        )
      res <- bridge_sampler(fit_stan, seed = i, return_always = TRUE, verbose = TRUE, cores = parallel::detectCores())
      iteration_successful <- TRUE
    }, error = function(e) {
      print("Iteration failed, retrying")
    })
  }
  results_bruteforce <- rbind(results_bruteforce, data.frame(logml = res$logml, 
                                       pareto_k_numi = res$pareto_k_numi, 
                                       pareto_k_deni = res$pareto_k_deni,
                                       mcse_logml = res$mcse_logml
                                       ))
}
write.csv(results_bruteforce, file = "bball_pathfinder_bruteforce.csv", row.names = FALSE)
results_bruteforce <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
for (i in 1:100) {
  print(paste("Iteration", i))
  iteration_successful <- FALSE
  while(!iteration_successful){
    tryCatch({
      fit_stan <- model_cmdstanr$sample(data = data,
                                        chains = 4, 
                                        iter_warmup = 1000, 
                                        iter_sampling = 4000, 
                                        thin = 1, 
                                        sig_figs=9
                                        )
      res <- bridge_sampler(fit_stan, seed = i, return_always = TRUE, verbose = TRUE, cores = parallel::detectCores(), pareto_smoothing_all = TRUE)
      iteration_successful <- TRUE
    }, error = function(e) {
      print("Iteration failed, retrying")
    })
  }
  results_bruteforce <- rbind(results_bruteforce, data.frame(logml = res$logml, 
                                       pareto_k_numi = res$pareto_k_numi, 
                                       pareto_k_deni = res$pareto_k_deni,
                                       mcse_logml = res$mcse_logml
                                       ))
}
write.csv(results_bruteforce, file = "bball_pathfinder_bruteforce_smoothed.csv", row.names = FALSE)
