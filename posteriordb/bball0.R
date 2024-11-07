##Note: The bridgesampling version of CmdstanR must be installed, comment the line below if already installed
remotes::install_github("stan-dev/cmdstanr@bridge_sampler-method")
cmdstanr::cmdstan_make_local(cpp_options=list(STAN_THREADS=TRUE),append=TRUE)
cmdstanr::rebuild_cmdstan()
setwd("../posteriordb/")
Sys.setenv(GITHUB_PAT = "YOUR_TOKEN")
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
po <- posterior("bball_drive_event_0-hmm_drive_0", pdb = pd)
sc <- stan_code(po)
data <- get_data(po)
model <- stan_model(model_code = sc)
write_stan_file(sc, dir = getwd(), basename = "bball0.stan")
model_cmdstanr <- cmdstan_model("bball0.stan", force_recompile = TRUE)
fit_stan <- model_cmdstanr$sample(data = data,
                                  chains = 4, 
                                  iter_warmup = 1000, 
                                  iter_sampling = 4000, 
                                  thin = 1,
                                  sig_figs=9,
                                  seed = 1)
print("Finished fitting the model")
res <- bridge_sampler(fit_stan, num_splits = 6, total_perms = 100, seed = 1, return_always = TRUE, verbose = TRUE, cores = parallel::detectCores())[[1]]
results <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
for (j in 1:length(res)) {
  results <- rbind(results, data.frame(logml = res[[j]]$logml, 
                                       pareto_k_numi = res[[j]]$pareto_k_numi, 
                                       pareto_k_deni = res[[j]]$pareto_k_deni,
                                       mcse_logml = res[[j]]$mcse_logml
                                       ))
}

write.csv(results, file = "bball0_pathfinder.csv", row.names = FALSE)
res <- bridge_sampler(fit_stan, num_splits = 6, total_perms = 100, seed = 1, return_always = TRUE, verbose = TRUE, cores = parallel::detectCores(), pareto_smoothing_all = TRUE)[[1]]
results <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
for (j in 1:length(res)) {
  results <- rbind(results, data.frame(logml = res[[j]]$logml, 
                                       pareto_k_numi = res[[j]]$pareto_k_numi, 
                                       pareto_k_deni = res[[j]]$pareto_k_deni,
                                       mcse_logml = res[[j]]$mcse_logml
                                       ))
}

write.csv(results, file = "bball0_pathfinder_smoothed.csv", row.names = FALSE)

results_bruteforce <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
for (i in 1:100) {
  print(paste("Iteration", i))
  iteration_successful <- FALSE
  while(!iteration_successful){
    try({
      fit_stan <- model_cmdstanr$sample(data = data,
                                        chains = 4, 
                                        iter_warmup = 1000, 
                                        iter_sampling = 4000, 
                                        thin = 1, 
                                        sig_figs=9
                                        )

      res <- bridge_sampler(fit_stan, return_always = TRUE, verbose = TRUE, cores = parallel::detectCores())
      if (is.infinite(res$pareto_k_numi[[1]]$khat)) {
        ##Make it go to the next iteration
        stop("Infinite khat")
      }
      if (is.infinite(res$pareto_k_deni[[1]]$khat)) {
        ##Make it go to the next iteration
        stop("Infinite khat")
      }
      results_bruteforce <- rbind(results_bruteforce, data.frame(logml = res$logml, 
                                       pareto_k_numi = res$pareto_k_numi, 
                                       pareto_k_deni = res$pareto_k_deni,
                                       mcse_logml = res$mcse_logml
                                       ))
      iteration_successful <- TRUE
    }, silent = TRUE)
  }
}
write.csv(results_bruteforce, file = "bball0_pathfinder_bruteforce.csv", row.names = FALSE)



results_bruteforce <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
for (i in 1:100) {
  print(paste("Iteration", i))
  iteration_successful <- FALSE
  while(!iteration_successful){
    try({
      fit_stan <- model_cmdstanr$sample(data = data,
                                        chains = 4, 
                                        iter_warmup = 1000, 
                                        iter_sampling = 4000, 
                                        thin = 1, 
                                        sig_figs=9
                                        )
      res <- bridge_sampler(fit_stan, return_always = TRUE, verbose = TRUE, cores = parallel::detectCores(), pareto_smoothing_all = TRUE)
      if (is.infinite(res$pareto_k_numi[[1]]$khat)) {
        stop("Infinite khat")
      }
      if (is.infinite(res$pareto_k_deni[[1]]$khat)) {
        stop("Infinite khat")
      }
      results_bruteforce <- rbind(results_bruteforce, data.frame(logml = res$logml, 
                                       pareto_k_numi = res$pareto_k_numi, 
                                       pareto_k_deni = res$pareto_k_deni,
                                       mcse_logml = res$mcse_logml
                                       ))
      iteration_successful <- TRUE
    }, silent = TRUE)
  }
}
write.csv(results_bruteforce, file = "bball0_pathfinder_bruteforce_smoothed.csv", row.names = FALSE)
