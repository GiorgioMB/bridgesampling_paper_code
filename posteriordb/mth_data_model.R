##Note: The bridgesampling version of CmdstanR must be installed, comment the line below if already installed
remotes::install_github("stan-dev/cmdstanr@bridge_sampler-method")
Sys.setenv(GITHUB_PAT = "YOUR_TOKEN")
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
pd <- pdb_github()
po <- posterior("Mth_data-Mth_model", pdb = pd)
sc <- stan_code(po)
data <- get_data(po)
model <- stan_model(model_code = sc)
write_stan_file(sc, dir = getwd(), basename = "Mth_data-Mth_model.stan")
model_cmdstanr <- cmdstan_model("Mth_data-Mth_model.stan", force_recompile = TRUE)
init_val <- model_cmdstanr$pathfinder(data = data, 
                                      num_paths = 10, 
                                      single_path_draws = 40, 
                                      draws = 400, 
                                      history_size = 50, 
                                      max_lbfgs_iters = 100, 
                                      psis_resample = FALSE)
fit_stan <- model_cmdstanr$sample(data = data,
                                  chains = 4,
                                  parallel_chains = 4, 
                                  iter_warmup = 1000, 
                                  iter_sampling = 9000, 
                                  thin = 1, 
                                  init = init_val,
                                  seed = 1)
print("Finished fitting the model")
res <- bridge_sampler(fit_stan, 
                      num_splits = 6, 
                      total_perms = 100, 
                      seed = 1, 
                      return_always = TRUE, 
                      verbose = TRUE, 
                      cores = parallel::detectCores(),
                      calculate_covariance = TRUE)[[1]]
print("Finished bridge sampling")
results <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
for (j in 1:length(res)) {
  results <- rbind(results, data.frame(logml = res[[j]]$logml, 
                                       pareto_k_numi = res[[j]]$pareto_k_numi, 
                                       pareto_k_deni = res[[j]]$pareto_k_deni,
                                       mcse_logml = res[[j]]$mcse_logml
                                       ))
}

write.csv(results, file = "mth_data_model_pathfinder.csv", row.names = FALSE)

res <- bridge_sampler(fit_stan, 
                      num_splits = 6, 
                      total_perms = 100, 
                      seed = 1, 
                      return_always = TRUE, 
                      verbose = TRUE, 
                      cores = parallel::detectCores(),
                      calculate_covariance = FALSE)[[1]]
results <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
for (j in 1:length(res)) {
  results <- rbind(results, data.frame(logml = res[[j]]$logml, 
                                       pareto_k_numi = res[[j]]$pareto_k_numi, 
                                       pareto_k_deni = res[[j]]$pareto_k_deni,
                                       mcse_logml = res[[j]]$mcse_logml
                                       ))
}

write.csv(results, file = "mth_data_model_pathfinder_no_cov.csv", row.names = FALSE)



res <- bridge_sampler(fit_stan,
                      num_splits = 6, 
                      total_perms = 100, 
                      seed = 1, 
                      return_always = TRUE, 
                      verbose = TRUE, 
                      cores = parallel::detectCores(), 
                      pareto_smoothing_all = TRUE,
                      calculate_covariance = TRUE)[[1]]
results <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
for (j in 1:length(res)) {
  results <- rbind(results, data.frame(logml = res[[j]]$logml, 
                                       pareto_k_numi = res[[j]]$pareto_k_numi, 
                                       pareto_k_deni = res[[j]]$pareto_k_deni,
                                       mcse_logml = res[[j]]$mcse_logml
                                       ))
}

write.csv(results, file = "mth_data_model_pathfinder_smoothed.csv", row.names = FALSE)

res <- bridge_sampler(fit_stan,
                      num_splits = 6, 
                      total_perms = 100, 
                      seed = 1, 
                      return_always = TRUE, 
                      verbose = TRUE, 
                      cores = parallel::detectCores(), 
                      pareto_smoothing_all = TRUE,
                      calculate_covariance = FALSE)[[1]]
results <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
for (j in 1:length(res)) {
  results <- rbind(results, data.frame(logml = res[[j]]$logml, 
                                       pareto_k_numi = res[[j]]$pareto_k_numi, 
                                       pareto_k_deni = res[[j]]$pareto_k_deni,
                                       mcse_logml = res[[j]]$mcse_logml
                                       ))
}

write.csv(results, file = "mth_data_model_pathfinder_smoothed_no_cov.csv", row.names = FALSE)

attempt_fit <- function(calculate_covariance = TRUE,
                        pareto_smoothing_all = FALSE,
                        sig_figs = 12) {
  repeat {
    try_res <- tryCatch({
      model_cmdstanr <- cmdstan_model(
        "Mth_data-Mth_model.stan",
        force_recompile = TRUE
      )

      init_val <- model_cmdstanr$pathfinder(
        data              = data,
        sig_figs          = sig_figs,
        num_paths         = 10,
        single_path_draws = 40,
        draws             = 400,
        history_size      = 50,
        max_lbfgs_iters   = 100,
        psis_resample     = FALSE
      )

      fit_stan <- model_cmdstanr$sample(
        data             = data,
        chains           = 4,
        parallel_chains  = 4,
        iter_warmup      = 1000,
        iter_sampling    = 9000,
        thin             = 1,
        sig_figs         = sig_figs,
        init             = init_val
      )

      res <- bridge_sampler(
        fit_stan,
        return_always        = TRUE,
        verbose              = TRUE,
        cores                = parallel::detectCores(),
        calculate_covariance = calculate_covariance,
        pareto_smoothing_all = pareto_smoothing_all
      )

      list(res = res)  # success path
    }, error = function(e) {
      message("Attempt failed: ", e$message)
      NULL  
    })

    if (!is.null(try_res)) return(try_res$res)
    message("Retrying...")
  }
}

results_bruteforce <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
numi <- numeric(100)
deni <- numeric(100)
for (i in 1:100) {
  print(paste("Iteration", i))
  res <- attempt_fit(
    calculate_covariance = TRUE,
    pareto_smoothing_all = FALSE,
    sig_figs = 12
  )
  numi[[i]] <- res$numi
  deni[[i]] <- res$deni
  results_bruteforce <- rbind(results_bruteforce, data.frame(logml = res$logml, 
                                       pareto_k_numi = res$pareto_k_numi, 
                                       pareto_k_deni = res$pareto_k_deni,
                                       mcse_logml = res$mcse_logml
                                       ))
}
write.csv(results_bruteforce, file = "mth_data_model_pathfinder_bruteforce.csv", row.names = FALSE)


results_bruteforce <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
numi <- numeric(100)
deni <- numeric(100)
for (i in 1:100) {
  print(paste("Iteration", i))
  res <- attempt_fit(
    calculate_covariance = FALSE,
    pareto_smoothing_all = FALSE,
    sig_figs = 12
  )
  numi[[i]] <- res$numi
  deni[[i]] <- res$deni
  results_bruteforce <- rbind(results_bruteforce, data.frame(logml = res$logml, 
                                       pareto_k_numi = res$pareto_k_numi, 
                                       pareto_k_deni = res$pareto_k_deni,
                                       mcse_logml = res$mcse_logml
                                       ))
}
write.csv(results_bruteforce, file = "mth_data_model_pathfinder_bruteforce_no_cov.csv", row.names = FALSE)


results_bruteforce <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
numi <- numeric(100)
deni <- numeric(100)
for (i in 1:100) {
  print(paste("Iteration", i))
  res <- attempt_fit(
    calculate_covariance = TRUE,
    pareto_smoothing_all = TRUE,
    sig_figs = 12
  )
  numi[[i]] <- res$numi
  deni[[i]] <- res$deni
  results_bruteforce <- rbind(results_bruteforce, data.frame(logml = res$logml, 
                                       pareto_k_numi = res$pareto_k_numi, 
                                       pareto_k_deni = res$pareto_k_deni,
                                       mcse_logml = res$mcse_logml
                                       ))
}
write.csv(results_bruteforce, file = "mth_data_model_pathfinder_bruteforce_smoothed.csv", row.names = FALSE)


results_bruteforce <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
numi <- numeric(100)
deni <- numeric(100)
for (i in 1:100) {
  print(paste("Iteration", i))
  res <- attempt_fit(
    calculate_covariance = FALSE,
    pareto_smoothing_all = TRUE,
    sig_figs = 12
  )
  numi[[i]] <- res$numi
  deni[[i]] <- res$deni
  results_bruteforce <- rbind(results_bruteforce, data.frame(logml = res$logml, 
                                       pareto_k_numi = res$pareto_k_numi, 
                                       pareto_k_deni = res$pareto_k_deni,
                                       mcse_logml = res$mcse_logml
                                       ))
}
write.csv(results_bruteforce, file = "mth_data_model_pathfinder_bruteforce_smoothed_no_cov.csv", row.names = FALSE)
