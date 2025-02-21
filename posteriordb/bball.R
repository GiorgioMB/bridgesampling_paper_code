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
po <- posterior("bball_drive_event_1-hmm_drive_1", pdb = pd)
sc <- stan_code(po)
data <- get_data(po)
model <- stan_model(model_code = sc)
write_stan_file(sc, dir = getwd(), basename = "bball.stan")
model_cmdstanr <- cmdstan_model("bball.stan", force_recompile = TRUE)
#init_val <- model_cmdstanr$pathfinder(data = data, 
#                                      num_paths = 10, 
#                                      single_path_draws = 40, 
#                                      draws = 400, 
#                                      history_size = 50, 
#                                      max_lbfgs_iters = 100, 
#                                      psis_resample = FALSE)
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
    res <- bridge_sampler(fit_stan, 
                          num_splits = 6, 
                          total_perms = 100, 
                          seed = 1, 
                          return_always = TRUE, 
                          verbose = TRUE, 
                          cores = parallel::detectCores(),
                          calculate_covariance = TRUE)[[1]]
    iteration_successful <- TRUE
  }, error = function(e) {
    print("Iteration failed, retrying")
  })
}
print("Iteration successful")
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
write.csv(results, file = "bball_pathfinder.csv", row.names = FALSE)

res <- bridge_sampler(fit_stan, 
                      num_splits = 6, 
                      total_perms = 100, 
                      seed = 1, 
                      return_always = TRUE, 
                      verbose = TRUE, 
                      cores = parallel::detectCores(),
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
write.csv(results, file = "bball_pathfinder_no_cov.csv", row.names = FALSE)



res <- bridge_sampler(fit_stan, 
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
write.csv(results, file = "bball_pathfinder_smoothed.csv", row.names = FALSE)

res <- bridge_sampler(fit_stan, 
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
write.csv(results, file = "bball_pathfinder_smoothed_no_cov.csv", row.names = FALSE)


results_bruteforce <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
numi <- numeric(100)
deni <- numeric(100)
for (i in 1:100) {
  print(paste("Iteration", i))
  #init_val <- model_cmdstanr$pathfinder(data = data, 
  #                                      num_paths = 10, 
  #                                      single_path_draws = 40, 
  #                                      draws = 400, 
  #                                      history_size = 50, 
  #                                      max_lbfgs_iters = 100, 
  #                                      psis_resample = FALSE)
  iteration_successful <- FALSE
  while(!iteration_successful){
    tryCatch({
      fit_stan <- model_cmdstanr$sample(data = data,
                                        chains = 4, 
                                        iter_warmup = 1000, 
                                        iter_sampling = 4000, 
                                        thin = 1, 
                                        sig_figs=9
                                        #init = init_val,
                                        )
      res <- bridge_sampler(fit_stan, 
                            seed = i, 
                            return_always = TRUE, 
                            verbose = TRUE, 
                            cores = parallel::detectCores(),
                            calculate_covariance = TRUE)
      iteration_successful <- TRUE
    }, error = function(e) {
      print("Iteration failed, retrying")
    })
  }
  numi[[i]] <- res$numi
  deni[[i]] <- res$deni
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
write.csv(results_bruteforce, file = "bball_pathfinder_bruteforce.csv", row.names = FALSE)


results_bruteforce <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
numi <- numeric(100)
deni <- numeric(100)
for (i in 1:100) {
  print(paste("Iteration", i))
  #init_val <- model_cmdstanr$pathfinder(data = data, 
  #                                      num_paths = 10, 
  #                                      single_path_draws = 40, 
  #                                      draws = 400, 
  #                                      history_size = 50, 
  #                                      max_lbfgs_iters = 100, 
  #                                      psis_resample = FALSE)
  iteration_successful <- FALSE
  while(!iteration_successful){
    tryCatch({
      fit_stan <- model_cmdstanr$sample(data = data,
                                        chains = 4, 
                                        iter_warmup = 1000, 
                                        iter_sampling = 4000, 
                                        thin = 1, 
                                        sig_figs=9
                                        #init = init_val,
                                        )
      res <- bridge_sampler(fit_stan, 
                            seed = i, 
                            return_always = TRUE, 
                            verbose = TRUE, 
                            cores = parallel::detectCores(),
                            calculate_covariance = FALSE)
      iteration_successful <- TRUE
    }, error = function(e) {
      print("Iteration failed, retrying")
    })
  }
  numi[[i]] <- res$numi
  deni[[i]] <- res$deni
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
write.csv(results_bruteforce, file = "bball_pathfinder_bruteforce_no_cov.csv", row.names = FALSE)


results_bruteforce <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
numi <- numeric(100)
deni <- numeric(100)
for (i in 1:100) {
  print(paste("Iteration", i))
  #init_val <- model_cmdstanr$pathfinder(data = data, 
  #                                      num_paths = 10, 
  #                                      single_path_draws = 40, 
  #                                      draws = 400, 
  #                                      history_size = 50, 
  #                                      max_lbfgs_iters = 100, 
  #                                      psis_resample = FALSE)
  iteration_successful <- FALSE
  while(!iteration_successful){
    tryCatch({
      fit_stan <- model_cmdstanr$sample(data = data,
                                        chains = 4, 
                                        iter_warmup = 1000, 
                                        iter_sampling = 4000, 
                                        thin = 1, 
                                        sig_figs=9
                                        #init = init_val,
                                        )
      res <- bridge_sampler(fit_stan, 
                            seed = i, 
                            return_always = TRUE, 
                            verbose = TRUE, 
                            cores = parallel::detectCores(),
                            pareto_smoothing_all = TRUE,
                            calculate_covariance = TRUE)
      iteration_successful <- TRUE
    }, error = function(e) {
      print("Iteration failed, retrying")
    })
  }
  numi[[i]] <- res$numi
  deni[[i]] <- res$deni
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
write.csv(results_bruteforce, file = "bball_pathfinder_bruteforce_smoothed.csv", row.names = FALSE)


results_bruteforce <- data.frame(logml = numeric(), pareto_k_numi = numeric(), pareto_k_deni = numeric(), mcse_logml = numeric())
numi <- numeric(100)
deni <- numeric(100)
for (i in 1:100) {
  print(paste("Iteration", i))
  #init_val <- model_cmdstanr$pathfinder(data = data, 
  #                                      num_paths = 10, 
  #                                      single_path_draws = 40, 
  #                                      draws = 400, 
  #                                      history_size = 50, 
  #                                      max_lbfgs_iters = 100, 
  #                                      psis_resample = FALSE)
  iteration_successful <- FALSE
  while(!iteration_successful){
    tryCatch({
      fit_stan <- model_cmdstanr$sample(data = data,
                                        chains = 4, 
                                        iter_warmup = 1000, 
                                        iter_sampling = 4000, 
                                        thin = 1, 
                                        sig_figs=9
                                        #init = init_val,
                                        )
      res <- bridge_sampler(fit_stan, 
                            seed = i, 
                            return_always = TRUE, 
                            verbose = TRUE, 
                            cores = parallel::detectCores(),
                            pareto_smoothing_all = TRUE,
                            calculate_covariance = FALSE)
      iteration_successful <- TRUE
    }, error = function(e) {
      print("Iteration failed, retrying")
    })
  }
  numi[[i]] <- res$numi
  deni[[i]] <- res$deni
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
write.csv(results_bruteforce, file = "bball_pathfinder_bruteforce_smoothed_no_cov.csv", row.names = FALSE)
