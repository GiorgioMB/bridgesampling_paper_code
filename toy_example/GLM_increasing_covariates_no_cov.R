#devtools::install_github("GiorgioMB/bridgesamplingparetok", force = TRUE)
library(rstanarm)
library(tidyr)
library(dplyr)
library(devtools)
library(ggplot2)
library(bridgesampling)
set.seed(1)
n <- 100
results <- data.frame(k=integer(), logml_reshuffling = list(), pareto_k_numi_reshuffling = list(), pareto_k_deni_reshuffling = list(), mcse_logml_reshuffling = list(), logml_brute = list(), pareto_k_numi_brute = list(), pareto_k_deni_brute = list(), mcse_logml_brute = list(), numi_split = list(), deni_split = list(), numi_brute = list(), deni_brute = list())
for (k in 10:103) {
  print(paste("Number of covariates", k, "started"))
  set.seed(k)
  # Generate data
  X <- rnorm(n * k, mean=0, sd=1)
  X <- matrix(data=X, nrow=n, ncol=k) / sqrt(n * k)
  x <- rnorm(k, mean=0, sd=1); x[k] = 0
  y <- as.vector(X %*% x) + rnorm(n, mean=0, sd=2) / sqrt(n)
  num_cols <- k
  df <- data.frame(X=X[, 1:num_cols], y=y)
  # Fit model
  set.seed(k)
  samples <- stan_glm(y ~ X, 
                      data=df, 
                      mean_PPD=FALSE, 
                      prior=hs(), 
                      seed=1, 
                      cores=parallel::detectCores(), 
                      diagnostic_file=file.path(tempdir(), "df2.csv"))
  results_mine <- bridge_sampler(samples, 
                                 method = "normal", 
                                 num_splits = 6, 
                                 total_perms = 100, 
                                 return_always = TRUE, 
                                 pareto_smoothing_all = TRUE, 
                                 cores = parallel::detectCores(),
                                 calculate_covariance = FALSE)[[1]]
  logml_reshuffle <- lapply(results_mine, function(x) x$logml)
  pareto_k_numi_reshuffle <- lapply(results_mine, function(x) x$pareto_k_numi[[1]]$khat)
  pareto_k_deni_reshuffle <- lapply(results_mine, function(x) x$pareto_k_numi[[1]]$khat)
  mcse_logml_reshuffling <- lapply(results_mine, function(x) x$mcse_logml)
  numi_split <- lapply(results_mine, function(x) x$numi)
  deni_split <- lapply(results_mine, function(x) x$deni)

  logml_old <- numeric(100)
  pareto_k_numi_old <- numeric(100)
  pareto_k_deni_old <- numeric(100)
  mcse_logml_old <- numeric(100)
  numi_brute <- numeric(100)
  deni_brute <- numeric(100)

  for (i in 1:100) {
    set.seed(i)
    results_old <- bridge_sampler(samples, 
                                  method = "normal", 
                                  cores=parallel::detectCores(),
                                  calculate_covariance = FALSE)
    logml_old[[i]] <- results_old$logml
    pareto_k_numi_old[[i]] <- results_old$pareto_k_numi[[1]]$khat
    pareto_k_deni_old[[i]] <- results_old$pareto_k_deni[[1]]$khat
    mcse_logml_old[[i]] <- results_old$mcse_logml
    numi_brute[[i]] <- results_old$numi
    deni_brute[[i]] <- results_old$deni

    if(i>1){
      set.seed(k)
      samples <- stan_glm(y ~ X, 
                          data=df, 
                          mean_PPD=FALSE, 
                          prior=hs(), 
                          seed=1, 
                          cores=parallel::detectCores(), 
                          diagnostic_file=file.path(tempdir(), "df2.csv"))
    }
  }
  # Save into results
  new_row <- data.frame(
    k = k,
    logml_reshuffling = I(list(logml_reshuffle)),
    pareto_k_numi_reshuffling = I(list(pareto_k_numi_reshuffle)),
    pareto_k_deni_reshuffling = I(list(pareto_k_deni_reshuffle)),
    mcse_logml_reshuffling = I(list(mcse_logml_reshuffling)),
    logml_brute = I(list(logml_old)),
    pareto_k_numi_brute = I(list(pareto_k_numi_old)),
    pareto_k_deni_brute = I(list(pareto_k_deni_old)),
    mcse_logml_brute = I(list(mcse_logml_old)),
    numi_split = I(list(numi_split)),
    deni_split = I(list(deni_split)),
    numi_brute = I(list(numi_brute)),
    deni_brute = I(list(deni_brute))
  )
  results <- rbind(results, new_row)
  print("Saving progress")
  saveRDS(results, "/scratch/work/micaleg1/bridge_sampling/glm_storing_everything_no_cov.rds")
  print("Progress saved")
  print(paste("Number of covariates", k, "finished"))
}
