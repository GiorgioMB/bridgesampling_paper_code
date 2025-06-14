##Note: The bridgesampling version of CmdstanR must be installed, comment the line below if already installed
remotes::install_github("stan-dev/cmdstanr@bridge_sampler-method")
Sys.setenv(GITHUB_PAT = "YOUR_TOKEN")         
cmdstanr::cmdstan_make_local(cpp_options = list(STAN_THREADS = TRUE), append = TRUE)
cmdstanr::rebuild_cmdstan()
setwd("posteriordb")    
suppressPackageStartupMessages({
  library(cmdstanr)
  library(rstan);
  rstan_options(auto_write = TRUE)
  library(bridgesampling)
  library(parallel)
  library(foreach)
  library(bayesplot)
  library(posteriordb)
  library(posterior)
})
options(mc.cores = parallel::detectCores())
source("./utils/sim_pf.R")
source("./utils/lp_utils.R")
set.seed(123)
pd  <- pdb_github()
po  <- posterior("Mtbh_data-Mtbh_model", pdb = pd)
sc  <- stan_code(po)
write_stan_file(sc, dir = getwd(), basename = "Mtbh_data-Mtbh_model.stan")
model_cmdstanr <- cmdstan_model("Mtbh_data-Mtbh_model.stan", force_recompile = TRUE)
data <- get_data(po)


init_val <- model_cmdstanr$pathfinder(
  data               = data,
  num_paths          = 10,
  single_path_draws  = 40,
  draws              = 400,
  sig_figs           = 12,
  history_size       = 50,
  max_lbfgs_iters    = 100,
  psis_resample      = FALSE
)


fit_stan <- model_cmdstanr$sample(
  data             = data,
  chains           = 4,
  parallel_chains  = 4,
  iter_warmup      = 1000,
  iter_sampling    = 4000,
  thin             = 1,
  sig_figs         = 12,
  init             = init_val,
  seed             = 1
)
print("Finished fitting the model")

run_bridge <- function(fit,
                       calculate_covariance = TRUE,
                       pareto_smoothing_all = FALSE,
                       use_ess              = FALSE,
                       file_stub            = "mtbh_data_model_pathfinder") {
  res <- bridge_sampler(
    fit,
    num_splits           = 6,
    total_perms          = 100,
    seed                 = 1,
    return_always        = TRUE,
    verbose              = TRUE,
    use_ess              = use_ess,
    cores                = parallel::detectCores(),
    calculate_covariance = calculate_covariance,
    pareto_smoothing_all = pareto_smoothing_all
  )[[1]]

  numi_split <- I(list(lapply(res, function(x) x$numi)))
  deni_split <- I(list(lapply(res, function(x) x$deni)))
  kh10 <- extract_khat(numi_split, n_draws = 10)[[1]]
  kh10d <- extract_khat(deni_split, n_draws = 10)[[1]]
  kh20 <- extract_khat(numi_split, n_draws = 20)[[1]]
  kh20d <- extract_khat(deni_split, n_draws = 20)[[1]]

  results <- data.frame()
  for (j in seq_along(res)) {
    results <- rbind(results, data.frame(
      logml            = res[[j]]$logml,
      pareto_k_numi    = res[[j]]$pareto_k_numi,
      pareto_k_deni    = res[[j]]$pareto_k_deni,
      mcse_logml       = res[[j]]$mcse_logml,
      pareto_k_numi10  = kh10[[j]],
      pareto_k_deni10  = kh10d[[j]],
      pareto_k_numi20  = kh20[[j]],
      pareto_k_deni20  = kh20d[[j]]
    ))
  }

  suffix <- paste0(
    if (pareto_smoothing_all) "_smoothed" else "",
    if (!calculate_covariance) "_no_cov" else "",
    if (use_ess) "_ess" else "",
    ".csv"
  )
  out_file <- paste0(file_stub, suffix)
  write.csv(results, file = out_file, row.names = FALSE)
  message("Saved results to: ", out_file)
  invisible(results)
}

settings <- expand.grid(calculate_covariance = c(TRUE, FALSE),
                        pareto_smoothing_all = c(FALSE, TRUE),
                        use_ess              = c(FALSE, TRUE),
                        KEEP.OUT.ATTRS       = FALSE)
apply(settings, 1, function(ss) {
  run_bridge(fit_stan,
             calculate_covariance = ss[["calculate_covariance"]],
             pareto_smoothing_all = ss[["pareto_smoothing_all"]],
             use_ess              = ss[["use_ess"]])
})

attempt_fit <- function(calculate_covariance = TRUE,
                        pareto_smoothing_all = FALSE,
                        use_ess              = FALSE,
                        seed                 = 1,
                        sig_figs             = 12) {
  repeat {
    try_res <- tryCatch({
      init_val <- model_cmdstanr$pathfinder(
        data               = data,
        num_paths          = 10,
        single_path_draws  = 40,
        draws              = 400,
        history_size       = 50,
        max_lbfgs_iters    = 100,
        sig_figs           = sig_figs,
        psis_resample      = FALSE
      )
      fit_stan_tmp <- model_cmdstanr$sample(
        data             = data,
        chains           = 4,
        parallel_chains  = 4,
        iter_warmup      = 1000,
        iter_sampling    = 4000,
        thin             = 1,
        sig_figs         = sig_figs,
        init             = init_val
      )
      bridge_sampler(
        fit_stan_tmp,
        seed                 = seed,
        return_always        = TRUE,
        verbose              = TRUE,
        use_ess              = use_ess,
        cores                = parallel::detectCores(),
        calculate_covariance = calculate_covariance,
        pareto_smoothing_all = pareto_smoothing_all
      )
    }, error = function(e) {
      message("Attempt failed: ", e$message)
      NULL
    })
    if (!is.null(try_res)) return(try_res)
    message("Retrying...")
  }
}

run_bruteforce <- function(calculate_covariance = TRUE,
                           pareto_smoothing_all = FALSE,
                           use_ess              = FALSE,
                           file_stub            = "mtbh_data_model_pathfinder",
                           n_bruteforce_iter    = 100,
                           sig_figs             = 12) {

  results_bruteforce <- data.frame()
  numi <- vector("list", n_bruteforce_iter)
  deni <- vector("list", n_bruteforce_iter)

  for (i in seq_len(n_bruteforce_iter)) {
    message(sprintf("Brute-force iteration %d / %d", i, n_bruteforce_iter))

    res <- attempt_fit(
      calculate_covariance = calculate_covariance,
      pareto_smoothing_all = pareto_smoothing_all,
      use_ess              = use_ess,
      seed                 = i,
      sig_figs             = sig_figs
    )

    numi[[i]] <- res$numi
    deni[[i]] <- res$deni
    results_bruteforce <- rbind(
      results_bruteforce,
      data.frame(
        logml         = res$logml,
        pareto_k_numi = res$pareto_k_numi,
        pareto_k_deni = res$pareto_k_deni,
        mcse_logml    = res$mcse_logml
      )
    )
  }

  numi_pack <- I(list(numi))
  deni_pack <- I(list(deni))
  results_bruteforce$pareto_k_numi10 <- extract_khat(numi_pack, n_draws = 10)[[1]]
  results_bruteforce$pareto_k_deni10 <- extract_khat(deni_pack, n_draws = 10)[[1]]
  results_bruteforce$pareto_k_numi20 <- extract_khat(numi_pack, n_draws = 20)[[1]]
  results_bruteforce$pareto_k_deni20 <- extract_khat(deni_pack, n_draws = 20)[[1]]

  suffix <- paste0(
    if (pareto_smoothing_all) "_smoothed" else "",
    if (!calculate_covariance) "_no_cov" else "",
    if (use_ess) "_ess" else "",
    "_bruteforce.csv"
  )
  out_file <- paste0(file_stub, suffix)
  write.csv(results_bruteforce, file = out_file, row.names = FALSE)
  message("Brute-force results saved to: ", out_file)

  invisible(results_bruteforce)
}

settings_bf <- expand.grid(calculate_covariance = c(TRUE, FALSE),
                           pareto_smoothing_all = c(FALSE, TRUE),
                           use_ess              = c(FALSE, TRUE),
                           KEEP.OUT.ATTRS       = FALSE)

apply(settings_bf, 1, function(ss) {
  run_bruteforce(
    calculate_covariance = ss[["calculate_covariance"]],
    pareto_smoothing_all = ss[["pareto_smoothing_all"]],
    use_ess              = ss[["use_ess"]]
  )
})
