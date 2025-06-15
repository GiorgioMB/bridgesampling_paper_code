### ---- Package-bootstrap helper ---------------------------------------------
devtools::install_github("GiorgioMB/bridgesamplingparetok", force = TRUE)
required_pkgs <- c(
  "rstanarm",
  "bridgesampling",
  "dplyr",
  "tidyr",
  "purrr",
  "future.apply",
  "furrr",
  "glue"
)

# 1. figure out which are missing
to_install <- setdiff(required_pkgs, rownames(installed.packages()))

# 2. install the missing ones (if any), stop if it fails
if (length(to_install)) {
  message("Installing missing packages: ",
          paste(to_install, collapse = ", "))
  tryCatch(
    install.packages(to_install, dependencies = TRUE),
    error = function(e) {
      stop("Package installation failed: ", conditionMessage(e), call. = FALSE)
    }
  )
}

# 3. load everything, stopping if anything still can’t be found
invisible(
  lapply(required_pkgs, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package ", pkg, " is not available even after installation.",
           call. = FALSE)
    }
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
  })
)
### ---------------------------------------------------------------------------

setwd("/scratch/work/micaleg1/bridge_sampling")  
# ---------- Global options ----------------------------------------------------
SEED_BASE              <- 1L                  # reproducibility for data & model draws
N_OBS                  <- 100L                # sample size per simulated data set
K_GRID                 <- 10:104              # range of numbers of covariates
NUM_SPLITS             <- 6L                  # bridge‑sampling split count
TOTAL_PERMS            <- 100L                # permutations in split estimator
N_BRUTEFORCE_ITER      <- 100L                # independent bridges for baseline
RESULTS_DIR            <- "./bridge_results"  # one file per (k, cfg)

# ---- Parallel resource allocation -------------------------------------------
CORES_TOTAL <- parallel::detectCores()
CORES_PER_TASK <- as.integer(Sys.getenv("CORES_PER_TASK", unset = 4))
if (CORES_PER_TASK < 1L || CORES_PER_TASK > CORES_TOTAL)
  stop("CORES_PER_TASK must be between 1 and available cores (", CORES_TOTAL, ")")
MAX_WORKERS <- max(1L, floor(CORES_TOTAL / CORES_PER_TASK))

# future plan: one R process per grid row, each using CORES_PER_TASK internally
future::plan(future::multisession, workers = MAX_WORKERS)
options(mc.cores = CORES_PER_TASK)   # for rstanarm’s within‑process parallelism
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

# ---------- Helper: pull raw vectors & scalars -------------------------------
bridge_to_lists <- function(res) {
  list(
    logml      = res$logml,
    pareto_k_numi = res$pareto_k_numi[[1]]$khat,
    pareto_k_deni = res$pareto_k_deni[[1]]$khat,
    mcse_logml = res$mcse_logml,
    numi       = res$numi,         
    deni       = res$deni          
  )
}

# ---------- Helper: synthetic‑data generator ----------------------------------
make_regression_data <- function(k, n = N_OBS, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- matrix(rnorm(n * k), nrow = n, ncol = k) / sqrt(n * k)
  beta <- rnorm(k); beta[k] <- 0              # last coefficient zero
  y <- as.vector(X %*% beta) + rnorm(n, sd = 2) / sqrt(n)
  list(
    df   = data.frame(                      
             y = y,
             X = I(X)                       
           ),
    beta = beta
  )
}


# ---------- Helper: Stan‑glm fit ---------------------------------------------
fit_horseshoe <- function(data, seed) {
  stan_glm(
    y ~ X,                                   
    data            = data,
    mean_PPD        = FALSE,
    prior           = hs(),
    seed            = seed,
    diagnostic_file = file.path(tempdir(), "diag.csv"),
    refresh         = 0,
    cores           = CORES_PER_TASK
  )
}


# ---------- Helper: split / reshuffled estimator ------------------------------
bridge_split <- function(fit,
                         calculate_covariance = TRUE,
                         pareto_smoothing_all = FALSE,
                         use_ess              = FALSE) {
  bridge_sampler(
    fit,
    method               = "normal",
    num_splits           = NUM_SPLITS,
    total_perms          = TOTAL_PERMS,
    return_always        = TRUE,
    seed                 = SEED_BASE,
    calculate_covariance = calculate_covariance,
    pareto_smoothing_all = pareto_smoothing_all,
    use_ess              = use_ess,
    cores                = CORES_PER_TASK   # keep inner bridges single‑process
  )[[1]]
}

# ---------- Helper: brute‑force estimator -------------------------------------
bridge_bruteforce <- function(data_df,
                              seed_base,
                              calculate_covariance = TRUE,
                              pareto_smoothing_all = FALSE,
                              use_ess              = FALSE) {

  purrr::map(seq_len(N_BRUTEFORCE_ITER), function(i) {

    # ---- (re-)fit Stan ------------------------------------------------------
    fit_i <- fit_horseshoe(data_df, seed = seed_base + i)

    # ---- bridge -------------------------------------------------------------
    bridge_sampler(
      fit_i,
      method               = "normal",
      seed                 = seed_base + i,   # keep RNG streams independent
      calculate_covariance = calculate_covariance,
      pareto_smoothing_all = pareto_smoothing_all,
      use_ess              = use_ess,
      cores                = CORES_PER_TASK   # single-process inside worker
    )
  })
}


# ---------- Helper: tidy output ----------------------------------------------
bridge_list_to_df <- function(lst) {
  purrr::map_dfr(lst, function(res) {
    tibble(
      logml      = res$logml,
      khat_numi  = res$pareto_k_numi[[1]]$khat,
      khat_deni  = res$pareto_k_deni[[1]]$khat,
      mcse_logml = res$mcse_logml
    )
  })
}

# ---------- Run both estimators for one configuration -------------------------
run_one_cfg <- function(k,
                        calculate_covariance,
                        pareto_smoothing_all,
                        use_ess) {

  cfg_label <- glue(
    "cov{calculate_covariance}_smooth{pareto_smoothing_all}_ess{use_ess}"
  )
  message(glue("[k = {k}]  Starting {cfg_label} (cores/task = {CORES_PER_TASK})"))

  # --- simulate + fit --------------------------------------------------------
  dat_obj <- make_regression_data(k, seed = k)
  fit     <- fit_horseshoe(dat_obj$df, seed = k)

  # --- estimators ------------------------------------------------------------
  res_split <- bridge_split(
    fit,
    calculate_covariance = calculate_covariance,
    pareto_smoothing_all = pareto_smoothing_all,
    use_ess              = use_ess
  )

  res_brute <- bridge_bruteforce(
    dat_obj$df,
    seed_base           = k,
    calculate_covariance = calculate_covariance,
    pareto_smoothing_all = pareto_smoothing_all,
    use_ess              = use_ess
  )

  # ---------- NEW: pack exactly what the second script stores ----------------
  # lists of scalars / vectors -----------------------------------------------
  split_lists <- lapply(res_split,  bridge_to_lists)
  brute_lists <- lapply(res_brute,  bridge_to_lists)

  # vectors (length = #splits  or  100) --------------------------------------
  logml_reshuffling         <- lapply(split_lists, `[[`, "logml")
  pareto_k_numi_reshuffling <- lapply(split_lists, `[[`, "pareto_k_numi")
  pareto_k_deni_reshuffling <- lapply(split_lists, `[[`, "pareto_k_deni")
  mcse_logml_reshuffling    <- lapply(split_lists, `[[`, "mcse_logml")
  numi_split                <- lapply(split_lists, `[[`, "numi")
  deni_split                <- lapply(split_lists, `[[`, "deni")

  logml_brute               <- lapply(brute_lists, `[[`, "logml")
  pareto_k_numi_brute       <- lapply(brute_lists, `[[`, "pareto_k_numi")
  pareto_k_deni_brute       <- lapply(brute_lists, `[[`, "pareto_k_deni")
  mcse_logml_brute          <- lapply(brute_lists, `[[`, "mcse_logml")
  numi_brute                <- lapply(brute_lists, `[[`, "numi")
  deni_brute                <- lapply(brute_lists, `[[`, "deni")

  # --- pack + save -----------------------------------------------------------
  out <- list(
    k   = k,
    cfg = list(
      calculate_covariance = calculate_covariance,
      pareto_smoothing_all = pareto_smoothing_all,
      use_ess              = use_ess,
      cores_per_task       = CORES_PER_TASK
    ),

    # --- second-script payload ----------------------------------------------
    logml_reshuffling         = logml_reshuffling,
    pareto_k_numi_reshuffling = pareto_k_numi_reshuffling,
    pareto_k_deni_reshuffling = pareto_k_deni_reshuffling,
    mcse_logml_reshuffling    = mcse_logml_reshuffling,
    logml_brute               = logml_brute,
    pareto_k_numi_brute       = pareto_k_numi_brute,
    pareto_k_deni_brute       = pareto_k_deni_brute,
    mcse_logml_brute          = mcse_logml_brute,
    numi_split                = numi_split,
    deni_split                = deni_split,
    numi_brute                = numi_brute,
    deni_brute                = deni_brute,

    # keep the old tidy data frames too (nice for quick ggplotting)
    split_df = bridge_list_to_df(res_split),
    brute_df = bridge_list_to_df(res_brute)
  )

  out_file <- file.path(RESULTS_DIR, glue("k{k}_", cfg_label, ".rds"))
  saveRDS(out, out_file)
  message(glue("[k = {k}]  Saved -> {out_file}"))
  invisible(out)
}


# ---------- Parameter grid ----------------------------------------------------
param_grid <- expand.grid(k                     = K_GRID,
                          calculate_covariance = c(TRUE, FALSE),
                          pareto_smoothing_all = c(FALSE, TRUE),
                          use_ess              = c(FALSE, TRUE),
                          KEEP.OUT.ATTRS       = FALSE)

# ---------- Parallel execution ------------------------------------------------
furrr::future_pmap(param_grid, ~ run_one_cfg(..1, ..2, ..3, ..4))

message("All configurations complete")
