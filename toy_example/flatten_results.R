library(purrr)
library(dplyr)
setwd("../toy_example")  

RESULTS_DIR <- "./bridge_results"
OUT_CSV     <- file.path(RESULTS_DIR, "bridges_all_flat.csv")

nz <- function(x) if (is.null(x)) NA else x

# helper: flatten into a long tibble
flatten_component <- function(lst, estimator, measure) {
  if (is.null(lst)) return(NULL)
  purrr::imap_dfr(lst, function(x, idx) {
    xv <- as.numeric(unlist(x, recursive = TRUE, use.names = FALSE))
    if (!length(xv)) return(NULL)
    if (length(xv) == 1L) {
      tibble::tibble(
        estimator = estimator,
        replicate = as.integer(idx),
        measure   = measure,
        position  = NA_integer_,
        value     = xv[1]
      )
    } else {
      tibble::tibble(
        estimator = estimator,
        replicate = as.integer(idx),
        measure   = measure,
        position  = seq_along(xv),
        value     = xv
      )
    }
  })
}

# flatten one .rds object to rows
flatten_out <- function(out, file = NA_character_) {
  cfg <- out$cfg; if (is.null(cfg)) cfg <- list()
  parts <- list(
    flatten_component(out$logml_reshuffling,         "split", "logml"),
    flatten_component(out$pareto_k_numi_reshuffling, "split", "pareto_k_numi"),
    flatten_component(out$pareto_k_deni_reshuffling, "split", "pareto_k_deni"),
    flatten_component(out$mcse_logml_reshuffling,    "split", "mcse_logml"),
    flatten_component(out$numi_split,                "split", "numi"),
    flatten_component(out$deni_split,                "split", "deni"),
    flatten_component(out$logml_brute,               "brute", "logml"),
    flatten_component(out$pareto_k_numi_brute,       "brute", "pareto_k_numi"),
    flatten_component(out$pareto_k_deni_brute,       "brute", "pareto_k_deni"),
    flatten_component(out$mcse_logml_brute,          "brute", "mcse_logml"),
    flatten_component(out$numi_brute,                "brute", "numi"),
    flatten_component(out$deni_brute,                "brute", "deni")
  )
  df <- dplyr::bind_rows(parts)
  if (!nrow(df)) return(df)
  dplyr::mutate(
    df,
    k                     = nz(out$k),
    calculate_covariance  = nz(cfg$calculate_covariance),
    pareto_smoothing_all  = nz(cfg$pareto_smoothing_all),
    use_ess               = nz(cfg$use_ess),
    cores_per_task        = nz(cfg$cores_per_task),
    file                  = file,
    .before = 1
  )
}

choose_writer <- function() {
  if (requireNamespace("data.table", quietly = TRUE)) {
    function(df, path, append) data.table::fwrite(df, path, append = append)
  } else if (requireNamespace("readr", quietly = TRUE)) {
    function(df, path, append) readr::write_csv(df, path, append = append)
  } else {
    function(df, path, append)
      utils::write.table(
        df, file = path, sep = ",",
        row.names = FALSE, col.names = !append, append = append, qmethod = "double"
      )
  }
}

write_long_csv <- function(results_dir = RESULTS_DIR, out_csv = OUT_CSV) {
  files <- list.files(results_dir, pattern = "\\.rds$", full.names = TRUE)
  if (!length(files)) stop("No .rds files found in: ", results_dir)

  writer   <- choose_writer()
  appended <- FALSE
  if (file.exists(out_csv)) file.remove(out_csv)

  for (f in files) {
    out <- tryCatch(readRDS(f), error = function(e) {
      message("Skipping ", basename(f), " (readRDS failed): ", conditionMessage(e))
      NULL
    })
    if (is.null(out)) next

    df <- tryCatch(flatten_out(out, file = basename(f)), error = function(e) {
      message("Skipping ", basename(f), " (flatten failed): ", conditionMessage(e))
      NULL
    })
    if (!is.null(df) && nrow(df)) {
      writer(df, out_csv, append = appended)
      appended <- TRUE
    }
    rm(out, df); gc()
  }

  message("Wrote: ", out_csv)
  invisible(out_csv)
}

write_long_csv()
