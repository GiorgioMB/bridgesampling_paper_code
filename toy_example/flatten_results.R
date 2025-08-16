library(purrr)
library(dplyr)
setwd("../toy_example")  

RESULTS_DIR <- "./bridge_results"
OUT_CSV     <- file.path(RESULTS_DIR, "bridges_all_flat.csv")


collect_numeric <- function(x) {
  if (is.null(x)) return(numeric(0))
  if (is.atomic(x) && (is.numeric(x) || is.integer(x) || is.logical(x))) {
    return(as.numeric(x))
  }
  if (is.matrix(x) || is.array(x)) {
    return(as.numeric(x))
  }
  if (is.list(x)) {
    out <- numeric(0)
    for (elt in x) out <- c(out, collect_numeric(elt))
    return(out)
  }
  numeric(0)
}

as_scalar <- function(x) {
  xs <- collect_numeric(x)
  if (!length(xs)) return(NA_real_)
  xs[[1L]]
}

get_ith <- function(lst, i) {
  if (!is.null(lst) && length(lst) >= i) lst[[i]] else NULL
}

assemble_method <- function(method_label,
                            logml_list,
                            k_numi_list,
                            k_deni_list,
                            mcse_list) {
  R <- max(length(logml_list), length(k_numi_list),
           length(k_deni_list), length(mcse_list))
  if (R == 0L) return(NULL)

  rows <- vector("list", R)
  for (i in seq_len(R)) {
    rows[[i]] <- data.frame(
      method         = method_label,
      replicate      = i,
      logml          = as_scalar(get_ith(logml_list, i)),
      pareto_k_numi  = as_scalar(get_ith(k_numi_list, i)),
      pareto_k_deni  = as_scalar(get_ith(k_deni_list, i)),
      mcse_logml     = as_scalar(get_ith(mcse_list, i)),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

flatten_out <- function(out, file = NA_character_) {
  cfg <- out$cfg; if (is.null(cfg)) cfg <- list()

  df_reshuffling <- assemble_method(
    "reshuffling",
    out$logml_reshuffling,
    out$pareto_k_numi_reshuffling,
    out$pareto_k_deni_reshuffling,
    out$mcse_logml_reshuffling
  )

  df_mcmc <- assemble_method(
    "mcmc",
    out$logml_brute,
    out$pareto_k_numi_brute,
    out$pareto_k_deni_brute,
    out$mcse_logml_brute
  )

  parts <- Filter(Negate(is.null), list(df_reshuffling, df_mcmc))
  if (!length(parts)) return(data.frame())

  df <- do.call(rbind, parts)

  # add metadata (leftmost)
  meta <- data.frame(
    file                 = file,
    k                    = if (is.null(out$k)) NA else out$k,
    calculate_covariance = if (is.null(cfg$calculate_covariance)) NA else cfg$calculate_covariance,
    pareto_smoothing_all = if (is.null(cfg$pareto_smoothing_all)) NA else cfg$pareto_smoothing_all,
    use_ess              = if (is.null(cfg$use_ess)) NA else cfg$use_ess,
    stringsAsFactors     = FALSE
  )
  cbind(meta[rep(1L, nrow(df)), , drop = FALSE], df)
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
    cat("Processing file:", f, "\n")
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

bridges_all_flat <- read.csv(OUT_CSV, stringsAsFactors = FALSE)
unique_k_values <- unique(bridges_all_flat$k)
sorted_k_values <- sort(unique_k_values, na.last = TRUE)
if (all(10:104 %in% sorted_k_values)) {
  message("All k values from 10 to 104 are present.")
} else {
  missing_k <- setdiff(10:104, sorted_k_values)
  message("Missing k values: ", paste(missing_k, collapse = ", "))
}
