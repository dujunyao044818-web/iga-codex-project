run_external_bulk_ml <- function(cfg) {
  out <- cfg$output_dir
  tab_dir <- file.path(out, "tables")
  ml_dir <- file.path(out, "supplementary", "ml_external_bulk")
  rep_dir <- file.path(out, "reports")
  dir.create(ml_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)

  notes <- character()
  add_note <- function(x) notes <<- c(notes, x)
  safe_csv <- function(x, path) tryCatch(utils::write.csv(x, path, row.names = FALSE), error = function(e) add_note(paste("Could not write", basename(path), conditionMessage(e))))
  safe_ggsave <- function(path, plot, width = 6, height = 5) tryCatch({ ggplot2::ggsave(path, plot, width = width, height = height); TRUE }, error = function(e) { add_note(paste("Skipped figure", basename(path), conditionMessage(e))); FALSE })
  bind_rows_fill <- function(rows) {
    if (length(rows) == 0) return(data.frame())
    all_cols <- unique(unlist(lapply(rows, names)))
    rows <- lapply(rows, function(x) {
      missing <- setdiff(all_cols, names(x))
      for (m in missing) x[[m]] <- NA
      x[, all_cols, drop = FALSE]
    })
    do.call(rbind, rows)
  }
  add_model_row <- function(dataset, model, threshold = "median", n = NA_integer_, n_high = NA_integer_, n_low = NA_integer_, n_features = NA_integer_, status = "skipped", reason = "", auc = NA_real_, accuracy_at_0_5 = NA_real_) {
    data.frame(dataset = dataset, model = model, threshold = threshold, n = n, n_high = n_high, n_low = n_low, n_features = n_features, status = status, reason = reason, auc = auc, accuracy_at_0_5 = accuracy_at_0_5, stringsAsFactors = FALSE)
  }
  auc_rank <- function(y, score) {
    y <- as.integer(y)
    ok <- is.finite(score) & !is.na(y)
    y <- y[ok]; score <- score[ok]
    if (length(unique(y)) != 2) return(NA_real_)
    n1 <- sum(y == 1); n0 <- sum(y == 0)
    if (n1 == 0 || n0 == 0) return(NA_real_)
    r <- rank(score)
    (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  }
  roc_df <- function(y, score) {
    ok <- is.finite(score) & !is.na(y)
    y <- as.integer(y[ok]); score <- score[ok]
    if (length(unique(y)) != 2) return(data.frame())
    ord <- order(score, decreasing = TRUE)
    y <- y[ord]
    tp <- cumsum(y == 1); fp <- cumsum(y == 0)
    data.frame(fpr = c(0, fp / sum(y == 0), 1), tpr = c(0, tp / sum(y == 1), 1))
  }

  score_file <- file.path(tab_dir, "external_validation_signature_scores.csv")
  pathway_file <- file.path(tab_dir, "external_validation_pathway_signature_scores.csv")
  if (!file.exists(score_file) || !file.exists(pathway_file)) {
    msg <- "External bulk ML skipped because external validation score/pathway tables are not available. Run external validation first."
    writeLines(c("# External bulk ML report", "", msg), file.path(rep_dir, "ml_external_bulk_report.md"))
    safe_csv(add_model_row(NA_character_, "all", status = "skipped", reason = msg), file.path(tab_dir, "ml_external_bulk_model_summary.csv"))
    return(invisible(NULL))
  }

  scores <- utils::read.csv(score_file, stringsAsFactors = FALSE, check.names = FALSE)
  pathway <- utils::read.csv(pathway_file, stringsAsFactors = FALSE, check.names = FALSE)
  required_score_cols <- c("dataset", "sample", "inferred_group", "small_cluster_signature_score")
  required_path_cols <- c("dataset", "sample", "signature", "score")
  if (!all(required_score_cols %in% colnames(scores)) || !all(required_path_cols %in% colnames(pathway))) {
    msg <- "External bulk ML skipped because required columns are missing from score/pathway tables."
    writeLines(c("# External bulk ML report", "", msg), file.path(rep_dir, "ml_external_bulk_report.md"))
    safe_csv(add_model_row(NA_character_, "all", status = "skipped", reason = msg), file.path(tab_dir, "ml_external_bulk_model_summary.csv"))
    return(invisible(NULL))
  }

  pathway$key <- paste(pathway$dataset, pathway$sample, sep = "__")
  wide <- tryCatch(reshape2::dcast(pathway, key + dataset + sample ~ signature, value.var = "score", fun.aggregate = mean), error = function(e) data.frame())
  if (nrow(wide) == 0) {
    msg <- "External bulk ML skipped because pathway score matrix could not be reshaped."
    writeLines(c("# External bulk ML report", "", msg), file.path(rep_dir, "ml_external_bulk_report.md"))
    safe_csv(add_model_row(NA_character_, "all", status = "skipped", reason = msg), file.path(tab_dir, "ml_external_bulk_model_summary.csv"))
    return(invisible(NULL))
  }
  scores$key <- paste(scores$dataset, scores$sample, sep = "__")
  dat <- merge(scores, wide, by = c("key", "dataset", "sample"), all = FALSE)
  dat <- dat[dat$inferred_group == "IgAN" & is.finite(dat$small_cluster_signature_score), , drop = FALSE]
  if (nrow(dat) == 0) {
    msg <- "External bulk ML skipped because no finite IgAN external CERI scores were available."
    writeLines(c("# External bulk ML report", "", msg), file.path(rep_dir, "ml_external_bulk_report.md"))
    safe_csv(add_model_row(NA_character_, "all", status = "skipped", reason = msg), file.path(tab_dir, "ml_external_bulk_model_summary.csv"))
    return(invisible(NULL))
  }

  feature_cols <- setdiff(colnames(wide), c("key", "dataset", "sample"))
  model_rows <- list()
  selected_rows <- list()
  pred_rows <- list()

  for (ds in unique(dat$dataset)) {
    d <- dat[dat$dataset == ds, , drop = FALSE]
    d <- d[is.finite(d$small_cluster_signature_score), , drop = FALSE]
    if (nrow(d) < 12) {
      add_note(paste(ds, "ML skipped: fewer than 12 IgAN samples with finite CERI scores."))
      model_rows[[length(model_rows) + 1]] <- add_model_row(ds, "all", n = nrow(d), status = "skipped", reason = "too_few_samples")
      next
    }
    med <- stats::median(d$small_cluster_signature_score, na.rm = TRUE)
    d$label <- ifelse(d$small_cluster_signature_score >= med, 1L, 0L)
    n_pos <- sum(d$label == 1); n_neg <- sum(d$label == 0)
    if (min(n_pos, n_neg) < 5) {
      add_note(paste(ds, "ML skipped: one CERI-high/low class has fewer than 5 samples after median split."))
      model_rows[[length(model_rows) + 1]] <- add_model_row(ds, "all", n = nrow(d), n_high = n_pos, n_low = n_neg, status = "skipped", reason = "class_too_small")
      next
    }
    x <- as.matrix(d[, feature_cols, drop = FALSE])
    x[!is.finite(x)] <- NA_real_
    keep_features <- colSums(!is.na(x)) >= max(4, floor(0.7 * nrow(x)))
    x <- x[, keep_features, drop = FALSE]
    if (ncol(x) < 2) {
      add_note(paste(ds, "ML skipped: fewer than two usable pathway/immune features."))
      model_rows[[length(model_rows) + 1]] <- add_model_row(ds, "all", n = nrow(d), n_high = n_pos, n_low = n_neg, status = "skipped", reason = "too_few_features")
      next
    }
    for (j in seq_len(ncol(x))) {
      nas <- is.na(x[, j]); if (any(nas)) x[nas, j] <- stats::median(x[, j], na.rm = TRUE)
    }
    y <- d$label
    folds <- min(5, min(table(y)))
    if (requireNamespace("glmnet", quietly = TRUE)) {
      for (alpha in c(1, 0.5)) {
        model_name <- ifelse(alpha == 1, "LASSO_glmnet", "ElasticNet_glmnet")
        fit <- tryCatch(glmnet::cv.glmnet(x, y, family = "binomial", alpha = alpha, nfolds = folds, type.measure = "deviance", keep = TRUE), error = function(e) e)
        if (inherits(fit, "error")) {
          add_note(paste(ds, model_name, "failed:", conditionMessage(fit)))
          model_rows[[length(model_rows) + 1]] <- add_model_row(ds, model_name, n = nrow(d), n_high = n_pos, n_low = n_neg, n_features = ncol(x), status = "failed", reason = conditionMessage(fit))
        } else {
          pred <- as.numeric(stats::predict(fit, newx = x, s = "lambda.min", type = "response"))
          auc <- auc_rank(y, pred)
          class <- ifelse(pred >= 0.5, 1L, 0L)
          acc <- mean(class == y)
          co <- as.matrix(stats::coef(fit, s = "lambda.min"))
          sel <- data.frame(dataset = ds, model = model_name, feature = rownames(co), coefficient = as.numeric(co[, 1]), stringsAsFactors = FALSE)
          sel <- sel[sel$feature != "(Intercept)" & sel$coefficient != 0, , drop = FALSE]
          if (nrow(sel) == 0) sel <- data.frame(dataset = ds, model = model_name, feature = "no_nonzero_feature_at_lambda_min", coefficient = 0, stringsAsFactors = FALSE)
          selected_rows[[length(selected_rows) + 1]] <- sel
          model_rows[[length(model_rows) + 1]] <- add_model_row(ds, model_name, n = nrow(d), n_high = n_pos, n_low = n_neg, n_features = ncol(x), status = "success", reason = "exploratory external CERI-high vs CERI-low stratification", auc = auc, accuracy_at_0_5 = acc)
          pred_rows[[length(pred_rows) + 1]] <- data.frame(dataset = ds, model = model_name, sample = d$sample, observed_CERI_high = y, predicted_probability = pred, stringsAsFactors = FALSE)
          r <- roc_df(y, pred)
          if (nrow(r) > 0) {
            p <- ggplot2::ggplot(r, ggplot2::aes(fpr, tpr)) + ggplot2::geom_line() + ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) + ggplot2::coord_equal() + ggplot2::theme_bw(base_size = 10) + ggplot2::labs(title = paste0(ds, " ", model_name, " ROC (exploratory)"), subtitle = paste0("AUC=", round(auc, 3)), x = "False positive rate", y = "True positive rate")
            safe_ggsave(file.path(ml_dir, paste0("Suppl_ML_", ds, "_", model_name, "_ROC.pdf")), p, 5, 5)
          }
        }
      }
    } else {
      add_note(paste(ds, "glmnet unavailable; LASSO and Elastic Net skipped."))
    }
    model_rows[[length(model_rows) + 1]] <- add_model_row(ds, "RandomForest", n = nrow(d), n_high = n_pos, n_low = n_neg, n_features = ncol(x), status = "skipped", reason = "randomForest/ranger dependency not required in GitHub Actions; optional supplementary model not run")
    model_rows[[length(model_rows) + 1]] <- add_model_row(ds, "SVM-RFE", n = nrow(d), n_high = n_pos, n_low = n_neg, n_features = ncol(x), status = "skipped", reason = "caret/e1071 dependency not required in GitHub Actions; optional supplementary model not run")
    model_rows[[length(model_rows) + 1]] <- add_model_row(ds, "XGBoost_or_GBM", n = nrow(d), n_high = n_pos, n_low = n_neg, n_features = ncol(x), status = "skipped", reason = "xgboost/gbm dependency not required in GitHub Actions; optional supplementary model not run")
  }

  model_summary <- if (length(model_rows) > 0) bind_rows_fill(model_rows) else add_model_row(NA_character_, "all", status = "skipped", reason = "No dataset eligible for external bulk ML.")
  selected <- if (length(selected_rows) > 0) bind_rows_fill(selected_rows) else data.frame(dataset = NA_character_, model = NA_character_, feature = NA_character_, coefficient = NA_real_)
  preds <- if (length(pred_rows) > 0) bind_rows_fill(pred_rows) else data.frame()
  safe_csv(model_summary, file.path(tab_dir, "ml_external_bulk_model_summary.csv"))
  safe_csv(selected, file.path(tab_dir, "ml_external_bulk_selected_genes.csv"))
  if (nrow(preds) > 0) safe_csv(preds, file.path(tab_dir, "ml_external_bulk_predictions.csv"))
  safe_csv(if (length(notes) == 0) data.frame(note = "External bulk ML finished without fatal errors; optional unavailable models were recorded as skipped.") else data.frame(note = notes), file.path(tab_dir, "ml_external_bulk_notes.csv"))

  report <- c(
    "# External bulk machine-learning report",
    "",
    "These analyses are exploratory CERI-like stratification models. They should not be described as definitive diagnostic or prognostic models.",
    "",
    "The internal GSE104948 23 vs 4 subtype split is not used as the sole basis for a strong predictive model. External bulk IgAN samples are stratified as CERI-high versus CERI-low using the median CERI score within each usable dataset.",
    "",
    "## Model summary",
    paste(capture.output(print(model_summary)), collapse = "\n"),
    "",
    "## Notes",
    paste0("- ", if (length(notes) == 0) "No fatal ML errors." else notes)
  )
  writeLines(report, file.path(rep_dir, "ml_external_bulk_report.md"))
  invisible(list(summary = model_summary, selected = selected))
}
