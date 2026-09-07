library(Seurat)        
library(ranger)        
library(glmnet)        
library(xgboost)       
library(dplyr)         
library(caret)         
library(Matrix)        
library(org.Hs.eg.db)  
library(AnnotationDbi)
set.seed(42)           

pipeline_start <- Sys.time()    # start timer
cat("\npipeline started:", format(pipeline_start), "\n\n")

seurat_obj <- readRDS("arterial_seurat_processed.rds")   # load processed object

X_full <- t(GetAssayData(seurat_obj, layer = "data"))     # expression matrix (cells × genes)
y <- seurat_obj$cell_type                                 # cell-type labels

keep <- y != "unknown"                                    # drop unknown cell types
X_full <- X_full[keep, ]
y <- factor(y[keep])

donors_vector <- seurat_obj$donor_id[keep]                # donor IDs for LODO
donors        <- unique(donors_vector)
all_levels    <- levels(y)

cat("Cells after removing unknowns:", nrow(X_full), "\n")
cat("Total genes available:", ncol(X_full), "\n")
cat("Cell types:", paste(all_levels, collapse = ", "), "\n")
cat("LODO folds (donors):", paste(donors, collapse = ", "), "\n\n")

N_VAR_GENES  <- 2000                                      # variance filter per fold
N_TOP_IMPORT <- 480                                       # top genes per fold for frequency tables

y_numeric <- as.integer(y) - 1                            # numeric labels for xgboost

map_ensembl_to_symbol <- function(gene_ids) {             # ensembl -> symbol mapping
  bare_ids <- sub("\\..*$", "", gene_ids)
  symbols  <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys      = bare_ids,
    column    = "SYMBOL",
    keytype   = "ENSEMBL",
    multiVals = "first"
  )
  data.frame(
    ensembl_gene_id = gene_ids,
    hgnc_symbol     = ifelse(is.na(symbols), gene_ids, as.character(symbols)),
    stringsAsFactors = FALSE
  )
}

map_to_symbol_vec <- function(genes) {                    # vector version for RF
  m <- map_ensembl_to_symbol(genes)
  setNames(m$hgnc_symbol, m$ensembl_gene_id)
}

gene_symbol_map <- map_to_symbol_vec(colnames(X_full))    # map all genes once

rf_results       <- list()                                # storage lists
rf_confusion     <- list()
lassoen_results  <- list()
lassoen_confusion <- list()
xgb_results      <- list()
xgb_shap_results <- list()
xgb_confusion    <- list()

for (donor in donors) {                                    # LODO loop
  
  fold_start <- Sys.time()                                # fold timer
  cat("\ndonor:", donor, "- fold started", format(fold_start), "\n")
  
  test_idx  <- donors_vector == donor                     # test split
  train_idx <- donors_vector != donor
  
  X_train_raw <- X_full[train_idx, ]                      # raw matrices
  X_test_raw  <- X_full[test_idx, ]
  y_train     <- y[train_idx]
  y_test      <- y[test_idx]
  y_train_num <- y_numeric[train_idx]
  y_test_num  <- y_numeric[test_idx]
  
  cat("  cells:", nrow(X_train_raw), "train /", nrow(X_test_raw), "test\n")
  
  train_vars <- Matrix::colMeans(X_train_raw^2) - Matrix::colMeans(X_train_raw)^2   # variance
  n_keep     <- min(N_VAR_GENES, ncol(X_train_raw))                                 # top 2000
  top_vars   <- order(train_vars, decreasing = TRUE)[1:n_keep]
  gene_names_fold <- colnames(X_full)[top_vars]
  
  X_train_unscaled <- X_train_raw[, top_vars]              # RF uses unscaled
  X_test_unscaled  <- X_test_raw[, top_vars]
  
  train_mean <- Matrix::colMeans(X_train_unscaled)         # scaling stats
  train_sd   <- sqrt(train_vars[top_vars])
  train_sd[train_sd == 0] <- 1
  
  X_train_scaled <- as.matrix(t((t(X_train_unscaled) - train_mean) / train_sd))     # scaled matrices
  X_test_scaled  <- as.matrix(t((t(X_test_unscaled)  - train_mean) / train_sd))
  colnames(X_train_scaled) <- gene_names_fold
  colnames(X_test_scaled)  <- gene_names_fold
  
  cat("\n  -> Random Forest...\n")                         # RF start
  
  df_train <- data.frame(y = y_train, as.matrix(X_train_unscaled))   # RF df
  
  rf_model <- ranger(                                      # RF model
    dependent.variable.name = "y",
    data        = df_train,
    num.trees   = 500,
    importance  = "permutation",
    scale.permutation.importance = TRUE
  )
  
  rf_pred <- predict(rf_model, data.frame(as.matrix(X_test_unscaled)))$predictions   # RF preds
  
  rf_cm <- confusionMatrix(                               # RF confusion
    factor(rf_pred, levels = all_levels),
    factor(y_test,  levels = all_levels),
    mode = "everything"
  )
  
  rf_balanced_acc <- mean(rf_cm$byClass[, "Sensitivity"], na.rm = TRUE)   # RF metrics
  rf_macro_f1     <- mean(rf_cm$byClass[, "F1"], na.rm = TRUE)
  
  cat("     RF balanced accuracy:", round(rf_balanced_acc, 3),
      "| macro F1:", round(rf_macro_f1, 3), "\n")
  
  rf_confusion[[donor]] <- rf_cm$table                    # save confusion
  write.csv(as.data.frame.matrix(rf_cm$table), paste0("rf_confusion_matrix_", donor, ".csv"))
  
  rf_importance_vec <- rf_model$variable.importance       # RF importance
  rf_genes_df <- data.frame(
    held_out_donor = donor,
    gene           = names(rf_importance_vec),
    gene_symbol    = gene_symbol_map[names(rf_importance_vec)],
    importance     = as.numeric(rf_importance_vec)
  ) %>% arrange(desc(importance))
  
  rf_results[[donor]] <- list(                            # store RF fold
    donor        = donor,
    balanced_acc = rf_balanced_acc,
    macro_f1     = rf_macro_f1,
    genes        = rf_genes_df
  )
  
  rm(df_train, rf_model, rf_pred, rf_cm)                  # cleanup
  gc()
  
  train_donors <- donors_vector[train_idx]                # glmnet fold ids
  fold_ids     <- as.integer(factor(train_donors))
  
  for (alpha_val in c(1, 0.5)) {                          # LASSO + EN loop
    
    model_name <- ifelse(alpha_val == 1, "LASSO", "ElasticNet")   # model label
    cat("\n  ->", model_name, "...\n")
    
    glmnet_model <- cv.glmnet(                            # glmnet model
      X_train_scaled,
      y_train,
      family  = "multinomial",
      alpha   = alpha_val,
      foldid  = fold_ids
    )
    
    glmnet_pred <- predict(glmnet_model, newx = X_test_scaled, s = "lambda.min", type = "class")  # preds
    
    glmnet_cm <- confusionMatrix(                         # confusion
      factor(glmnet_pred, levels = all_levels),
      factor(y_test,      levels = all_levels),
      mode = "everything"
    )
    
    glmnet_balanced_acc <- mean(glmnet_cm$byClass[, "Sensitivity"], na.rm = TRUE)   # metrics
    glmnet_macro_f1     <- mean(glmnet_cm$byClass[, "F1"], na.rm = TRUE)
    
    cat("    ", model_name, "balanced accuracy:", round(glmnet_balanced_acc, 3),
        "| macro F1:", round(glmnet_macro_f1, 3), "\n")
    
    lassoen_confusion[[paste(donor, model_name, sep = "_")]] <- glmnet_cm$table      # save confusion
    write.csv(as.data.frame.matrix(glmnet_cm$table),
              paste0("lodo_confusion_matrix_", model_name, "_", donor, ".csv"))
    
    coefficients <- coef(glmnet_model, s = "lambda.min")   # extract coefs
    coef_gene_names <- rownames(coefficients[[1]])
    
    selected_genes_df <- do.call(rbind, lapply(names(coefficients), function(ct) {   # non-zero genes
      coefs   <- as.numeric(coefficients[[ct]])
      nonzero <- coefs != 0 & coef_gene_names != "(Intercept)"
      if (any(nonzero)) {
        data.frame(
          held_out_donor = donor,
          model          = model_name,
          cell_type      = ct,
          gene           = coef_gene_names[nonzero],
          coefficient    = coefs[nonzero]
        )
      }
    }))
    
    lassoen_results[[paste(donor, model_name, sep = "_")]] <- list(                  # store fold
      donor        = donor,
      model        = model_name,
      balanced_acc = glmnet_balanced_acc,
      macro_f1     = glmnet_macro_f1,
      genes        = selected_genes_df
    )
    
    rm(glmnet_model, glmnet_pred, glmnet_cm)               # cleanup
  }
  gc()
  
  cat("\n  -> XGBoost (gain + SHAP)...\n")                 # xgboost start
  
  dtrain <- xgb.DMatrix(data = X_train_scaled, label = y_train_num)   # train matrix
  dtest  <- xgb.DMatrix(data = X_test_scaled,  label = y_test_num)    # test matrix
  
  xgb_model <- xgb.train(                                  # xgboost model
    params = list(
      objective        = "multi:softprob",
      num_class        = length(all_levels),
      eta              = 0.1,
      max_depth        = 6,
      subsample        = 0.8,
      colsample_bytree = 0.8,
      eval_metric      = "mlogloss"
    ),
    data    = dtrain,
    nrounds = 500,
    evals   = list(train = dtrain),
    verbose = 0
  )
  
  pred_prob <- predict(xgb_model, dtest)                   # preds
  pred_prob <- matrix(pred_prob, nrow = length(y_test_num),
                      ncol = length(all_levels), byrow = FALSE)
  pred_num  <- max.col(pred_prob) - 1
  pred_fac  <- factor(all_levels[pred_num + 1], levels = all_levels)
  
  xgb_cm <- confusionMatrix(                               # confusion
    pred_fac,
    factor(y_test, levels = all_levels),
    mode = "everything"
  )
  
  xgb_balanced_acc <- mean(xgb_cm$byClass[, "Sensitivity"], na.rm = TRUE)   # metrics
  xgb_macro_f1     <- mean(xgb_cm$byClass[, "F1"], na.rm = TRUE)
  
  cat("     XGBoost balanced accuracy:", round(xgb_balanced_acc, 3),
      "| macro F1:", round(xgb_macro_f1, 3), "\n")
  
  xgb_confusion[[donor]] <- xgb_cm$table                   # save confusion
  write.csv(as.data.frame.matrix(xgb_cm$table),
            paste0("xgb_confusion_matrix_", donor, ".csv"))
  
  gain_importance <- xgb.importance(feature_names = gene_names_fold, model = xgb_model) %>%   # gain
    arrange(desc(Gain))
  
  xgb_gain_genes_df <- data.frame(                         # gain df
    held_out_donor = donor,
    model          = "XGBoost",
    cell_type      = "all",
    gene           = gain_importance$Feature,
    gain           = gain_importance$Gain
  )
  
  xgb_results[[donor]] <- list(                            # store fold
    donor        = donor,
    balanced_acc = xgb_balanced_acc,
    macro_f1     = xgb_macro_f1,
    genes        = xgb_gain_genes_df
  )
  
  set.seed(42)                                             # shap subsample
  shap_idx <- sample(nrow(X_train_scaled), min(1000, nrow(X_train_scaled)))
  X_shap   <- X_train_scaled[shap_idx, , drop = FALSE]
  
  shap_contrib <- predict(xgb_model, X_shap, predcontrib = TRUE)      # shap contribs
  shap_array   <- shap_contrib[, , 1:length(gene_names_fold)]
  mean_shap_per_gene <- apply(abs(shap_array), 3, mean)
  
  xgb_shap_genes_df <- data.frame(                         # shap df
    held_out_donor = donor,
    model          = "XGBoost_SHAP",
    cell_type      = "all",
    gene           = gene_names_fold,
    mean_abs_SHAP  = mean_shap_per_gene
  ) %>% arrange(desc(mean_abs_SHAP))
  
  xgb_shap_results[[donor]] <- list(                       # store fold
    donor        = donor,
    balanced_acc = xgb_balanced_acc,
    macro_f1     = xgb_macro_f1,
    genes        = xgb_shap_genes_df
  )
  
  rm(dtrain, dtest, xgb_model, X_shap, shap_contrib, shap_array,      # cleanup
     X_train_unscaled, X_test_unscaled, X_train_scaled, X_test_scaled,
     X_train_raw, X_test_raw)
  gc()
  
  fold_end <- Sys.time()                                   # fold end
  cat("\nDONOR:", donor, "finished (took", format(difftime(fold_end, fold_start)), ")\n")
}

cat("\nall folds done, writing combined outputs\n\n")

all_genes_rf <- do.call(rbind, lapply(rf_results, function(x) x$genes))   # RF combine
write.csv(all_genes_rf, "genes_per_fold_rf.csv", row.names = FALSE)

top_genes_rf_per_fold <- all_genes_rf %>%                                 # RF top480
  group_by(held_out_donor) %>%
  slice_max(order_by = importance, n = N_TOP_IMPORT, with_ties = FALSE) %>%
  ungroup()
write.csv(top_genes_rf_per_fold, "genes_per_fold_rf_top480.csv", row.names = FALSE)

gene_freq_rf <- top_genes_rf_per_fold %>%                                 # RF freq
  group_by(gene, gene_symbol) %>%
  summarise(
    n_folds_selected = n_distinct(held_out_donor),
    mean_importance  = mean(importance),
    .groups = "drop"
  ) %>%
  arrange(desc(n_folds_selected), desc(mean_importance))
write.csv(gene_freq_rf, "gene_frequency_rf.csv", row.names = FALSE)

metrics_rf <- do.call(rbind, lapply(rf_results, function(x) {             # RF metrics
  data.frame(donor = x$donor, balanced_acc = x$balanced_acc, macro_f1 = x$macro_f1)
}))
write.csv(metrics_rf, "rf_lodo_metrics.csv", row.names = FALSE)

rf_results$confusion_matrices <- rf_confusion                            # RF save
saveRDS(rf_results, file = "rf_lodo_results.rds")

cat("Saved RF outputs.\n")

all_genes_lassoen <- do.call(rbind, lapply(lassoen_results, function(x) x$genes))   # glmnet combine

gene_freq_lassoen <- all_genes_lassoen %>%                                # glmnet freq
  group_by(model, gene) %>%
  summarise(
    n_folds_selected = n_distinct(held_out_donor),
    mean_coefficient = mean(abs(coefficient)),
    .groups = "drop"
  ) %>%
  arrange(model, desc(n_folds_selected), desc(mean_coefficient))

performance_lassoen <- do.call(rbind, lapply(lassoen_results, function(x) {         # glmnet metrics
  data.frame(donor = x$donor, model = x$model, balanced_acc = x$balanced_acc, macro_f1 = x$macro_f1)
}))
write.csv(all_genes_lassoen,       "genes_per_fold.csv",            row.names = FALSE)
write.csv(gene_freq_lassoen,       "gene_frequency_lassoen.csv",    row.names = FALSE)
write.csv(performance_lassoen,     "performance_table.csv",         row.names = FALSE)
saveRDS(lassoen_confusion, "confusion_matrices.rds")

lassoen_gene_map <- map_ensembl_to_symbol(unique(all_genes_lassoen$gene))           # glmnet annotate

all_genes_lassoen_annotated <- merge(                                               # annotate full
  all_genes_lassoen, lassoen_gene_map,
  by.x = "gene", by.y = "ensembl_gene_id", all.x = TRUE
)
write.csv(all_genes_lassoen_annotated, "genes_per_fold_annotated.csv", row.names = FALSE)

gene_freq_lassoen_annotated <- merge(                                               # annotate freq
  gene_freq_lassoen, lassoen_gene_map,
  by.x = "gene", by.y = "ensembl_gene_id", all.x = TRUE
)
gene_freq_lassoen_annotated <- gene_freq_lassoen_annotated[
  order(gene_freq_lassoen_annotated$model,
        -gene_freq_lassoen_annotated$n_folds_selected,
        -gene_freq_lassoen_annotated$mean_coefficient),
  c("gene", "hgnc_symbol", "model", "n_folds_selected", "mean_coefficient")
]
write.csv(gene_freq_lassoen_annotated, "gene_frequency_annotated.csv", row.names = FALSE)

n_folds_total <- n_distinct(all_genes_lassoen$held_out_donor)                        # total folds

lasso_all_folds <- gene_freq_lassoen_annotated %>%                                   # robust LASSO
  filter(model == "LASSO", n_f