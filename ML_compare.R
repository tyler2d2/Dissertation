library(dplyr) 
library(tidyr)  

N_TOP_IMPORT <- 480 

rf_folds    <- read.csv("genes_per_fold_rf.csv")            # RF fold-level genes
lasso_folds <- read.csv("genes_per_fold_annotated.csv")     # LASSO/EN fold-level genes
xgb_folds   <- read.csv("genes_per_fold_xgb.csv")           # XGB gain fold-level genes
shap_folds  <- read.csv("genes_per_fold_xgb_shap.csv")      # XGB SHAP fold-level genes

rf_symbol_col <- intersect(c("gene_symbol", "symbol", "hgnc_symbol"), names(rf_folds))   # symbol col
if (length(rf_symbol_col) == 0) stop("no symbol column found")
rf_symbol_col <- rf_symbol_col[1]

n_folds_total <- n_distinct(c(rf_folds$held_out_donor,      # total folds
                              lasso_folds$held_out_donor,
                              xgb_folds$held_out_donor,
                              shap_folds$held_out_donor))

cat("loaded\n")

rf_std <- rf_folds %>%                                      # RF top480 per fold
  group_by(held_out_donor) %>%
  slice_max(order_by = importance, n = N_TOP_IMPORT, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(gene, symbol = .data[[rf_symbol_col]], model = "RandomForest",
            held_out_donor, importance = importance)

lasso_std <- lasso_folds %>%                                # LASSO/EN top480 per fold
  group_by(gene, hgnc_symbol, model, held_out_donor) %>%
  summarise(importance = max(abs(coefficient)), .groups = "drop") %>%
  group_by(model, held_out_donor) %>%
  slice_max(order_by = importance, n = N_TOP_IMPORT, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(gene, symbol = hgnc_symbol, model, held_out_donor, importance)

xgb_gain_ranked <- xgb_folds %>%                             # XGB gain ranks
  group_by(held_out_donor) %>%
  mutate(gain_rank = rank(-gain, ties.method = "min")) %>%
  ungroup() %>%
  transmute(gene, symbol = hgnc_symbol, held_out_donor, gain_rank)

shap_ranked <- shap_folds %>%                                # XGB SHAP ranks
  group_by(held_out_donor) %>%
  mutate(shap_rank = rank(-mean_abs_SHAP, ties.method = "min")) %>%
  ungroup() %>%
  transmute(gene, symbol = hgnc_symbol, held_out_donor, shap_rank)

xgb_combined_std <- full_join(                               # combine gain + SHAP
  xgb_gain_ranked, shap_ranked,
  by = c("gene", "symbol", "held_out_donor")
) %>%
  rowwise() %>%
  mutate(combined_rank = mean(c(gain_rank, shap_rank), na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(held_out_donor) %>%
  slice_min(order_by = combined_rank, n = N_TOP_IMPORT, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(gene, symbol, model = "XGBoost_combined", held_out_donor,
            importance = -combined_rank)

all_models_long <- bind_rows(rf_std, lasso_std, xgb_combined_std)   # merge all models

stopifnot(all(c("RandomForest", "LASSO", "ElasticNet", "XGBoost_combined") %in%
                unique(all_models_long$model)))                     # sanity check

cat("standardised\n")

per_model_gene <- all_models_long %>%                               # per-gene stats
  group_by(model, held_out_donor) %>%
  mutate(fold_rank = if (all(is.na(importance))) NA_real_ else rank(-importance, ties.method = "min")) %>%
  ungroup() %>%
  group_by(model, gene, symbol) %>%
  summarise(
    n_folds_selected = n_distinct(held_out_donor),
    mean_importance   = mean(importance, na.rm = TRUE),
    mean_fold_rank    = mean(fold_rank, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(per_model_gene, "per_model_gene_summary.csv", row.names = FALSE)   # save summary

consensus_wide <- per_model_gene %>%                                       # wide format
  pivot_wider(
    id_cols     = gene,
    names_from  = model,
    values_from = c(n_folds_selected, mean_importance, mean_fold_rank),
    names_glue  = "{model}_{.value}"
  )

gene_symbols <- all_models_long %>%                                        # symbol lookup
  distinct(gene, symbol) %>%
  dplyr::filter(!is.na(symbol)) %>%
  group_by(gene) %>%
  summarise(symbol = dplyr::first(symbol), .groups = "drop")

consensus_wide <- consensus_wide %>%                                       # add symbols
  left_join(gene_symbols, by = "gene") %>%
  relocate(symbol, .after = gene)

folds_cols <- grep("_n_folds_selected$", names(consensus_wide), value = TRUE)   # fill NAs
consensus_wide[folds_cols] <- lapply(consensus_wide[folds_cols], function(x) ifelse(is.na(x), 0, x))

consensus_model_names <- c("RandomForest", "LASSO", "ElasticNet", "XGBoost_combined")   # model names
n_consensus_models    <- length(consensus_model_names)
consensus_folds_cols  <- paste0(consensus_model_names, "_n_folds_selected")

detail_model_names <- c("RandomForest", "LASSO", "ElasticNet", "XGBoost_combined")      # rank cols

consensus_wide <- consensus_wide %>%                                       # consensus stats
  rowwise() %>%
  mutate(
    n_models_selected_all_folds   = sum(c_across(all_of(consensus_folds_cols)) == n_folds_total),
    n_models_selected_5plus_folds = sum(c_across(all_of(consensus_folds_cols)) >= 5),
    n_models_selected_any_fold    = sum(c_across(all_of(consensus_folds_cols)) > 0),
    models_selected_all_folds     = paste(
      consensus_model_names[which(c_across(all_of(consensus_folds_cols)) == n_folds_total)],
      collapse = ", "
    ),
    mean_rank_across_models = mean(
      c_across(all_of(paste0(detail_model_names, "_mean_fold_rank"))),
      na.rm = TRUE
    )
  ) %>%
  ungroup() %>%
  arrange(desc(n_models_selected_all_folds), desc(n_models_selected_5plus_folds), mean_rank_across_models)

write.csv(consensus_wide, "consensus_gene_table.csv", row.names = FALSE)   # save consensus

cat("consensus done\n")

selected_all_folds_by_model <- lapply(consensus_model_names, function(m) {   # per-model robust genes
  col <- paste0(m, "_n_folds_selected")
  consensus_wide$gene[!is.na(consensus_wide[[col]]) & consensus_wide[[col]] == n_folds_total]
})
names(selected_all_folds_by_model) <- consensus_model_names

overlap_matrix <- matrix(                                                   # overlap matrix
  NA_integer_, nrow = n_consensus_models, ncol = n_consensus_models,
  dimnames = list(consensus_model_names, consensus_model_names)
)
for (i in seq_len(n_consensus_models)) {
  for (j in seq_len(n_consensus_models)) {
    overlap_matrix[i, j] <- length(intersect(
      selected_all_folds_by_model[[i]],
      selected_all_folds_by_model[[j]]
    ))
  }
}

write.csv(overlap_matrix, "model_overlap_matrix.csv")                       # save overlap

panel_480 <- consensus_wide %>%                                             # top 480 consensus
  arrange(mean_rank_across_models) %>%
  slice_head(n = 480)

write.csv(panel_480, "consensus_panel_480.csv", row.names = FALSE)          # save panel

cat("480-gene panel saved\n")
cat("saved\n")
