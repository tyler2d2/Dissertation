library(dplyr)      
library(ggplot2)  
library(tidyr)      
library(stringr)    
library(forcats)    

donor_labels <- c(patient1 = "Donor 1", patient2 = "Donor 2", patient3 = "Donor 3",
                  patient4 = "Donor 4", patient5 = "Donor 5", patient6 = "Donor 6")   # nicer labels

dir.create("figures", showWarnings = FALSE)               # ensure output dir exists

rf_metrics      <- read.csv("rf_lodo_metrics.csv")       %>% mutate(model = "Random Forest")   # RF metrics
xgb_metrics     <- read.csv("performance_table_xgb.csv") %>% mutate(model = "XGBoost")         # XGB metrics
lassoen_metrics <- read.csv("performance_table.csv")                                          # LASSO/EN metrics

perf_by_donor <- bind_rows(                              # combine metrics
  rf_metrics      %>% dplyr::select(donor, model, balanced_acc, macro_f1),
  lassoen_metrics %>% dplyr::select(donor, model, balanced_acc, macro_f1),
  xgb_metrics     %>% dplyr::select(donor, model, balanced_acc, macro_f1)
) %>%
  pivot_longer(cols = c(balanced_acc, macro_f1), names_to = "metric", values_to = "value") %>%   # long format
  mutate(metric = recode(metric, balanced_acc = "Balanced accuracy", macro_f1 = "Macro F1"),
         model  = factor(model, levels = c("Random Forest", "LASSO", "ElasticNet", "XGBoost")))  # model order

model_colours <- c("Random Forest" = "#2C7BB6", "LASSO" = "#FDAE61",
                   "ElasticNet" = "#ABD9E9", "XGBoost" = "#D7191C")                              # colours

perf_by_donor <- perf_by_donor %>%
  mutate(donor_id    = donor,                                   # raw donor id
         donor_label = factor(donor_labels[donor_id],
                              levels = donor_labels))            # pretty donor label

fig_perf <- ggplot(perf_by_donor, aes(x = donor, y = value, colour = model, group = model)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.5)) +   # points
  geom_line(aes(group = model), position = position_dodge(width = 0.5), alpha = 0.4) +   # lines
  facet_wrap(~metric, nrow = 2, scales = "free_y") +                 # two metrics
  scale_colour_manual(values = model_colours) +
  labs(x = "Held-out donor", y = NULL, colour = NULL,
       title = "Model performance across held-out donors (LODO)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        panel.grid.minor = element_blank())

ggsave("figures/fig_performance_by_donor.pdf", fig_perf, width = 7.5, height = 6)   # save plot

mean_bal_acc <- perf_by_donor %>%                           # mean balanced accuracy per model
  dplyr::filter(metric == "Balanced accuracy") %>%
  group_by(model) %>%
  summarise(mean_balanced_acc = mean(value), .groups = "drop") %>%
  arrange(desc(mean_balanced_acc))

best_model <- as.character(mean_bal_acc$model[1])           # best model name

cm_file_pattern <- function(model_name, donor) {            # confusion matrix file naming
  switch(model_name,
         "Random Forest" = paste0("rf_confusion_matrix_", donor, ".csv"),
         "XGBoost"       = paste0("xgb_confusion_matrix_", donor, ".csv"),
         "LASSO"         = paste0("lodo_confusion_matrix_LASSO_", donor, ".csv"),
         "ElasticNet"    = paste0("lodo_confusion_matrix_ElasticNet_", donor, ".csv"),
         stop("Unknown model: ", model_name))
}

donor_ids <- unique(perf_by_donor$donor)                    # donor list
cm_files  <- vapply(donor_ids, function(d) cm_file_pattern(best_model, d), character(1))   # file names
cm_files  <- cm_files[file.exists(cm_files)]               # keep existing files

cm_list <- lapply(cm_files, function(f) {                  # load matrices
  m <- read.csv(f, row.names = 1, check.names = FALSE)
  as.matrix(m)
})

cm_sum  <- Reduce(`+`, cm_list)                            # sum confusion matrices
cm_norm <- sweep(cm_sum, 1, rowSums(cm_sum), FUN = "/")    # normalise rows

cm_df <- as.data.frame(cm_norm) %>%                        # long format
  mutate(true_class = rownames(cm_norm)) %>%
  pivot_longer(-true_class, names_to = "predicted_class", values_to = "proportion")

fig_cm <- ggplot(cm_df, aes(x = predicted_class, y = true_class, fill = proportion)) +
  geom_tile(colour = "white") +                             # heatmap tiles
  geom_text(aes(label = sprintf("%.2f", proportion)), size = 3.2,
            colour = ifelse(cm_df$proportion > 0.5, "white", "black")) +   # text colour
  scale_fill_gradient(low = "#F7FBFF", high = "#08519C", limits = c(0, 1), name = "Proportion") +
  labs(x = "Predicted cell type", y = "True cell type",
       title = paste0("Normalised confusion matrix \u2014 ", best_model),
       subtitle = "Best model by mean balanced accuracy; pooled across folds") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 13),
        plot.subtitle = element_text(size = 10),
        panel.grid = element_blank())

ggsave("figures/fig_confusion_matrix_best_model.pdf", fig_cm, width = 9, height = 6)   # save plot

shap_folds <- read.csv("genes_per_fold_xgb_shap.csv")       # SHAP table
n_folds_total_shap <- n_distinct(shap_folds$held_out_donor) # total folds

top_shap_genes <- shap_folds %>%                            # top SHAP genes
  group_by(gene, hgnc_symbol) %>%
  summarise(
    n_folds_selected  = n_distinct(held_out_donor),
    mean_abs_SHAP     = mean(mean_abs_SHAP, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abs_SHAP)) %>%
  head(30) %>%
  mutate(label = ifelse(is.na(hgnc_symbol) | hgnc_symbol == "", gene, hgnc_symbol),
         label = fct_reorder(factor(label), mean_abs_SHAP))

fig_shap <- ggplot(top_shap_genes, aes(x = mean_abs_SHAP, y = label,
                                       colour = n_folds_selected,
                                       size = n_folds_selected)) +
  geom_point() +
  scale_colour_gradient(low = "#FDBB84", high = "#B30000",
                        name = paste0("Folds present\n(of ", n_folds_total_shap, ")")) +
  scale_size_continuous(range = c(2, 5), guide = "none") +
  labs(x = "Mean |SHAP value|", y = NULL,
       title = "Top genes by mean absolute SHAP value (XGBoost)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 9))

ggsave("figures/fig_shap_top_genes.pdf", fig_shap, width = 7, height = 8)   # save plot

consensus_wide <- read.csv("consensus_gene_table.csv")      # consensus table

consensus_model_names <- c("RandomForest", "LASSO", "ElasticNet", "XGBoost_combined")   # model names
rank_cols <- paste0(consensus_model_names, "_mean_fold_rank")                           # rank cols

n_top_consensus <- 25                                        # top N

top_consensus <- consensus_wide %>%                          # top consensus genes
  arrange(desc(n_models_selected_all_folds), desc(n_models_selected_5plus_folds),
          mean_rank_across_models) %>%
  head(n_top_consensus) %>%
  mutate(label = ifelse(is.na(symbol) | symbol == "", gene, symbol))

consensus_long <- top_consensus %>%                          # long format
  dplyr::select(label, mean_rank_across_models, all_of(rank_cols)) %>%
  pivot_longer(cols = all_of(rank_cols), names_to = "model", values_to = "mean_fold_rank") %>%
  mutate(model = str_remove(model, "_mean_fold_rank$"),
         model = recode(model,
                        "RandomForest"     = "Random Forest",
                        "XGBoost_combined" = "XGBoost"),
         model = factor(model, levels = c("Random Forest", "LASSO", "ElasticNet", "XGBoost")),
         label = fct_reorder(factor(label), -mean_rank_across_models))

fig_consensus <- ggplot(consensus_long, aes(x = model, y = label, fill = mean_fold_rank)) +
  geom_tile(colour = "white") +
  scale_fill_gradient(low = "#08519C", high = "#F7FBFF",
                      name = "Mean rank\n(lower = stronger)") +
  labs(x = NULL, y = NULL,
       title = paste0("Top ", n_top_consensus, " consensus genes across all four methods"),
       subtitle = "Selected in every fold by every model; shaded by each model's mean rank") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

ggsave("figures/fig_consensus_top_genes.pdf", fig_consensus, width = 7, height = 8)   # save plot

consensus_ranked <- consensus_wide %>%                        # full ranking
  arrange(desc(n_models_selected_all_folds), desc(n_models_selected_5plus_folds),
          mean_rank_across_models) %>%
  mutate(rank_position = row_number())

panel <- read.csv("new_480_gene_panel.csv")                  # final panel

panel_ranked <- panel %>%                                    # merge panel with ranks
  left_join(consensus_ranked %>% dplyr::select(gene, rank_position), by = "gene")

n_panel_total  <- nrow(panel_ranked)                         # panel size
n_not_in_model <- sum(is.na(panel_ranked$rank_position))     # non-negotiables

fig_panel_compare <- ggplot(panel_ranked %>% dplyr::filter(!is.na(rank_position)),
                            aes(x = rank_position, fill = inclusion)) +
  geom_histogram(binwidth = 250, boundary = 0, alpha = 0.85, position = "identity") +
  scale_fill_manual(values = c("non_negotiable" = "#D7191C", "ranked_fill" = "#2C7BB6"),
                    labels = c("non_negotiable" = "Non-negotiable (forced in)",
                               "ranked_fill"     = "Ranked fill (top consensus)"),
                    name = NULL) +
  labs(x = "Rank position in full consensus gene table", y = "Number of panel genes",
       title = "480-gene panel vs pure consensus ranking",
       subtitle = paste0(n_not_in_model, " of ", n_panel_total,
                         " panel genes (all non-negotiable) never appeared in the model output")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

ggsave("figures/fig_panel_vs_consensus.pdf", fig_panel_compare, width = 7.5, height = 5.5)   # save plot

cat("best model by mean balanced accuracy:", best_model, "\n")   # print best model
