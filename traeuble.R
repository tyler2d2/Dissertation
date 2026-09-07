library(Seurat)  
library(dplyr)     
library(ggplot2)    
library(patchwork)   

set.seed(42)         # reproducibility

cat("[1/5] Loading Traeuble Supplementary Table 1 marker genes...\n")   # step 1

traeuble_markers <- read.csv(
  "/Users/tyleradams/Desktop/DISS/traeuble_supp_table1_markers.csv",
  stringsAsFactors = FALSE
)

traeuble_markers$gene <- toupper(traeuble_markers$gene)                 # uppercase

cat("Loaded", nrow(traeuble_markers), "marker rows across",
    length(unique(traeuble_markers$level1)), "Traeuble level1 categories.\n")
cat("[1/5] Done.\n\n")

cat("[2/5] Loading consensus_panel_480.csv...\n")                       # step 2

consensus <- read.csv("/Users/tyleradams/Desktop/DISS/consensus_panel_480.csv")
colnames(consensus)[1] <- "gene"
consensus$gene <- toupper(consensus$gene)

cat("Loaded", nrow(consensus), "consensus genes\n")
cat("[2/5] Done.\n\n")

cat("[3/5] Loading precomputed Seurat object...\n")                     # step 3

seurat_obj <- readRDS("/Users/tyleradams/Desktop/DISS/arterial_seurat_processed.rds")
cat("Seurat object loaded OK\n")

cat("Computing marker genes directly per cell_type...\n")

Idents(seurat_obj) <- seurat_obj$cell_type                             # set identities

markers <- FindAllMarkers(                                              # marker detection
  seurat_obj,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.1
)

cat("Loaded", nrow(markers), "marker rows across",
    length(unique(markers$cluster)), "cell types\n")

gene_category_raw <- markers %>%                                        # strongest marker per gene
  group_by(gene) %>%
  slice_max(avg_log2FC, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    gene = toupper(gene),
    category_raw = cluster
  )

cat("Built gene -> category table:", nrow(gene_category_raw), "genes\n")

rm(markers, seurat_obj)                                                 # cleanup
gc()
cat("[3/5] Done.\n\n")

cat("[4/5] Harmonising category names...\n")                            # step 4

category_map <- c(
  "fibroblast"                                     = "Fibroblasts",
  "macrophage"                                     = "Macrophages",
  "blood vessel smooth muscle cell"                = "SMCs",
  "microcirculation associated smooth muscle cell" = "SMCs",
  "endothelial cell of artery"                     = "ECs",
  "lymphocyte"                                     = "Immune",
  "unknown"                                        = "Unknown"
)

gene_category_raw$category <- category_map[gene_category_raw$category_raw]   # map names
gene_category_raw$category[is.na(gene_category_raw$category)] <-
  gene_category_raw$category_raw[is.na(gene_category_raw$category)]

cat("[4/5] Done.\n\n")

cat("[5/5] Comparing consensus panel to Traeuble marker genes...\n")    # step 5

traeuble_gene_to_categories <- split(traeuble_markers$level1, traeuble_markers$gene)   # lookup

comparison_table <- merge(                                              # merge categories
  consensus["gene"],
  gene_category_raw[, c("gene", "category")],
  by = "gene",
  all.x = TRUE
)

comparison_table$category[is.na(comparison_table$category)] <- "No category assigned"   # fill NA

comparison_table$is_traeuble_marker <- comparison_table$gene %in% names(traeuble_gene_to_categories)

comparison_table$traeuble_categories <- sapply(comparison_table$gene, function(g) {     # list categories
  ct <- traeuble_gene_to_categories[[g]]
  if (is.null(ct)) "" else paste(sort(unique(ct)), collapse = "; ")
})

write.csv(comparison_table,
          "/Users/tyleradams/Desktop/DISS/consensus_vs_traeuble_markers.csv",
          row.names = FALSE)

cat("Saved consensus_vs_traeuble_markers.csv\n")

n_total <- nrow(comparison_table)                                      # stats
n_overlap <- sum(comparison_table$is_traeuble_marker)

cat("\n---- Traeuble Marker Overlap ----\n")
cat("Consensus panel size:", n_total, "\n")
cat("Overlap with Traeuble broad marker genes:", n_overlap,
    "(", round(100 * n_overlap / n_total, 2), "% )\n\n")

by_category <- traeuble_markers %>%                                    # overlap by cell type
  group_by(level1) %>%
  summarise(
    n_marker_genes = n(),
    n_overlap = sum(gene %in% consensus$gene),
    pct_overlap = round(100 * n_overlap / n_marker_genes, 2)
  ) %>%
  arrange(desc(n_overlap))

write.csv(by_category,
          "/Users/tyleradams/Desktop/DISS/consensus_vs_traeuble_by_category.csv",
          row.names = FALSE)

cat("Saved consensus_vs_traeuble_by_category.csv\n")
cat("[5/5] Done. Script complete.\n")

cat("Generating combined Traeuble overlap bar plot...\n")               # plotting

df_gene_bar <- data.frame(                                             # gene-level bar
  category = c("Overlap", "Consensus only"),
  count = c(
    sum(comparison_table$is_traeuble_marker),
    sum(!comparison_table$is_traeuble_marker)
  )
)

p_gene <- ggplot(df_gene_bar, aes(x = category, y = count, fill = category)) +
  geom_bar(stat = "identity") +
  theme_minimal(base_size = 14) +
  labs(title = "Gene-level Overlap", x = "", y = "Number of Genes") +
  scale_fill_manual(values = c("#0072B2", "#999999")) +
  theme(legend.position = "none")

df_cell_bar <- by_category %>%                                         # cell-type bar
  mutate(level1 = factor(level1, levels = level1[order(n_overlap, decreasing = TRUE)]))

p_cell <- ggplot(df_cell_bar, aes(x = level1, y = n_overlap, fill = level1)) +
  geom_bar(stat = "identity") +
  theme_minimal(base_size = 14) +
  labs(title = "Cell-type Overlap", x = "Traeuble Level 1 Cell Type", y = "Overlapping Genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set3") +
  theme(legend.position = "none")

combined_plot <- p_gene + p_cell + plot_layout(ncol = 2)               # combine

ggsave(
  "/Users/tyleradams/Desktop/DISS/figures/traeuble_combined_overlap_barplot.jpeg",
  plot = combined_plot,
  width = 14,
  height = 6,
  dpi = 300
)

cat("Saved traeuble_combined_overlap_barplot.jpeg\n")
