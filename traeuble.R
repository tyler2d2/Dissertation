library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

set.seed(42)   # keep runs reproducible

cat("[1/5] Loading Traeuble atlas level1 markers (filtered)...\n")

traeuble_markers <- read.csv(
  "/Users/tyleradams/Desktop/DISS/plaque_atlas_level1_markers_filtered.csv",
  stringsAsFactors = FALSE
)

colnames(traeuble_markers) <- tolower(colnames(traeuble_markers))

# normalise column names + uppercase genes
traeuble_markers <- traeuble_markers %>%
  transmute(
    gene = toupper(gene),
    level1 = cell_type_level1
  )

cat("Loaded", nrow(traeuble_markers), "Traeuble markers across",
    length(unique(traeuble_markers$level1)), "level1 categories.\n")
cat("[1/5] Done.\n\n")


cat("[2/5] Loading consensus panel...\n")

consensus <- read.csv("/Users/tyleradams/Desktop/DISS/consensus_panel_480.csv")
colnames(consensus)[1] <- "gene"
consensus$gene <- toupper(consensus$gene)

cat("Loaded", nrow(consensus), "consensus genes.\n")
cat("[2/5] Done.\n\n")


cat("[3/5] Loading Seurat object + computing markers...\n")

seurat_obj <- readRDS("/Users/tyleradams/Desktop/DISS/arterial_seurat_processed.rds")
cat("Seurat object loaded.\n")

Idents(seurat_obj) <- seurat_obj$cell_type   # use provided cell-type labels

markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.1
)

cat("Detected", nrow(markers), "marker rows across",
    length(unique(markers$cluster)), "cell types.\n")

# strongest marker per gene
gene_category_raw <- markers %>%
  group_by(gene) %>%
  slice_max(avg_log2FC, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    gene = toupper(gene),
    category_raw = cluster
  )

cat("Built gene→category table with", nrow(gene_category_raw), "genes.\n")

rm(markers, seurat_obj)
gc()
cat("[3/5] Done.\n\n")


cat("[4/5] Harmonising category names...\n")

category_map <- c(
  "fibroblast"                                     = "Fibroblasts",
  "macrophage"                                     = "Macrophages",
  "blood vessel smooth muscle cell"                = "SMCs",
  "microcirculation associated smooth muscle cell" = "SMCs",
  "endothelial cell of artery"                     = "ECs",
  "lymphocyte"                                     = "Immune",
  "unknown"                                        = "Unknown"
)

# map raw categories → harmonised categories
gene_category_raw$category <- category_map[gene_category_raw$category_raw]
gene_category_raw$category[is.na(gene_category_raw$category)] <-
  gene_category_raw$category_raw[is.na(gene_category_raw$category)]

cat("[4/5] Done.\n\n")


cat("[5/5] Comparing consensus panel to Traeuble markers...\n")

# lookup table: gene → list of level1 categories
traeuble_gene_to_categories <- split(traeuble_markers$level1, traeuble_markers$gene)

# merge consensus genes with Seurat-derived categories
comparison_table <- merge(
  consensus["gene"],
  gene_category_raw[, c("gene", "category")],
  by = "gene",
  all.x = TRUE
)

comparison_table$category[is.na(comparison_table$category)] <- "No category assigned"

comparison_table$is_traeuble_marker <- comparison_table$gene %in% names(traeuble_gene_to_categories)

# attach Traeuble level1 categories 
comparison_table$traeuble_categories <- sapply(comparison_table$gene, function(g) {
  ct <- traeuble_gene_to_categories[[g]]
  if (is.null(ct)) "" else paste(sort(unique(ct)), collapse = "; ")
})

write.csv(comparison_table,
          "/Users/tyleradams/Desktop/DISS/consensus_vs_traeuble_markers.csv",
          row.names = FALSE)

cat("Saved consensus_vs_traeuble_markers.csv\n")

n_total <- nrow(comparison_table)
n_overlap <- sum(comparison_table$is_traeuble_marker)

cat("\n---- Traeuble Marker Overlap ----\n")
cat("Consensus panel size:", n_total, "\n")
cat("Overlap:", n_overlap, " (", round(100 * n_overlap / n_total, 2), "% )\n\n")

# overlap by cell-type category
by_category <- traeuble_markers %>%
  group_by(level1) %>%
  summarise(
    n_marker_genes = n(),
    n_overlap = sum(gene %in% consensus$gene),
    pct_overlap = round(100 * n_overlap / n_marker_genes, 2)
  ) %>%
  arrange(desc(n_overlap))

write.csv(by_category,
          "/Users/tyleradams/Desktop/DISS/consensus_vs
