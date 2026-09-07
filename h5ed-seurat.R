library(BPCells)
library(dplyr)
library(zellkonverter)
library(Seurat)
library(SingleCellExperiment)
library(presto)
library(ggplot2)
set.seed(42)

# Step 1: Load h5ad
sce <- readH5AD("DATA/441fd653-78f4-4654-b3be-1b0bc79b5568.h5ad")

# Step 2: Save metadata FIRST before anything else
meta <- as.data.frame(colData(sce))

# Step 3: Extract counts and write to disk
counts_mat <- assay(sce, "X")
write_matrix_dir(mat = counts_mat, dir = "arterial_counts_bp")

# Step 4: Free RAM
rm(sce, counts_mat)
gc()

# Step 5: Load counts back as on-disk matrix
counts_bp <- open_matrix_dir("arterial_counts_bp")

# Step 6: Build Seurat object
seurat_obj <- CreateSeuratObject(counts = counts_bp, meta.data = meta)
rm(meta)
gc()

# Step 7: Check object
dim(seurat_obj)
head(seurat_obj@meta.data)
saveRDS(object = seurat_obj, file = "DATA/seurat_object.rds")

# QC plots
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
summary(seurat_obj$nFeature_RNA)
summary(seurat_obj$nCount_RNA)

# Normalisation + HVGs
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 5000)

top10 <- head(VariableFeatures(seurat_obj), 10)
plot1 <- VariableFeaturePlot(seurat_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scaling
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)

# PCA
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
ElbowPlot(seurat_obj, ndims = 50)

VizDimLoadings(seurat_obj, dims = 1:2, reduction = "pca")
DimPlot(seurat_obj, reduction = "pca")
DimHeatmap(seurat_obj, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(seurat_obj, dims = 16:30, cells = 500, balanced = TRUE)

# Clustering + UMAP
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "cell_type", label = TRUE)
p2 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE)
p1
p2

# Marker detection
if (!file.exists("arterial_markers.rds")) {
  seurat_obj.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)
  saveRDS(seurat_obj.markers, "arterial_markers.rds")
} else {
  seurat_obj.markers <- readRDS("arterial_markers.rds")
}

top_genes <- seurat_obj.markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 1) %>%
  pull(gene)

FeaturePlot(seurat_obj, features = top_genes)
p1 + p2

saveRDS(seurat_obj, "arterial_seurat_processed.rds")

# Cluster → cell-type mapping
table(seurat_obj$cell_type, seurat_obj$seurat_clusters)

smc_clusters      <- c(9, 10, 11, 13, 21)
stromal_clusters  <- c(3, 6, 7, 12, 16, 17, 20)
immune_clusters   <- c(0, 2, 22, 23, 24, 25, 26, 27, 28, 29, 30)

top_genes_smc <- seurat_obj.markers %>%
  filter(cluster %in% smc_clusters) %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 1) %>%
  pull(gene)

top_genes_stromal <- seurat_obj.markers %>%
  filter(cluster %in% stromal_clusters) %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 1) %>%
  pull(gene)

top_genes_immune <- seurat_obj.markers %>%
  filter(cluster %in% immune_clusters) %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 1) %>%
  pull(gene)

fp_smc     <- FeaturePlot(seurat_obj, features = top_genes_smc)
fp_stromal <- FeaturePlot(seurat_obj, features = top_genes_stromal)
fp_immune  <- FeaturePlot(seurat_obj, features = top_genes_immune)

# Save plots
ggsave("figures/01_qc_violin.jpeg", VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2), width = 10, height = 6, dpi = 300)
ggsave("figures/02_variable_features.jpeg", plot1 + plot2, width = 12, height = 6, dpi = 300)
ggsave("figures/03_elbow_plot.jpeg", ElbowPlot(seurat_obj, ndims = 50), width = 8, height = 6, dpi = 300)
ggsave("figures/04_pca_plot.jpeg", DimPlot(seurat_obj, reduction = "pca"), width = 8, height = 6, dpi = 300)
ggsave("figures/05_umap_cell_type.jpeg", p1, width = 10, height = 8, dpi = 300)
ggsave("figures/06_umap_clusters.jpeg", p2, width = 10, height = 8, dpi = 300)
ggsave("figures/07_umap_combined.jpeg", p1 + p2, width = 18, height = 8, dpi = 300)
ggsave("figures/08a_feature_plot_smc.jpeg", fp_smc, width = 16, height = 12, dpi = 300)
ggsave("figures/08b_feature_plot_stromal.jpeg", fp_stromal, width = 16, height = 12, dpi = 300)
ggsave("figures/08c_feature_plot_immune.jpeg", fp_immune, width = 16, height = 12, dpi = 300)

jpeg("figures/09_dim_heatmap_1_15.jpeg", width = 1200, height = 1600, res = 150)
DimHeatmap(seurat_obj, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

jpeg("figures/10_dim_heatmap_16_30.jpeg", width = 1200, height = 1600, res = 150)
DimHeatmap(seurat_obj, dims = 16:30, cells = 500, balanced = TRUE)
dev.off()

top10_markers <- seurat_obj.markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup()

jpeg("figures/11_marker_heatmap.jpeg", width = 2000, height = 1600, res = 150)
DoHeatmap(seurat_obj, features = top10_markers$gene) + NoLegend()
dev.off()

# Reload processed object
seurat_obj <- readRDS("arterial_seurat_processed.rds")

# Confirm integrity
dim(seurat_obj)
Reductions(seurat_obj)
head(seurat_obj@meta.data)
