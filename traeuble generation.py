import scanpy as sc
import pandas as pd

print("LOAD")
adata = sc.read_h5ad("/Users/tyleradams/Desktop/DISS/ad9f52c0-054b-4dc5-944c-46c376a544cf.h5ad")
adata.raw = adata
ct_col = "cell_type_level1"
print("OK")

print("DE")
sc.tl.rank_genes_groups(
    adata,
    groupby=ct_col,
    method="wilcoxon",
    use_raw=True,
    n_genes=adata.raw.var.shape[0]
)
rg = adata.uns["rank_genes_groups"]
print("OK")

print("TIDY")
groups = rg["names"].dtype.names

de = pd.DataFrame([
    {
        "cell_type_level1": g,
        "gene": gene,
        "p_adj": padj,
        "log2FC": lfc,
        "pct_1": pct1,
        "pct_2": pct2
    }
    for g in groups
    for gene, padj, lfc, pct1, pct2 in zip(
        rg["names"][g],
        rg["pvals_adj"][g],
        rg["logfoldchanges"][g],
        rg["pct.1"][g],
        rg["pct.2"][g]
    )
])
print("OK")

print("FILTER")
filtered = de[
    (de["p_adj"] < 0.05) &
    (de["log2FC"] >= 0.25) &
    ((de["pct_1"] - de["pct_2"]) >= 0.05)
]
print("OK")

print("COUNT")
marker_count_total = len(filtered)
marker_count_per_type = filtered.groupby("cell_type_level1").size()
print("OK")

print("SAVE")
filtered.to_csv("plaque_atlas_level1_markers_filtered.csv", index=False)
marker_count_per_type.to_csv("plaque_atlas_level1_marker_counts.csv")
print("DONE")
