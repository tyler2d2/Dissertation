import scanpy as sc
import numpy as np
import pandas as pd

print("Load")
adata = sc.read_h5ad("/Users/tyleradams/Desktop/DISS/ad9f52c0-054b-4dc5-944c-46c376a544cf.h5ad")
print("Loaded")

ct_col = "cell_type_level1"
adata.raw = adata

print("DE")
sc.tl.rank_genes_groups(
    adata,
    groupby=ct_col,
    method="t-test_overestim_var",
    use_raw=True,
    n_genes=adata.raw.var.shape[0]
)
print("DEdone")

print("Tidy")
rg = adata.uns["rank_genes_groups"]
groups = rg["names"].dtype.names

rows = []
for g in groups:
    for gene, padj, lfc in zip(rg["names"][g], rg["pvals_adj"][g], rg["logfoldchanges"][g]):
        rows.append({
            "cell_type_level1": g,
            "gene": gene,
            "p_adj": padj,
            "log2FC": lfc
        })

de = pd.DataFrame(rows)
print("Tidydone")

print("Freq")
X = adata.raw.X
genes = adata.raw.var_names
ct = adata.obs[ct_col]

print("Global")
global_counts = np.array((X > 0).sum(axis=0)).ravel()
freq_global = global_counts / X.shape[0]

print("PerType")
freq_target = {}
for celltype in ct.unique():
    idx = np.where(ct == celltype)[0]
    subX = X[idx, :]
    counts = np.array((subX > 0).sum(axis=0)).ravel()
    freq_target[celltype] = counts / len(idx)

print("Freqdone")

print("Attach")
records = []
for _, row in de.iterrows():
    g = row["cell_type_level1"]
    gene = row["gene"]

    if gene not in genes:
        continue

    gi = genes.get_loc(gene)
    ft = freq_target[g][gi]
    fg = freq_global[gi]
    diff = ft - fg

    records.append({**row, "freq_target": ft, "freq_global": fg, "freq_diff": diff})

de_freq = pd.DataFrame(records)

print("Filter")
filtered = de_freq[
    (de_freq["p_adj"] < 0.05) &
    (de_freq["log2FC"] >= 0.25) &
    (de_freq["freq_target"] >= 0.05) &
    (de_freq["freq_diff"] >= 0.05)
]
print("Filterdone")

print("Save")
filtered.to_csv("plaque_atlas_level1_markers_filtered.csv", index=False)
print("Done")
