# Methodology

## Spatial Biology Hackathon 2026 - PDAC Treatment Response Analysis

This document provides a detailed description of the methods used in our analysis of treatment response in pancreatic ductal adenocarcinoma (PDAC) spatial transcriptomics data.

---

## Table of Contents

1. [Day 1: Quality Control and Preprocessing](#day-1-quality-control-and-preprocessing)
2. [Day 2: Clustering, Annotation, and Spatial Analysis](#day-2-clustering-annotation-and-spatial-analysis)
3. [Day 3: Polymathic Cross-Domain Analysis](#day-3-polymathic-cross-domain-analysis)
4. [Statistical Methods](#statistical-methods)
5. [Software and Tools](#software-and-tools)

---

## Day 1: Quality Control and Preprocessing

### 1.1 Data Loading

Data was loaded from Visium spatial transcriptomics experiments using Scanpy's `read_visium()` function. Each sample included:
- Gene expression count matrix
- Spatial coordinates
- Tissue images (H&E)

```python
import scanpy as sc

adata = sc.read_visium(
    path=sample_dir,
    count_file="filtered_feature_bc_matrix.h5"
)
```

### 1.2 Quality Control Metrics

For each sample, we computed:

| Metric | Description | Typical Threshold |
|--------|-------------|-------------------|
| `n_genes_by_counts` | Number of genes detected per spot | > 200 |
| `total_counts` | Total UMI counts per spot | > 500 |
| `pct_counts_mt` | Percentage mitochondrial genes | < 20% |

```python
# Calculate QC metrics
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=['mt'],
    percent_top=None,
    log1p=False,
    inplace=True
)
```

### 1.3 Filtering

Spots and genes were filtered based on quality thresholds:

```python
# Filter cells
sc.pp.filter_cells(adata, min_genes=200)
adata = adata[adata.obs['pct_counts_mt'] < 20, :]

# Filter genes
sc.pp.filter_genes(adata, min_cells=10)
```

### 1.4 Normalization

We applied total-count normalization followed by log transformation:

```python
# Normalize to 10,000 counts per spot
sc.pp.normalize_total(adata, target_sum=1e4)

# Log transform
sc.pp.log1p(adata)

# Store raw counts
adata.raw = adata
```

### 1.5 Highly Variable Gene Selection

Highly variable genes (HVGs) were identified using the `seurat_v3` method:

```python
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=3000,
    flavor='seurat_v3',
    batch_key=None
)
```

---

## Day 2: Clustering, Annotation, and Spatial Analysis

### 2.1 Dimensionality Reduction

**PCA** was performed on HVGs:

```python
sc.tl.pca(adata, n_comps=50, use_highly_variable=True)
```

**UMAP** was computed for visualization:

```python
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)
```

### 2.2 Clustering

We used Leiden clustering with multi-resolution analysis:

```python
import squidpy as sq

# Build spatial neighbors graph
sq.gr.spatial_neighbors(adata, n_neighs=6)

# Leiden clustering
sc.tl.leiden(adata, resolution=0.5, key_added='leiden_0.5')
sc.tl.leiden(adata, resolution=1.0, key_added='leiden_1.0')
```

Resolution selection was guided by:
- Silhouette score
- Calinski-Harabasz index
- Biological interpretability

### 2.3 Cell Type Annotation

Marker-based annotation was performed using canonical PDAC markers:

| Cell Type | Markers |
|-----------|---------|
| Ductal Epithelial | KRT19, KRT7, EPCAM |
| Acinar | PRSS1, AMY2A, CPA1 |
| Endocrine | INS, GCG, SST |
| CAF_myCAF | ACTA2, TAGLN, MYH11 |
| CAF_iCAF | IL6, CXCL12, PDGFRA |
| Macrophage | CD68, CD163, CSF1R |
| T_cells | CD3D, CD8A, CD4 |
| NK_cells | NKG7, GNLY, KLRD1 |
| B_cells | CD79A, MS4A1 |

```python
def annotate_cell_types(adata, markers_dict):
    """Score cells for each cell type based on marker expression."""
    for cell_type, markers in markers_dict.items():
        present_markers = [m for m in markers if m in adata.var_names]
        if present_markers:
            sc.tl.score_genes(adata, present_markers, score_name=f'{cell_type}_score')
```

### 2.4 Spatial Analysis

#### Neighborhood Enrichment

Computed pairwise cell type spatial co-localization:

```python
sq.gr.nhood_enrichment(
    adata,
    cluster_key='cell_type',
    n_perms=1000
)
```

#### Spatially Variable Genes (SVGs)

Identified using Moran's I autocorrelation:

```python
sq.gr.spatial_autocorr(
    adata,
    mode='moran',
    genes=adata.var_names[:5000]
)
```

### 2.5 Differential Expression

Responder vs non-responder comparison using Wilcoxon rank-sum test:

```python
# Merge samples by response group
merged = anndata.concat([r_samples, nr_samples])
merged.obs['response'] = ...

# Run DE
sc.tl.rank_genes_groups(
    merged,
    groupby='response',
    method='wilcoxon',
    reference='NR'
)
```

---

## Day 3: Polymathic Cross-Domain Analysis

This is the key innovation of our approach - applying algorithms from domains outside biology to spatial transcriptomics.

### 3.1 Graph Centrality Analysis

**Rationale:** Social network analysis algorithms can identify "influential" nodes in a network. We hypothesized that cells with high centrality play important roles in tissue organization.

**Methods:**

1. **Betweenness Centrality** - Measures how often a cell lies on shortest paths between other cells:

```python
import networkx as nx

# Build spatial graph
sq.gr.spatial_neighbors(adata, n_neighs=6)
G = nx.from_scipy_sparse_array(adata.obsp['spatial_connectivities'])

# Compute betweenness
betweenness = nx.betweenness_centrality(G)
adata.obs['betweenness'] = [betweenness.get(i, 0) for i in range(adata.n_obs)]
```

2. **PageRank** - Originally for web page ranking, identifies cells connected to other "important" cells:

```python
pagerank = nx.pagerank(G)
adata.obs['pagerank'] = [pagerank.get(i, 0) for i in range(adata.n_obs)]
```

3. **Hub Cell Identification** - Top 5% by PageRank:

```python
threshold = np.percentile(adata.obs['pagerank'], 95)
adata.obs['is_hub_cell'] = adata.obs['pagerank'] > threshold
```

### 3.2 Persistent Homology (Topological Data Analysis)

**Rationale:** Persistent homology captures topological features (connected components, holes, voids) that persist across multiple scales. This can quantify tissue architecture in a way that's robust to noise and deformation.

**Methods:**

We used the Vietoris-Rips complex to compute persistence diagrams:

```python
from gtda.homology import VietorisRipsPersistence
from gtda.diagrams import BettiCurve

# Normalize spatial coordinates
coords = adata.obsm['spatial']
coords = (coords - coords.min(axis=0)) / (coords.max(axis=0) - coords.min(axis=0) + 1e-6)

# Compute persistence
VR = VietorisRipsPersistence(
    metric='euclidean',
    homology_dimensions=[0, 1],  # H0=components, H1=holes
    n_jobs=-1
)
diagrams = VR.fit_transform([coords])

# Extract Betti curves as features
betti = BettiCurve(n_bins=100)
features = betti.fit_transform(diagrams)
```

**Interpretation:**
- **H0 (Betti-0):** Number of connected components at each scale
- **H1 (Betti-1):** Number of 1-dimensional holes (loops) at each scale

### 3.3 Mutual Information for Gene Selection

**Rationale:** Mutual information captures non-linear dependencies that Pearson correlation misses. This is important for gene expression data where relationships are often non-linear.

**Methods:**

```python
from sklearn.feature_selection import mutual_info_classif

# Compute MI between genes and cell type
X = adata.X.toarray()
y = adata.obs['cell_type'].astype('category').cat.codes.values

mi_scores = mutual_info_classif(X, y, random_state=42)

# Rank genes by MI
mi_df = pd.DataFrame({
    'gene': adata.var_names,
    'mutual_information': mi_scores
}).sort_values('mutual_information', ascending=False)
```

---

## Statistical Methods

### Multiple Testing Correction

- Differential expression: Benjamini-Hochberg FDR correction
- Neighborhood enrichment: Permutation-based p-values (1000 permutations)

### Group Comparisons

- **Two-group continuous:** Mann-Whitney U test (non-parametric)
- **Cell type proportions:** Fisher's exact test or chi-squared

### Effect Size Metrics

- Log2 fold change for expression differences
- Percentage difference for centrality metrics

---

## Software and Tools

### Core Analysis Stack

| Tool | Version | Purpose |
|------|---------|---------|
| Python | 3.10 | Programming language |
| Scanpy | 1.10+ | Single-cell analysis |
| Squidpy | 1.4+ | Spatial transcriptomics |
| AnnData | 0.10+ | Data structure |
| NetworkX | 3.0+ | Graph algorithms |
| giotto-tda | 0.6+ | Topological data analysis |
| scikit-learn | 1.4+ | Mutual information |

### Visualization

| Tool | Purpose |
|------|---------|
| Matplotlib | Base plotting |
| Seaborn | Statistical visualization |
| Squidpy | Spatial plots |

### Knowledge Base

| Resource | Purpose |
|----------|---------|
| Polymath v4 | Algorithm discovery (50K+ algorithms) |
| Open Targets | Gene-disease validation |
| PubMed | Literature validation |

---

## Reproducibility

All random seeds were set for reproducibility:

```python
import numpy as np
np.random.seed(42)

# In TDA
VR = VietorisRipsPersistence(..., n_jobs=-1)  # Deterministic

# In MI
mi_scores = mutual_info_classif(X, y, random_state=42)
```

### Computational Resources

- **CPU:** 24 cores
- **RAM:** 196 GB
- **GPU:** NVIDIA RTX 5090 24GB (not required for this analysis)
- **OS:** Ubuntu 22.04 (WSL2)

---

## References

1. Wolf FA, et al. (2018). SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*.
2. Palla G, et al. (2022). Squidpy: a scalable framework for spatial omics analysis. *Nature Methods*.
3. Otter N, et al. (2017). A roadmap for the computation of persistent homology. *EPJ Data Science*.
4. Cover TM, Thomas JA. (2006). Elements of Information Theory. *Wiley*.

---

*Last updated: 2026-01-19*
