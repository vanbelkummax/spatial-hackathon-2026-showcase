# Spatial Biology Hackathon 2026: PDAC Treatment Response Analysis

[![Python 3.10+](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Scanpy](https://img.shields.io/badge/Scanpy-1.10+-green.svg)](https://scanpy.readthedocs.io/)
[![Squidpy](https://img.shields.io/badge/Squidpy-1.4+-orange.svg)](https://squidpy.readthedocs.io/)
[![Polymath](https://img.shields.io/badge/Polymath-v4-purple.svg)](https://github.com/vanbelkummax/polymath-v4)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**A polymathic approach to spatial transcriptomics analysis, applying cross-domain algorithms to identify treatment response biomarkers in pancreatic ductal adenocarcinoma (PDAC).**

---

## Overview

This repository contains the analysis code, results, and methodology from the Spatial Biology Hackathon 2026 (January 20-23). We analyzed spatial transcriptomics data from PDAC patients comparing **treatment responders (R)** versus **non-responders (NR)** using a combination of traditional bioinformatics methods and cross-domain algorithms identified through the Polymath knowledge base.

### Key Innovation

Rather than relying solely on standard spatial transcriptomics workflows, we leveraged the **Polymath v4 knowledge base** (2,193 papers, 50,528 algorithms, 578K code chunks) to identify and apply algorithms from **graph theory**, **topology**, and **information theory** to spatial biology problems.

---

## Key Findings

### 1. Graph Topology Differences

| Metric | Responders | Non-Responders | Interpretation |
|--------|------------|----------------|----------------|
| Mean Betweenness | 0.0091 | 0.0047 | **+94% higher** in R |
| Mean PageRank | 0.00081 | 0.00091 | Slightly lower in R |
| Hub Cell Fraction | 5.0% | 5.0% | Similar by definition |

**Key Finding:** Responders show significantly higher betweenness centrality, indicating more interconnected tissue architecture with cells serving as critical communication bridges.

### 2. Hub Cell Types Differ

| Response | Dominant Hub Cell Types |
|----------|------------------------|
| **Responders** | Macrophages, CAF_myCAF, Endocrine |
| **Non-Responders** | NK_cells, Ductal_Epithelial |

**Interpretation:** In responders, immune cells (macrophages) and stromal cells (CAFs) occupy central positions in the spatial network. Non-responders show epithelial and NK cells as hubs, suggesting different tissue organization patterns.

### 3. Top Discriminative Genes (Mutual Information)

| Response | Sample | Top MI Gene | MI Score | Biological Role |
|----------|--------|-------------|----------|-----------------|
| R | YP12A | **CD8A** | 0.097 | T-cell marker (immune infiltration) |
| R | YP12C | CLPS | 0.234 | Colipase (pancreatic function) |
| R | YP15A | DES | 0.216 | Desmin (muscle/stroma marker) |
| R | YP15C | **TAGLN** | 0.206 | Transgelin (CAF marker) |
| NR | YP03A | **MTA1** | 0.124 | Metastasis-associated gene |
| NR | YP03C | TTR | 0.208 | Transthyretin (carrier protein) |
| NR | YP04C | **KRT7** | 0.188 | Keratin 7 (epithelial marker) |

**Key Finding:** Responders show T-cell (CD8A) and stromal (TAGLN) markers as top discriminators, while non-responders show epithelial (KRT7) and metastasis (MTA1) markers.

### 4. Novel Biomarkers (Validated via Open Targets)

These genes were upregulated in responders and validated through external databases:

| Gene | Log2FC | Disease Associations | Cancer Links |
|------|--------|---------------------|--------------|
| **MS4A2** | +0.67 | 208 diseases | Breast, gastric, colorectal cancer |
| **NLRP7** | +0.67 | 132 diseases | Colorectal carcinoma, IBD |
| **CCZ1** | +1.22 | 22 diseases | Neurodegeneration, SCC |

**Hypothesis:** MS4A2 (mast cell marker) and NLRP7 (inflammasome component) suggest enhanced immune/inflammatory activity in treatment responders.

---

## Methodology

### Data

| Dataset | Samples | Cells | Platform |
|---------|---------|-------|----------|
| PDAC Treatment Response | 8 (7 analyzed) | 8,925 | Visium |

- **Responders (R):** YP12A, YP12C, YP15A, YP15C
- **Non-Responders (NR):** YP03A, YP03C, YP04C
- *YP04A excluded (134 cells after QC)*

### Analysis Pipeline

```
Day 1: QC & Preprocessing
    |-- Per-sample QC metrics
    |-- Filtering (genes, cells, mitochondrial %)
    |-- Normalization & HVG selection

Day 2: Standard Analysis
    |-- Leiden clustering (resolution optimization)
    |-- Cell type annotation (marker-based)
    |-- Spatial analysis (neighborhood enrichment, SVGs)
    |-- Differential expression (R vs NR)

Day 3: Polymathic Analysis (Cross-Domain)
    |-- Graph centrality (from social network analysis)
    |-- Persistent homology (from algebraic topology)
    |-- Mutual information (from information theory)
```

### Cross-Domain Algorithms Applied

| Algorithm | Origin Domain | Spatial Biology Application |
|-----------|---------------|---------------------------|
| **Betweenness Centrality** | Graph Theory / Social Networks | Hub cell identification |
| **PageRank** | Graph Theory / Web Search | Cell importance scoring |
| **Persistent Homology** | Algebraic Topology | Tissue architecture quantification |
| **Betti Curves** | Topological Data Analysis | R vs NR structural signatures |
| **Mutual Information** | Information Theory | Non-linear gene-phenotype relationships |

---

## Figure Gallery

### Figure 1: Sample Overview
*Spatial distribution of cell types across all PDAC samples*
![Sample Overview](figures/fig1_sample_overview.png)

### Figure 2: Cell Type Composition
*Differential cell type enrichment between responders and non-responders*
![Cell Type Composition](figures/fig2_cell_type_composition.png)

### Figure 3: Graph Centrality Analysis
*Betweenness and PageRank comparisons showing higher connectivity in responders*
![Centrality Analysis](figures/fig3_centrality_analysis.png)

### Figure 4: Spatial Hub Cells
*Hub cell spatial distribution highlighting key network positions*
![Spatial Hub Cells](figures/fig4_spatial_hub_cells.png)

### Figure 5: Persistent Homology
*Persistence diagrams capturing tissue architecture topology*
![Persistent Homology](figures/fig5_persistent_homology.png)

### Figure 6: Betti Curves
*Topological signatures comparing R vs NR tissue structure*
![Betti Curves](figures/fig6_betti_curves.png)

### Figure 7: Mutual Information Genes
*Top discriminative genes by mutual information score*
![MI Genes](figures/fig7_mi_genes.png)

### Figure 8: Summary Dashboard
*Comprehensive overview of all key findings*
![Summary Dashboard](figures/fig8_summary_dashboard.png)

---

## Repository Structure

```
spatial-hackathon-2026-showcase/
├── README.md                    # This file
├── LICENSE                      # MIT License
├── figures/                     # Publication-quality figures
│   ├── fig1_sample_overview.png
│   ├── fig2_cell_type_composition.png
│   ├── fig3_centrality_analysis.png
│   ├── fig4_spatial_hub_cells.png
│   ├── fig5_persistent_homology.png
│   ├── fig6_betti_curves.png
│   ├── fig7_mi_genes.png
│   └── fig8_summary_dashboard.png
├── data/
│   └── results_summary.csv      # Key results in tabular format
├── notebooks/
│   └── 01_analysis_walkthrough.ipynb  # Interactive analysis notebook
├── scripts/
│   ├── 07_comparative.py        # R vs NR comparative analysis
│   ├── 08_polymathic_analysis.py  # Cross-domain algorithms
│   └── 09_generate_showcase_figures.py  # Figure generation
└── docs/
    ├── METHODOLOGY.md           # Detailed methods description
    ├── FINDINGS.md              # Comprehensive findings report
    └── POLYMATH_INTEGRATION.md  # How Polymath KB was leveraged
```

---

## Quick Start

### Prerequisites

```bash
# Create conda environment
conda create -n spatial-hackathon python=3.10
conda activate spatial-hackathon

# Install core dependencies
pip install scanpy squidpy anndata pandas numpy matplotlib seaborn

# Install TDA library (for persistent homology)
pip install giotto-tda

# Install network analysis
pip install networkx
```

### Run Analysis

```bash
# Clone repository
git clone https://github.com/vanbelkummax/spatial-hackathon-2026-showcase.git
cd spatial-hackathon-2026-showcase

# Run polymathic analysis (requires processed data)
python scripts/08_polymathic_analysis.py

# Generate figures
python scripts/09_generate_showcase_figures.py
```

---

## Polymath Knowledge Base

This analysis was grounded in the **Polymath v4** knowledge base:

| Component | Count |
|-----------|-------|
| Papers | 2,193 |
| Algorithms | 50,528 |
| Code Chunks | 578,830 |
| Repositories | 1,881 |

### Key Algorithms from Registry

```python
# Algorithms tagged with spatial biology applications
from polymath.algo import search_algorithms

# Graph theory algorithms
graph_algos = search_algorithms(domain='graph', spatial_use=True)
# Returns: betweenness_centrality, pagerank, community_detection, ...

# Topology algorithms
topo_algos = search_algorithms(domain='topology', spatial_use=True)
# Returns: persistent_homology, betti_numbers, wasserstein_distance, ...
```

### Papers Supporting Methods

1. **CellMap (Liu et al., NAR 2026)** - Spatial transcriptomic cellular landscape mapping
2. **Singh et al. (2023)** - TDA applications in medical imaging
3. **Bulletin of Mathematical Biology (2026)** - Persistent homology parameter classification

---

## Authors

**Max Van Belkum**
MD-PhD Student, Vanderbilt University
Research Focus: Spatial transcriptomics, computational pathology

---

## Acknowledgments

- **Spatial Biology Hackathon 2026** organizers
- **Polymath v4** knowledge base development
- **Vanderbilt University** Department of Biomedical Informatics
- **scverse** community (Scanpy, Squidpy, AnnData)

---

## Citation

If you use this analysis or methodology, please cite:

```bibtex
@misc{vanbelkum2026pdac,
  author = {Van Belkum, Max},
  title = {Spatial Biology Hackathon 2026: PDAC Treatment Response Analysis},
  year = {2026},
  publisher = {GitHub},
  url = {https://github.com/vanbelkummax/spatial-hackathon-2026-showcase}
}
```

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

*Generated using Polymath v4 RRF search (lexical + semantic + Neo4j graph)*
*Algorithm Registry: 50,528 algorithms across 15+ domains with spatial biology applications*
