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

### 1. Graph Topology Differences (with Statistics)

| Metric | Responders | Non-Responders | Change | p-value |
|--------|------------|----------------|--------|---------|
| Mean Betweenness | 0.0091 | 0.0047 | **+96%** | 0.029* |
| Mean PageRank | 0.00081 | 0.00091 | -12% | 0.343 |
| Betti-0 AUC | 765.8 | 629.2 | +22% | 0.114 |
| Betti-1 AUC | 787.4 | 628.7 | **+25%** | 0.057 |

*Mann-Whitney U test, *p<0.05

**Key Finding:** Responders show significantly higher betweenness centrality (+96%, p=0.029), indicating more interconnected tissue architecture with cells serving as critical communication bridges.

### 2. Hub Cell Types Differ by Response

| Response | Dominant Hub Cell Types |
|----------|------------------------|
| **Responders** | Macrophages, CAF_myCAF, Endocrine |
| **Non-Responders** | NK_cells, Ductal_Epithelial |

**Interpretation:** In responders, immune cells (macrophages) and stromal cells (CAFs) occupy central positions in the spatial network. Non-responders show epithelial and NK cells as hubs, suggesting different tissue organization patterns.

### 3. Treatment Response Biomarkers (MI-based Discovery)

**NEW: Cross-sample mutual information analysis against R/NR status (8,925 cells pooled)**

| Rank | Gene | MI Score | Biological Function |
|------|------|----------|---------------------|
| 1 | **ACSS3** | 0.668 | Acetyl-CoA synthetase (metabolism) |
| 2 | **SOCS2** | 0.665 | JAK-STAT pathway suppressor |
| 3 | **STAT4** | 0.664 | Th1 differentiation, IFN-γ signaling |
| 4 | **CMKLR1** | 0.662 | Chemokine receptor |
| 5 | **SOBP** | 0.662 | Sine oculis binding protein |

**Key Finding:** STAT4 and SOCS2 are known immune regulators with established links to immunotherapy response, validating our cross-domain approach.

### 4. Topological Data Analysis Results

**Persistent Homology (from Algebraic Topology):**
- H0 (connected components): R samples show +22% higher AUC
- H1 (holes/voids): R samples show +25% higher AUC (p=0.057)
- **Interpretation:** Responder tissues have more complex spatial architecture with more "holes" potentially facilitating drug penetration and immune cell access

### 5. Cell Type Composition (Mann-Whitney U)

| Cell Type | R Mean% | NR Mean% | Effect | p-value |
|-----------|---------|----------|--------|---------|
| **Acinar** | 6.4% | 2.9% | **+3.5%** | **0.057†** |
| Ductal_Epithelial | 10.7% | 15.6% | -4.9% | 0.229 |
| NK_cells | 7.4% | 4.9% | +2.5% | 0.229 |

*Mann-Whitney U test; † p<0.1 (trending)

**Key Finding:** Acinar cells are the ONLY cell type showing a trending difference between responders and non-responders (2.2x higher in R, p=0.057). This suggests preserved exocrine function may associate with treatment response.

---

## Datasets Analyzed

### 1. PDAC Treatment Response (Primary)

| Dataset | Samples | Cells | Platform |
|---------|---------|-------|----------|
| PDAC Treatment Response | 7 analyzed | 8,925 | Visium |

- **Responders (R):** YP12A (865), YP12C (967), YP15A (2,159), YP15C (1,407)
- **Non-Responders (NR):** YP03A (908), YP03C (1,063), YP04C (1,556)
- *YP04A excluded (134 cells after QC)*

### 2. G4X Multimodal (Secondary - Gastric Cancer)

| Dataset | Samples | Cells | Platform |
|---------|---------|-------|----------|
| Gastric Cancer | 2 | 46,208 | Singular G4X |

- **A01:** 19,947 cells (RNA + 17 protein markers)
- **B01:** 26,261 cells (RNA + 17 protein markers)
- Rich immune infiltrate: CD8+ T (11-15%), CD4+ T (7%), Tregs (5%), Macrophages (3%)

---

## Methodology

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
    |-- MaxMin landmark sampling (topology-preserving)
```

### Cross-Domain Algorithms Applied

| Algorithm | Origin Domain | Spatial Biology Application | Implementation |
|-----------|---------------|---------------------------|----------------|
| **Betweenness Centrality** | Graph Theory / Social Networks | Hub cell identification | NetworkX |
| **PageRank** | Graph Theory / Web Search | Cell importance scoring | NetworkX |
| **Persistent Homology** | Algebraic Topology | Tissue architecture quantification | giotto-tda |
| **Betti Curves** | Topological Data Analysis | R vs NR structural signatures | giotto-tda |
| **Mutual Information** | Information Theory | Non-linear gene-phenotype relationships | sklearn |
| **MaxMin Sampling** | Computational Geometry | Topology-preserving subsampling | Custom |

---

## Figure Gallery

### Figure 1: Sample Overview
*Spatial distribution of cell types across all 7 PDAC samples with adaptive point sizing*
![Sample Overview](figures/fig1_sample_overview.png)

### Figure 2: Cell Type Composition
*Differential cell type enrichment with Mann-Whitney U test*
![Cell Type Composition](figures/fig2_cell_type_composition.png)

### Figure 3: Graph Centrality Analysis
*Betweenness and PageRank comparisons with Mann-Whitney U statistics*
![Centrality Analysis](figures/fig3_centrality_analysis.png)

### Figure 4: Spatial Hub Cells
*Hub cell spatial distribution across all samples with interpretation panel*
![Spatial Hub Cells](figures/fig4_spatial_hub_cells.png)

### Figure 5: Persistent Homology
*Persistence diagrams with R vs NR comparison and interpretation*
![Persistent Homology](figures/fig5_persistent_homology.png)

### Figure 6: Betti Curves
*Topological signatures with AUC statistics and explanation panel*
![Betti Curves](figures/fig6_betti_curves.png)

### Figure 7: Mutual Information Analysis
*Cell-type markers vs TRUE response biomarkers (cross-sample MI)*
![MI Genes](figures/fig7_mi_genes.png)

### Figure 8: Summary Dashboard
*Comprehensive overview of all key findings*
![Summary Dashboard](figures/fig8_summary_dashboard.png)

### Figure 9: Treatment Response Biomarkers (NEW)
*Dedicated biomarker figure with MI-based candidates and validation status*
![Biomarker Focus](figures/fig9_biomarker_focus.png)

### Figure 10: Sample Gallery (NEW)
*Traditional overview showing all samples with QC metrics*
![Sample Gallery](figures/fig10_sample_gallery.png)

### Figure 11: Acinar Cell Comparison (NEW)
*Dedicated analysis of Acinar cell enrichment in responders (MWU p=0.057)*
![Acinar Comparison](figures/fig11_acinar_comparison.png)

### G4X Multimodal Analysis
*Gastric cancer samples with protein-based cell type annotation*
![G4X A01](figures/g4x/A01_overview.png)
![G4X B01](figures/g4x/B01_overview.png)

---

## Repository Structure

```
spatial-hackathon-2026-showcase/
├── README.md                    # This file
├── LICENSE                      # MIT License
├── figures/                     # Publication-quality figures (11 PDAC + 2 G4X)
│   ├── fig1_sample_overview.png/pdf
│   ├── fig2_cell_type_composition.png/pdf
│   ├── fig3_centrality_analysis.png/pdf
│   ├── fig4_spatial_hub_cells.png/pdf
│   ├── fig5_persistent_homology.png/pdf
│   ├── fig6_betti_curves.png/pdf
│   ├── fig7_mi_genes.png/pdf
│   ├── fig8_summary_dashboard.png/pdf
│   ├── fig9_biomarker_focus.png/pdf
│   ├── fig10_sample_gallery.png/pdf
│   └── g4x/                     # G4X multimodal figures
├── data/
│   └── tables/
│       ├── mi_vs_response_biomarkers.csv   # Treatment response biomarkers
│       ├── polymathic_analysis_results.csv # Per-sample results
│       └── sample_overview_stats.csv       # QC metrics
├── scripts/
│   ├── 08_polymathic_analysis.py           # Cross-domain algorithms
│   ├── 09_generate_showcase_figures_fixed.py  # Figure generation
│   └── 10_g4x_multimodal_analysis.py       # G4X multimodal pipeline
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
pip install scanpy squidpy anndata pandas numpy matplotlib seaborn scipy

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
python scripts/09_generate_showcase_figures_fixed.py

# Run G4X multimodal analysis
python scripts/10_g4x_multimodal_analysis.py
```

---

## Polymath Knowledge Base Integration

This analysis was grounded in the **Polymath v4** knowledge base:

| Component | Count |
|-----------|-------|
| Papers | 2,193 |
| Algorithms | 50,528 |
| Code Chunks | 578,830 |
| Repositories | 1,882 (incl. G4X-viewer) |

### Key Algorithms from Registry

```bash
# Query Polymath for spatial biology algorithms
cd /home/user/polymath-v4
python scripts/algo.py --spatial --top 20

# Returns: optimal_transport, persistent_homology, community_detection,
#          pagerank, spectral_clustering, leiden_algorithm, ...

# Search for specific algorithms
python scripts/algo.py "wasserstein"
python scripts/q.py "optimal transport spatial transcriptomics" -n 5
```

### Papers Supporting Methods

1. **DeST-OT (2026)** - Spatiotemporal transcriptomics alignment via optimal transport
2. **Persistent Homology Classification (Bull. Math. Biol. 2026)** - TDA for parameter spaces
3. **MEcell (NAR 2026)** - Spatial transcriptomics simulation benchmarks
4. **CellMap (Liu et al.)** - Spatial cellular landscape mapping

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
- **Singular Genomics** (G4X platform and viewer)

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
*Last updated: 2026-01-19*
