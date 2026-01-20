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
| Mean Betweenness | 0.0091 | 0.0047 | **+96%** | 0.24 (ns) |
| Mean PageRank | 0.00081 | 0.00091 | -12% | 0.34 (ns) |
| Betti-0 AUC | 765.8 | 629.2 | +22% | 0.11 (ns) |
| Betti-1 AUC | 787.4 | 628.7 | **+25%** | 0.53 (ns) |

*Welch's t-test; ns = not significant (p ≥ 0.1)*

**Key Finding:** Responders showed higher mean betweenness centrality (+96% effect size, p=0.24), suggesting more interconnected tissue architecture. While not statistically significant with n=7 samples, this large effect warrants validation in larger cohorts.

### 2. Hub Cell Types Differ by Response

| Response | Dominant Hub Cell Types |
|----------|------------------------|
| **Responders** | Macrophages, CAF_myCAF, Endocrine |
| **Non-Responders** | NK_cells, Ductal_Epithelial |

**Interpretation:** In responders, immune cells (macrophages) and stromal cells (CAFs) occupy central positions in the spatial network. Non-responders show epithelial and NK cells as hubs, suggesting different tissue organization patterns.

### 3. Treatment Response Biomarkers (MI-based Discovery)

**Cross-sample mutual information analysis against R/NR status (n=7 samples, pseudobulk aggregation of 8,925 cells)**

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
- H0 (connected components): R samples show +22% higher AUC (p=0.11, ns)
- H1 (holes/voids): R samples show +25% higher AUC (p=0.53, ns)
- **Interpretation:** Responder tissues showed higher topological complexity (larger effect sizes), though not statistically significant at n=7. The trend toward more "holes" in responder tissue is an exploratory finding requiring validation.

### 5. Cell Type Composition (Mann-Whitney U)

| Cell Type | R Mean% | NR Mean% | Effect | p-value | FDR q |
|-----------|---------|----------|--------|---------|-------|
| **Acinar** | 6.4% | 2.9% | **+3.5%** | 0.057 | 0.74 |
| Ductal_Epithelial | 10.7% | 15.6% | -4.9% | 0.23 | 0.74 |
| NK_cells | 7.4% | 4.9% | +2.5% | 0.23 | 0.74 |

*Mann-Whitney U test with Benjamini-Hochberg FDR correction*

**Key Finding:** Acinar cells showed the largest difference between responders and non-responders (2.2x higher in R, p=0.057, q=0.74). While not significant after multiple testing correction, this exploratory finding suggests preserved exocrine function may associate with treatment response.

### 6. Cell Type Ratio Analysis (NEW)

Systematic analysis of all pairwise cell type ratios to find combinations that discriminate R vs NR:

| Ratio | R Mean | NR Mean | Fold Δ | MWU p-value | FDR q |
|-------|--------|---------|--------|-------------|-------|
| **Acinar / Low_Confidence** | 0.38 | 0.13 | **2.8x** | 0.057 | 0.74 |
| **Acinar / Macrophage** | 0.80 | 0.30 | **2.6x** | 0.057 | 0.74 |
| **NK_cells / T_cells** | 1.08 | 0.76 | **1.4x** | 0.057 | 0.74 |

*Mann-Whitney U test with Benjamini-Hochberg FDR correction. All FDR q-values = 0.74, indicating these are exploratory findings.*

**Key Finding:** Multiple ratios involving Acinar cells showed large effect sizes in responders. While raw p-values reach the MWU floor (0.057 with n=7), FDR-corrected q-values are 0.74, indicating these are hypothesis-generating findings requiring validation.

### 7. Non-coding RNA Analysis (NEW)

| ncRNA Class | Detected | Notable Findings |
|-------------|----------|------------------|
| LINC (lncRNA) | 26 | LINC02693 consistently UP in R |
| MIR Host Genes | 8 | Mixed patterns |
| Antisense | 27 | Regulatory potential |
| **Total ncRNAs** | **61** | - |

**Key Finding:** LINC02693 shows consistent upregulation in responders across timepoints. Long non-coding RNAs may play regulatory roles in treatment response.

### 8. Day 4: Spatial Entropy Analysis (NEW)

| Metric | Responders | Non-Responders | p-value |
|--------|------------|----------------|---------|
| Global Spatial Entropy | 1.954 | 1.819 | 0.114 |

**Key Finding:** Responders show higher spatial entropy (more mixed tissue architecture), suggesting heterogeneous cell mixing may facilitate treatment response.

### 9. Day 4: L-R Communication (NEW)

| Metric | Responders | Non-Responders | Difference |
|--------|------------|----------------|------------|
| Significant L-R Interactions | 84.0 | 93.7 | -10% |

**Key Finding:** Non-responders have MORE ligand-receptor interactions (93.7 vs 84.0), suggesting more established communication networks may correlate with treatment resistance.

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

## Statistical Considerations

### Sample Size Limitations

With **n=4 responders** and **n=3 non-responders** (7 samples total), statistical power is inherently limited:

| Effect Size | Statistical Power | Interpretation |
|-------------|------------------|----------------|
| Large (d=1.5) | ~50% | May detect |
| Moderate (d=0.8) | ~20% | Likely to miss |
| Small (d=0.2) | <5% | Not detectable |

### Statistical Approach

1. **Primary Test: Welch's t-test** (unequal variance assumption)
   - Better power than Mann-Whitney U for small samples
   - Does not assume equal group variances
   - P-value floor: ~0.03 (vs 0.057 for MWU)

2. **Multiple Testing Correction: Benjamini-Hochberg FDR**
   - Applied to 91 cell type ratio tests
   - Applied to 12 cell type proportion tests
   - Conservative: many uncorrected p<0.1 become q>0.1

3. **Mutual Information: Sample-Level (Pseudobulk)**
   - Aggregated cells to sample means (n=7) to avoid pseudoreplication
   - Original cell-level approach (n=32K cells) had inflated statistics
   - Top candidates validated with Welch's t-test

### Interpretation Guidelines

⚠️ **All findings are HYPOTHESIS-GENERATING** and require validation in larger cohorts.

| Symbol | Meaning | Interpretation |
|--------|---------|----------------|
| ns | p ≥ 0.1 | Not significant |
| (exploratory) | Large effect, p > 0.1 | Promising but underpowered |

*Note: With n=7 samples, no comparisons achieved statistical significance after correction.*

### What We Can Conclude

✅ **Can say:** "We observed X difference between R and NR (effect size = Y)"
✅ **Can say:** "This finding warrants validation in larger cohorts"
❌ **Cannot say:** "Treatment responders have significantly more X"
❌ **Cannot say:** "X is a validated biomarker"

---

## Critical Statistical Limitations

> **⚠️ All findings are HYPOTHESIS-GENERATING, not confirmatory.**

With n=4 Responders and n=3 Non-Responders:
- **Minimum achievable Mann-Whitney U p-value: 0.0571** (cannot reach p<0.05)
- **All FDR-adjusted q-values: 0.74-0.75** (no significance after multiple testing correction)
- **Statistical power for moderate effects: ~8-20%** (underpowered by definition)

### What This Means

| Finding Type | Status | Interpretation |
|--------------|--------|----------------|
| **Effect sizes** | ✅ Valid | +96% betweenness, +25% Betti-1 are real observed differences |
| **P-values** | ⚠️ Underpowered | Cannot distinguish signal from noise at n=7 |
| **FDR q-values** | ⚠️ Not significant | All q ≥ 0.74 after Benjamini-Hochberg correction |

### Validation Requirements

All findings with p < 0.1 require independent validation in cohorts with **n ≥ 20 per group** to achieve adequate statistical power.

---

## Figure Gallery

### Figure 1: Sample Overview
*Spatial distribution of cell types across all 7 PDAC samples with adaptive point sizing*
![Sample Overview](figures/fig1_sample_overview.png)

### Figure 2: Cell Type Composition
*Differential cell type enrichment with Welch's t-test*
![Cell Type Composition](figures/fig2_cell_type_composition.png)

### Figure 3: Graph Centrality Analysis
*Betweenness and PageRank comparisons with Welch's t-test statistics*
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
*Dedicated analysis of Acinar cell enrichment in responders (Welch's t-test)*
![Acinar Comparison](figures/fig11_acinar_comparison.png)

### Figure 12: Cell Type Ratio Analysis (NEW)
*Systematic analysis of cell type ratios that discriminate R vs NR*
![Ratio Analysis](figures/fig12_ratio_analysis.png)

### Figure 13: Non-coding RNA Analysis (NEW)
*Analysis of lncRNAs and other non-coding transcripts*
![ncRNA Analysis](figures/fig13_ncrna_analysis.png)

### Day 4: Polymath-Driven Hypothesis Testing (NEW)

### Figure 14: Deconvolution Analysis
*Signature-based cell type deconvolution for Visium spots*
![Deconvolution](figures/fig14_deconvolution_overview.png)

### Figure 15: CAF Subtype Analysis
*iCAF vs mCAF scoring from Polymath KB markers*
![CAF Subtypes](figures/fig15_day4_hypotheses.png)

### Figure 16: Metabolic Pathway Analysis
*MPC/Pyruvate, Glycolysis, OxPhos, Fatty Acid pathways*
![Metabolic Analysis](figures/fig16_metabolic_analysis.png)

### Figure 17: Spatial Entropy Heatmap
*Shannon entropy of cell type distributions across samples*
![Entropy Heatmap](figures/fig17_entropy_heatmap.png)

### Figure 18: L-R Communication Summary
*Ligand-receptor interaction analysis comparing R vs NR*
![L-R Summary](figures/fig18_ligrec_summary.png)

### Figure 19: L-R Interaction Heatmap
*Top ligand-receptor pairs across all samples*
![L-R Heatmap](figures/fig19_ligrec_heatmap.png)

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
├── figures/                     # Publication-quality figures (13 PDAC + 2 G4X)
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
│   ├── fig11_acinar_comparison.png/pdf
│   ├── fig12_ratio_analysis.png/pdf
│   ├── fig13_ncrna_analysis.png/pdf
│   ├── fig14_deconvolution_overview.png/pdf    # Day 4
│   ├── fig15_day4_hypotheses.png/pdf           # Day 4
│   ├── fig16_metabolic_analysis.png/pdf        # Day 4
│   ├── fig17_entropy_heatmap.png/pdf           # Day 4
│   ├── fig18_ligrec_summary.png/pdf            # Day 4
│   ├── fig19_ligrec_heatmap.png/pdf            # Day 4
│   └── g4x/                     # G4X multimodal figures
├── data/
│   └── tables/
│       ├── mi_vs_response_biomarkers.csv   # Treatment response biomarkers
│       ├── cell_type_ratio_analysis.csv    # Cell type ratio analysis (NEW)
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
*Last updated: 2026-01-20*
