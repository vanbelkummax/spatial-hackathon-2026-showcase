# Next Session: Final Audit & Expansion Plan

**Date:** 2026-01-20+
**Goal:** Complete senior review fixes, expand analysis to all datasets, prepare for Choi samples

---

## 🚨 Priority 1: Remaining Code Fixes

### 1.1 Verify TDA Fix is Active

**Status:** Code written but needs verification

```bash
# Run the figure generation and check Figure 5 output
cd /home/user/spatial-hackathon-2026
source ~/miniforge3/etc/profile.d/conda.sh && conda activate enact

# Generate just Figure 5 and inspect
python -c "
from scripts.09_generate_showcase_figures_fixed import fig5_persistent_homology
fig5_persistent_homology()
"

# Visually inspect the PDF - should show scatter plots with H0 (red) and H1 (blue) dots
# NOT error messages or blank panels
```

**Verification Criteria:**
- [ ] All 7 sample panels show persistence diagrams (not "Error" text)
- [ ] Diagonal line (birth=death) visible
- [ ] H0 (red) and H1 (blue) points distinguishable
- [ ] PDF file size > 40KB (broken version was ~22KB)

---

### 1.2 Fix MI Logic - Critical Scientific Error

**Current Problem:** MI computed against `cell_type` but presented as "Treatment Response Biomarkers"

**Two Options:**

#### Option A: Rename and Clarify (Minimal Change)
- Rename Figure 7 title: "Top Cell Type Discriminative Genes" → "Cell Type Marker Stability Across Response Groups"
- Update FINDINGS.md to clarify these are NOT response biomarkers
- Keep DE genes (MS4A2, NLRP7, CCZ1) as the actual response biomarkers

#### Option B: Compute Cross-Sample MI Against Response (Recommended)

```python
# NEW FUNCTION for scripts/08_polymathic_analysis.py

def compute_cross_sample_mi_response(all_adatas: Dict[str, anndata.AnnData],
                                      metadata: Dict, top_n: int = 50) -> pd.DataFrame:
    """
    Compute mutual information between genes and TREATMENT RESPONSE.

    This is the correct approach for finding treatment response biomarkers,
    as it compares gene expression ACROSS samples (pooled) against R/NR labels.
    """
    from sklearn.feature_selection import mutual_info_classif
    import anndata as ad

    # Merge all samples
    adatas_list = []
    for sample, adata in all_adatas.items():
        adata_copy = adata.copy()
        adata_copy.obs['sample'] = sample
        adata_copy.obs['response'] = metadata[sample]['response']
        adatas_list.append(adata_copy)

    merged = ad.concat(adatas_list, join='inner')

    # Get expression and response labels
    X = merged.X.toarray() if hasattr(merged.X, 'toarray') else merged.X
    y = (merged.obs['response'] == 'R').astype(int).values  # 1=R, 0=NR

    # Compute MI
    mi_scores = mutual_info_classif(X, y, random_state=42)

    return pd.DataFrame({
        'gene': merged.var_names,
        'mi_vs_response': mi_scores
    }).sort_values('mi_vs_response', ascending=False).head(top_n)
```

**Query Polymath for MI implementations:**
```bash
cd /home/user/polymath-v4
python scripts/q.py "mutual information feature selection treatment response" -n 10
python scripts/q.py "information theoretic gene selection cancer" --funcs -n 10
```

---

### 1.3 Upgrade TDA Sampling: MaxMin Landmark Selection

**Current Problem:** Random subsampling destroys topological features

**Solution:** Use MaxMin (farthest point) sampling to preserve shape

**Query Polymath for implementations:**
```bash
cd /home/user/polymath-v4

# Search for landmark sampling methods
python scripts/q.py "landmark selection maxmin farthest point sampling" -n 10
python scripts/q.py "topology preserving subsampling point cloud" --funcs -n 10
python scripts/algo.py "farthest point sampling"
python scripts/algo.py "landmark selection"

# Search for TDA best practices
python scripts/q.py "persistent homology subsampling strategies" -n 10
```

**Implementation (add to 08_polymathic_analysis.py):**

```python
def maxmin_subsample(coords: np.ndarray, n_landmarks: int, seed: int = 42) -> np.ndarray:
    """
    MaxMin (farthest point) landmark selection.

    Preserves topological features by selecting points that maximize coverage
    of the point cloud shape, rather than random selection which can break loops.

    Algorithm:
    1. Start with random point
    2. Repeatedly add the point farthest from all current landmarks
    3. Continue until n_landmarks reached

    References:
    - De Silva & Carlsson (2004) "Topological estimation using witness complexes"
    - Polymath KB: Search "landmark selection topology"
    """
    rng = np.random.default_rng(seed)
    n = len(coords)

    if n <= n_landmarks:
        return np.arange(n)

    # Start with random point
    landmarks = [rng.integers(n)]

    # Distance to nearest landmark for each point
    min_dists = np.full(n, np.inf)

    for _ in range(n_landmarks - 1):
        # Update distances to nearest landmark
        last_landmark = coords[landmarks[-1]]
        dists_to_last = np.linalg.norm(coords - last_landmark, axis=1)
        min_dists = np.minimum(min_dists, dists_to_last)

        # Add farthest point
        next_landmark = np.argmax(min_dists)
        landmarks.append(next_landmark)

    return np.array(landmarks)


def compute_topology_features(adata, max_cells: int = 5000) -> Optional[Dict[str, float]]:
    """Updated to use MaxMin sampling instead of random."""
    # ... existing code ...

    if len(coords) > max_cells:
        # OLD: idx = rng.choice(len(coords), max_cells, replace=False)
        # NEW: Topology-preserving landmark selection
        idx = maxmin_subsample(coords, max_cells, seed=42)
        coords = coords[idx]
        logger.info(f"  - MaxMin subsampled to {max_cells} cells (preserves topology)")

    # ... rest of function ...
```

---

## 🔬 Priority 2: Expand Analysis to ALL Datasets

### 2.1 G4X Multimodal Data (H&E + RNA + Protein)

**Location:** `~/data/hackathon/multimodal/` (35GB)

**What's in G4X:**
- H&E images
- Spatial gene expression
- Spatial protein expression (antibody panel)

**Analysis Plan:**

```bash
# Explore G4X data structure
ls -la ~/data/hackathon/multimodal/

# Query Polymath for multimodal integration methods
cd /home/user/polymath-v4
python scripts/q.py "multimodal spatial transcriptomics proteomics integration" -n 15
python scripts/q.py "H&E RNA protein joint embedding" -n 10
python scripts/q.py "CITE-seq spatial integration" --repos -n 10

# Cross-domain algorithms for multimodal data
python scripts/algo.py "canonical correlation analysis"
python scripts/algo.py "multi-view learning"
python scripts/algo.py --bridges  # Find cross-domain transfers
```

**New Analysis Scripts Needed:**

1. `10_g4x_qc_preprocessing.py` - QC the multimodal data
2. `11_g4x_protein_analysis.py` - Analyze protein channels
3. `12_g4x_multimodal_integration.py` - Joint RNA+Protein analysis
4. `13_g4x_polymathic.py` - Apply graph centrality + TDA to multimodal

**G4X-Specific Questions:**
- Do hub cells in RNA space also have high protein expression?
- Does TDA on protein coordinates match TDA on RNA coordinates?
- Can we find RNA-Protein discordance (post-transcriptional regulation)?

---

### 2.2 Complete PDAC Analysis (All Choi Samples)

**Current:** 7 samples analyzed (YP03A, YP03C, YP04C, YP12A, YP12C, YP15A, YP15C)
**Excluded:** YP04A (134 cells - too few after QC)

**When Remaining Choi Samples Arrive:**

Expected additional samples from Choi lab:
- More patients with known response status
- Potentially different treatment timepoints
- Higher quality samples (better cell yield)

**Protocol for New Samples:**

```bash
# 1. QC Check (before full analysis)
python scripts/01_qc_check.py --input /path/to/new/samples --min-cells 500

# 2. If passes QC, run full pipeline
python scripts/02_preprocessing.py --samples NEW_SAMPLE_IDS
python scripts/03_clustering.py --samples NEW_SAMPLE_IDS
python scripts/04_annotation.py --samples NEW_SAMPLE_IDS
python scripts/05_spatial.py --samples NEW_SAMPLE_IDS
python scripts/06_de.py --samples NEW_SAMPLE_IDS

# 3. Run polymathic analysis
python scripts/08_polymathic_analysis.py --samples NEW_SAMPLE_IDS

# 4. Regenerate all figures with combined data
python scripts/09_generate_showcase_figures_fixed.py
```

**Statistical Power Considerations:**
- Current: n=4 patients (2 R, 2 NR) → weak statistics
- With 10+ patients: Can use proper statistical tests
- With 20+ patients: Can train simple ML classifiers

---

## 📚 Priority 3: Polymath-Driven Method Expansion

### 3.1 Additional Cross-Domain Algorithms to Apply

**Query Polymath for more algorithms:**

```bash
cd /home/user/polymath-v4

# Traditional spatial biology methods we might have missed
python scripts/q.py "ripley K function spatial point pattern" -n 10
python scripts/q.py "spatial autocorrelation local indicators LISA" -n 10
python scripts/q.py "co-localization analysis spatial transcriptomics" -n 10

# Cross-domain algorithms with spatial biology potential
python scripts/algo.py --domain physics --spatial
python scripts/algo.py --domain information_theory --spatial
python scripts/algo.py --domain control_theory --spatial

# Specific algorithms to explore
python scripts/algo.py "ising model"           # Physics → cell interactions
python scripts/algo.py "markov random field"   # Statistics → spatial patterns
python scripts/algo.py "optimal transport"     # Math → cell distribution comparison
python scripts/algo.py "spectral clustering"   # Graph theory → spatial communities
```

### 3.2 New Methods to Implement

| Method | Origin | Application | Polymath Query |
|--------|--------|-------------|----------------|
| **Ripley's K** | Point process stats | Spatial clustering significance | `ripley K spatial` |
| **Wasserstein Distance** | Optimal transport | Compare R vs NR cell distributions | `wasserstein distance spatial` |
| **Spatial Entropy** | Information theory | Tissue heterogeneity quantification | `spatial entropy single cell` |
| **Ising Model** | Statistical physics | Cell-cell interaction energy | `ising model cell` |
| **Spectral Clustering** | Graph theory | Identify spatial cell communities | `spectral clustering spatial` |

### 3.3 Code Repositories to Ingest

**Query for relevant repos:**

```bash
cd /home/user/polymath-v4

# TDA repos
python scripts/q.py "giotto-tda spatial" --repos
python scripts/q.py "scikit-tda persistent homology" --repos

# Spatial statistics repos
python scripts/q.py "squidpy spatial statistics" --repos
python scripts/q.py "spatialdata analysis" --repos

# Multimodal integration repos
python scripts/q.py "muon multimodal" --repos
python scripts/q.py "totalVI protein RNA" --repos
```

**Priority Repos to Ingest:**

1. `scikit-tda/scikit-tda` - Python TDA toolkit
2. `giotto-ai/giotto-tda` - Already using, ensure full indexing
3. `scverse/spatialdata` - Unified spatial data structures
4. `scverse/muon` - Multimodal single-cell analysis
5. `theislab/cell2location` - Spatial deconvolution

---

## 📋 Session Checklist

### Before Starting
- [ ] Activate enact environment: `conda activate enact`
- [ ] Verify Polymath is running: `python scripts/system_report.py --quick`
- [ ] Check disk space: `python scripts/fs.py stats`

### Code Fixes (Priority 1)
- [ ] Verify Figure 5 TDA fix is working (visual inspection)
- [ ] Implement Option B: Cross-sample MI against response variable
- [ ] Implement MaxMin landmark selection for TDA
- [ ] Run full figure regeneration
- [ ] Push fixes to GitHub

### Analysis Expansion (Priority 2)
- [ ] Explore G4X multimodal data structure
- [ ] Write G4X preprocessing script
- [ ] Apply polymathic analysis to G4X
- [ ] Document protocol for incoming Choi samples

### Polymath Enhancement (Priority 3)
- [ ] Query for additional cross-domain algorithms
- [ ] Implement at least 2 new methods (Ripley's K, Wasserstein)
- [ ] Ingest priority repos (scikit-tda, spatialdata)
- [ ] Update POLYMATH_INTEGRATION.md with new methods

### Final Audit
- [ ] All figures render without errors
- [ ] MI analysis correctly labeled (cell type vs response)
- [ ] MaxMin sampling verified with before/after topology comparison
- [ ] G4X analysis added to showcase
- [ ] README updated with new methods

---

## 🎯 Success Metrics

| Metric | Current | Target |
|--------|---------|--------|
| Samples analyzed | 7 PDAC | 7 PDAC + G4X |
| Cross-domain algorithms | 3 (centrality, TDA, MI) | 6+ |
| Figure 5 errors | Fixed | 0 errors |
| MI analysis | Against cell_type | Against response (new) |
| TDA subsampling | Random | MaxMin landmark |
| Repos in Polymath | 1,881 | 1,886+ |

---

## 📝 Commands Reference

```bash
# Environment
conda activate enact
cd /home/user/polymath-v4

# Polymath queries
python scripts/q.py "your query" -n 10           # Papers
python scripts/q.py "your query" --repos         # Code repos
python scripts/q.py "your query" --funcs         # Functions
python scripts/algo.py "algorithm name"          # Algorithm registry
python scripts/algo.py --domain DOMAIN --spatial # Spatial algorithms

# Analysis
cd /home/user/spatial-hackathon-2026
python scripts/08_polymathic_analysis.py
python scripts/09_generate_showcase_figures_fixed.py

# Git
cd /home/user/spatial-hackathon-2026-showcase
git add -A && git commit -m "message" && git push origin master
```

---

*Created: 2026-01-19*
*For: Spatial Biology Hackathon 2026 Final Audit*
