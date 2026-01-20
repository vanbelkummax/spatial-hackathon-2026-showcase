#!/usr/bin/env python3
"""
Day 3: Polymathic Analysis - Cross-Domain Algorithms
=====================================================

Apply algorithms from topology, graph theory, and information theory
to spatial transcriptomics data, as identified through Polymath KB.

Methods:
1. Graph Centrality (from Graph Theory) - Hub cell identification
2. Persistent Homology (from Topology) - Tissue architecture
3. Mutual Information (from Information Theory) - Gene selection

Author: Max Van Belkum
Date: 2026-01-19
"""

import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd
import networkx as nx
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import warnings
import logging

warnings.filterwarnings('ignore')

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
logger = logging.getLogger(__name__)

# Check for optional TDA imports
try:
    from gtda.homology import VietorisRipsPersistence
    from gtda.diagrams import BettiCurve
    HAS_TDA = True
except ImportError:
    HAS_TDA = False
    logger.warning("giotto-tda not available. Install with: pip install giotto-tda")

# =============================================================================
# Configuration
# =============================================================================

PROJECT_ROOT = Path(__file__).parent.parent
OUTPUT_DIR = PROJECT_ROOT / "outputs"
ADATA_DIR = OUTPUT_DIR / "adata" / "annotated"
POLYMATHIC_DIR = OUTPUT_DIR / "adata" / "polymathic"
FIG_DIR = OUTPUT_DIR / "figures" / "polymathic"
TABLE_DIR = OUTPUT_DIR / "tables"

POLYMATHIC_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Sample metadata (PDAC only)
PDAC_METADATA = {
    "YP03A": {"patient": "YP03", "timepoint": "Pre", "response": "NR"},
    "YP03C": {"patient": "YP03", "timepoint": "Post", "response": "NR"},
    "YP04C": {"patient": "YP04", "timepoint": "Post", "response": "NR"},  # YP04A excluded
    "YP12A": {"patient": "YP12", "timepoint": "Pre", "response": "R"},
    "YP12C": {"patient": "YP12", "timepoint": "Post", "response": "R"},
    "YP15A": {"patient": "YP15", "timepoint": "Pre", "response": "R"},
    "YP15C": {"patient": "YP15", "timepoint": "Post", "response": "R"},
}


# =============================================================================
# Graph Centrality Analysis (from Graph Theory)
# =============================================================================

def compute_graph_centrality(adata, n_neighs: int = 6) -> Dict[str, float]:
    """
    Compute hub cell statistics using graph centrality measures.

    Cross-domain application: Graph centrality algorithms from social network
    analysis applied to spatial cell interaction networks.

    Returns dict with centrality statistics.
    """
    logger.info("Computing spatial graph centrality...")

    # Ensure spatial graph exists
    if 'spatial_connectivities' not in adata.obsp:
        sq.gr.spatial_neighbors(adata, n_neighs=n_neighs)

    # Convert to networkx
    G = nx.from_scipy_sparse_array(adata.obsp['spatial_connectivities'])

    # Compute centrality measures
    logger.info("  - Betweenness centrality...")
    betweenness = nx.betweenness_centrality(G)

    logger.info("  - PageRank...")
    pagerank = nx.pagerank(G)

    logger.info("  - Degree centrality...")
    degree = nx.degree_centrality(G)

    # Store in adata
    adata.obs['betweenness'] = [betweenness.get(i, 0) for i in range(adata.n_obs)]
    adata.obs['pagerank'] = [pagerank.get(i, 0) for i in range(adata.n_obs)]
    adata.obs['degree_centrality'] = [degree.get(i, 0) for i in range(adata.n_obs)]

    # Identify hub cells (top 5% by pagerank)
    threshold = np.percentile(adata.obs['pagerank'], 95)
    adata.obs['is_hub_cell'] = adata.obs['pagerank'] > threshold

    # Calculate statistics
    stats = {
        'n_cells': adata.n_obs,
        'n_hub_cells': adata.obs['is_hub_cell'].sum(),
        'hub_fraction': adata.obs['is_hub_cell'].mean(),
        'mean_betweenness': adata.obs['betweenness'].mean(),
        'max_betweenness': adata.obs['betweenness'].max(),
        'mean_pagerank': adata.obs['pagerank'].mean(),
        'max_pagerank': adata.obs['pagerank'].max(),
        'mean_degree': adata.obs['degree_centrality'].mean(),
    }

    # Hub cell type distribution
    if 'cell_type' in adata.obs.columns:
        hub_cells = adata[adata.obs['is_hub_cell']]
        hub_type_counts = hub_cells.obs['cell_type'].value_counts()
        stats['dominant_hub_type'] = hub_type_counts.index[0] if len(hub_type_counts) > 0 else 'Unknown'

    return stats


# =============================================================================
# Persistent Homology Analysis (from Topology)
# =============================================================================

def compute_topology_features(adata, max_cells: int = 5000) -> Optional[Dict[str, float]]:
    """
    Compute persistent homology features for tissue architecture.

    Cross-domain application: Topological data analysis from algebraic topology
    applied to quantify tissue structure complexity.

    H0 = connected components (tissue clusters)
    H1 = holes/voids (tissue gaps/boundaries)
    """
    if not HAS_TDA:
        logger.warning("Skipping topology - giotto-tda not installed")
        return None

    logger.info("Computing persistent homology...")

    coords = adata.obsm['spatial']

    # Subsample if too large (TDA is O(n^3))
    if len(coords) > max_cells:
        idx = np.random.choice(len(coords), max_cells, replace=False)
        coords = coords[idx]
        logger.info(f"  - Subsampled to {max_cells} cells")

    # Normalize coordinates
    coords = (coords - coords.min(axis=0)) / (coords.max(axis=0) - coords.min(axis=0) + 1e-6)

    # Compute Vietoris-Rips persistence
    VR = VietorisRipsPersistence(
        metric='euclidean',
        homology_dimensions=[0, 1],
        n_jobs=-1
    )
    diagrams = VR.fit_transform([coords])

    # Compute Betti curves
    betti = BettiCurve()
    features = betti.fit_transform(diagrams)

    stats = {
        'betti_0_mean': features[0, :, 0].mean(),  # Connected components
        'betti_0_max': features[0, :, 0].max(),
        'betti_1_mean': features[0, :, 1].mean(),  # Holes/voids
        'betti_1_max': features[0, :, 1].max(),
        'topology_complexity': features[0, :, 0].mean() + features[0, :, 1].mean(),
    }

    return stats


# =============================================================================
# Mutual Information Analysis (from Information Theory)
# =============================================================================

def compute_mutual_information_genes(adata, target_col: str = 'cell_type', top_n: int = 50) -> pd.DataFrame:
    """
    Compute mutual information between genes and a target variable.

    Cross-domain application: Information theory from signal processing/statistics
    applied to identify non-linearly associated genes.

    MI captures dependencies that correlation misses.
    """
    from sklearn.feature_selection import mutual_info_classif

    logger.info(f"Computing mutual information for genes vs {target_col}...")

    if target_col not in adata.obs.columns:
        logger.warning(f"  - {target_col} not found, skipping MI analysis")
        return pd.DataFrame()

    # Get expression matrix and target
    X = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X
    y = adata.obs[target_col].astype('category').cat.codes.values

    # Remove cells with unknown labels
    mask = y >= 0
    X = X[mask]
    y = y[mask]

    # Compute MI (this can be slow for many genes)
    n_genes = min(5000, X.shape[1])  # Limit for speed
    if X.shape[1] > n_genes:
        # Use most variable genes
        gene_var = np.var(X, axis=0)
        top_var_idx = np.argsort(gene_var)[-n_genes:]
        X = X[:, top_var_idx]
        gene_names = adata.var_names[top_var_idx]
    else:
        gene_names = adata.var_names

    logger.info(f"  - Computing MI for {len(gene_names)} genes...")
    mi_scores = mutual_info_classif(X, y, random_state=42)

    # Create results DataFrame
    mi_df = pd.DataFrame({
        'gene': gene_names,
        'mutual_information': mi_scores
    }).sort_values('mutual_information', ascending=False)

    return mi_df.head(top_n)


# =============================================================================
# Main Analysis
# =============================================================================

def analyze_sample(sample_name: str, adata_path: Path) -> Dict:
    """Run full polymathic analysis on a single sample."""
    logger.info(f"\n{'='*60}")
    logger.info(f"Analyzing {sample_name}")
    logger.info(f"{'='*60}")

    adata = sc.read_h5ad(adata_path)
    logger.info(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")

    # Get metadata
    meta = PDAC_METADATA.get(sample_name, {})

    results = {
        'sample': sample_name,
        'n_cells': adata.n_obs,
        'n_genes': adata.n_vars,
        **meta
    }

    # 1. Graph Centrality
    centrality_stats = compute_graph_centrality(adata)
    results.update({f'centrality_{k}': v for k, v in centrality_stats.items()})

    # 2. Topology Features
    topology_stats = compute_topology_features(adata)
    if topology_stats:
        results.update({f'topology_{k}': v for k, v in topology_stats.items()})

    # 3. Mutual Information (only for cell type)
    mi_df = compute_mutual_information_genes(adata, 'cell_type')
    if len(mi_df) > 0:
        results['top_mi_gene'] = mi_df.iloc[0]['gene']
        results['top_mi_score'] = mi_df.iloc[0]['mutual_information']

        # Save MI results
        mi_path = TABLE_DIR / f"{sample_name}_mi_genes.csv"
        mi_df.to_csv(mi_path, index=False)

    # Save updated adata
    out_path = POLYMATHIC_DIR / f"{sample_name}_polymathic.h5ad"
    adata.write_h5ad(out_path)
    logger.info(f"Saved to {out_path}")

    return results


def main():
    """Run polymathic analysis on all PDAC samples."""
    logger.info("="*60)
    logger.info("DAY 3: POLYMATHIC ANALYSIS")
    logger.info("Cross-domain algorithms for spatial transcriptomics")
    logger.info("="*60)

    all_results = []

    # Analyze PDAC samples
    for sample_name in PDAC_METADATA.keys():
        adata_path = ADATA_DIR / f"{sample_name}_annotated.h5ad"

        if not adata_path.exists():
            logger.warning(f"Skipping {sample_name} - file not found")
            continue

        try:
            results = analyze_sample(sample_name, adata_path)
            all_results.append(results)
        except Exception as e:
            logger.error(f"Error analyzing {sample_name}: {e}")
            continue

    # Combine results
    if all_results:
        results_df = pd.DataFrame(all_results)
        results_path = TABLE_DIR / "polymathic_analysis_results.csv"
        results_df.to_csv(results_path, index=False)
        logger.info(f"\nSaved combined results to {results_path}")

        # Print summary comparison
        logger.info("\n" + "="*60)
        logger.info("R vs NR COMPARISON (Polymathic Features)")
        logger.info("="*60)

        r_samples = results_df[results_df['response'] == 'R']
        nr_samples = results_df[results_df['response'] == 'NR']

        metrics = ['centrality_hub_fraction', 'centrality_mean_pagerank',
                   'centrality_mean_betweenness']

        if HAS_TDA:
            metrics.extend(['topology_betti_0_mean', 'topology_betti_1_mean',
                           'topology_complexity'])

        for metric in metrics:
            if metric in results_df.columns:
                r_mean = r_samples[metric].mean()
                nr_mean = nr_samples[metric].mean()
                diff_pct = ((r_mean - nr_mean) / (nr_mean + 1e-6)) * 100
                logger.info(f"  {metric:35s}: R={r_mean:.4f}, NR={nr_mean:.4f} ({diff_pct:+.1f}%)")

    logger.info("\n" + "="*60)
    logger.info("ANALYSIS COMPLETE")
    logger.info("="*60)


if __name__ == '__main__':
    main()
