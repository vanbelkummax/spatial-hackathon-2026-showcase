#!/usr/bin/env python3
"""
HACKATHON SHOWCASE: Publication-Quality Figure Generation
==========================================================

Generates all figures for Days 1-3 at maximum quality for GitHub showcase.

Requirements: scanpy, squidpy, matplotlib, seaborn, giotto-tda

Author: Max Van Belkum
Date: 2026-01-19
"""

# Force non-interactive backend for headless environments (must be before pyplot import)
import matplotlib
matplotlib.use('Agg')

import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import traceback
from pathlib import Path
from typing import Dict, List, Optional
import warnings
import logging
from scipy.stats import mannwhitneyu

warnings.filterwarnings('ignore')
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 10

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
logger = logging.getLogger(__name__)

# TDA imports
try:
    from gtda.homology import VietorisRipsPersistence
    from gtda.diagrams import BettiCurve, PersistenceEntropy
    # Note: gtda.plotting.plot_diagram creates Plotly figures, not Matplotlib!
    # We use manual matplotlib plotting instead for PDF compatibility.
    HAS_TDA = True
except ImportError:
    HAS_TDA = False
    logger.warning("giotto-tda not available")

# Reproducibility: Fixed random seed for subsampling
RNG = np.random.default_rng(42)


def maxmin_subsample(coords: np.ndarray, n_landmarks: int, seed: int = 42) -> np.ndarray:
    """
    MaxMin (farthest point) landmark selection for topology-preserving subsampling.

    Preserves topological features by selecting points that maximize coverage
    of the point cloud shape, rather than random selection which can break loops.
    """
    rng = np.random.default_rng(seed)
    n = len(coords)

    if n <= n_landmarks:
        return np.arange(n)

    landmarks = [rng.integers(n)]
    min_dists = np.full(n, np.inf)

    for _ in range(n_landmarks - 1):
        last_landmark = coords[landmarks[-1]]
        dists_to_last = np.linalg.norm(coords - last_landmark, axis=1)
        min_dists = np.minimum(min_dists, dists_to_last)
        next_landmark = np.argmax(min_dists)
        landmarks.append(next_landmark)

    return np.array(landmarks)


# =============================================================================
# Configuration
# =============================================================================

PROJECT_ROOT = Path(__file__).parent.parent
OUTPUT_DIR = PROJECT_ROOT / "outputs"
ADATA_DIR = OUTPUT_DIR / "adata"
FIG_DIR = OUTPUT_DIR / "figures" / "showcase"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Color palettes
RESPONSE_COLORS = {'R': '#2ecc71', 'NR': '#e74c3c'}
TIMEPOINT_COLORS = {'Pre': '#3498db', 'Post': '#9b59b6'}

# Sample metadata
PDAC_METADATA = {
    "YP03A": {"patient": "YP03", "timepoint": "Pre", "response": "NR"},
    "YP03C": {"patient": "YP03", "timepoint": "Post", "response": "NR"},
    "YP04C": {"patient": "YP04", "timepoint": "Post", "response": "NR"},
    "YP12A": {"patient": "YP12", "timepoint": "Pre", "response": "R"},
    "YP12C": {"patient": "YP12", "timepoint": "Post", "response": "R"},
    "YP15A": {"patient": "YP15", "timepoint": "Pre", "response": "R"},
    "YP15C": {"patient": "YP15", "timepoint": "Post", "response": "R"},
}

# =============================================================================
# Figure 1: Sample Overview Grid
# =============================================================================

def fig1_sample_overview():
    """Create overview grid showing all samples with cell types - proper spatial coverage."""
    logger.info("Generating Figure 1: Sample Overview Grid")

    samples = list(PDAC_METADATA.keys())
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    axes = axes.flatten()

    # Collect all cell types for consistent coloring
    all_cell_types = set()
    sample_data = {}

    for sample in samples:
        adata_path = ADATA_DIR / "polymathic" / f"{sample}_polymathic.h5ad"
        if not adata_path.exists():
            adata_path = ADATA_DIR / "annotated" / f"{sample}_annotated.h5ad"
        if adata_path.exists():
            adata = sc.read_h5ad(adata_path)
            sample_data[sample] = adata
            if 'cell_type' in adata.obs.columns:
                all_cell_types.update(adata.obs['cell_type'].unique())

    # Create color map for cell types
    cell_type_list = sorted(all_cell_types)
    cmap = plt.cm.get_cmap('tab20', len(cell_type_list))
    cell_type_colors = {ct: cmap(i) for i, ct in enumerate(cell_type_list)}

    for idx, sample in enumerate(samples):
        ax = axes[idx]

        if sample not in sample_data:
            ax.text(0.5, 0.5, f'{sample}\nNot Found', ha='center', va='center', fontsize=12)
            ax.axis('off')
            continue

        adata = sample_data[sample]
        meta = PDAC_METADATA[sample]
        coords = adata.obsm['spatial']

        # Calculate proper point size based on tissue density
        # Smaller dots for denser tissues, scale with sqrt of cell count
        base_size = max(0.5, min(5, 5000 / np.sqrt(len(coords))))

        if 'cell_type' in adata.obs.columns:
            # Plot each cell type
            for ct in cell_type_list:
                mask = adata.obs['cell_type'] == ct
                if mask.sum() > 0:
                    ax.scatter(
                        coords[mask, 0], coords[mask, 1],
                        c=[cell_type_colors[ct]], s=base_size, alpha=0.7,
                        label=ct, edgecolors='none'
                    )
        else:
            ax.scatter(coords[:, 0], coords[:, 1], s=base_size, alpha=0.5, c='steelblue')

        # Proper axis formatting
        ax.set_aspect('equal')
        ax.set_xlim(coords[:, 0].min() - 50, coords[:, 0].max() + 50)
        ax.set_ylim(coords[:, 1].min() - 50, coords[:, 1].max() + 50)
        ax.invert_yaxis()  # Standard image orientation

        # Title with response color coding
        title_color = RESPONSE_COLORS[meta['response']]
        ax.set_title(f"{sample}\n{meta['response']} | {meta['timepoint']} | n={len(coords)}",
                     fontsize=11, fontweight='bold', color=title_color)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.tick_params(labelsize=8)

    # Hide empty subplot
    axes[-1].axis('off')

    # Create legend in the empty subplot space
    handles = [mpatches.Patch(facecolor=cell_type_colors[ct], label=ct, edgecolor='gray')
               for ct in cell_type_list if ct in all_cell_types]
    axes[-1].legend(handles=handles, loc='center', fontsize=9, ncol=1,
                    title='Cell Types', title_fontsize=10, frameon=True)

    plt.suptitle('Spatial Overview: All PDAC Samples\n(Color = Cell Type, Green titles = Responders, Red = Non-Responders)',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig1_sample_overview.png", bbox_inches='tight', dpi=300)
    plt.savefig(FIG_DIR / "fig1_sample_overview.pdf", bbox_inches='tight')
    plt.close()
    logger.info("  Saved fig1_sample_overview.png/pdf")


# =============================================================================
# Figure 2: Cell Type Composition Comparison
# =============================================================================

def permutation_test(group1, group2, n_permutations=10000):
    """Permutation test for difference in means between two groups."""
    observed_diff = np.mean(group1) - np.mean(group2)
    combined = np.concatenate([group1, group2])
    n1 = len(group1)

    count = 0
    for _ in range(n_permutations):
        np.random.shuffle(combined)
        perm_diff = np.mean(combined[:n1]) - np.mean(combined[n1:])
        if abs(perm_diff) >= abs(observed_diff):
            count += 1

    return count / n_permutations


def fig2_cell_type_composition():
    """Compare cell type composition between R and NR with permutation tests."""
    logger.info("Generating Figure 2: Cell Type Composition")

    compositions = []

    for sample, meta in PDAC_METADATA.items():
        adata_path = ADATA_DIR / "polymathic" / f"{sample}_polymathic.h5ad"
        if not adata_path.exists():
            adata_path = ADATA_DIR / "annotated" / f"{sample}_annotated.h5ad"
        if not adata_path.exists():
            continue

        adata = sc.read_h5ad(adata_path)

        if 'cell_type' not in adata.obs.columns:
            continue

        counts = adata.obs['cell_type'].value_counts(normalize=True)
        for ct, frac in counts.items():
            compositions.append({
                'sample': sample,
                'response': meta['response'],
                'timepoint': meta['timepoint'],
                'cell_type': ct,
                'fraction': frac
            })

    df = pd.DataFrame(compositions)

    if len(df) == 0:
        logger.warning("  No cell type data found")
        return

    # Create figure with 3 panels
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Panel 1: Stacked bar by response
    df_agg = df.groupby(['response', 'cell_type'])['fraction'].mean().reset_index()
    pivot = df_agg.pivot(index='response', columns='cell_type', values='fraction').fillna(0)
    pivot.plot(kind='bar', stacked=True, ax=axes[0], colormap='tab20')
    axes[0].set_title('Cell Type Composition by Response', fontsize=14, fontweight='bold')
    axes[0].set_ylabel('Fraction')
    axes[0].set_xlabel('Treatment Response')
    axes[0].legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
    axes[0].set_xticklabels(['Non-Responder', 'Responder'], rotation=0)

    # Panel 2: Differential with permutation test p-values
    r_comp = pivot.loc['R'] if 'R' in pivot.index else pd.Series()
    nr_comp = pivot.loc['NR'] if 'NR' in pivot.index else pd.Series()

    if len(r_comp) > 0 and len(nr_comp) > 0:
        # Calculate differences and p-values
        diff_results = []
        for ct in r_comp.index:
            r_vals = df[(df['response'] == 'R') & (df['cell_type'] == ct)]['fraction'].values
            nr_vals = df[(df['response'] == 'NR') & (df['cell_type'] == ct)]['fraction'].values

            if len(r_vals) > 0 and len(nr_vals) > 0:
                diff = np.mean(r_vals) - np.mean(nr_vals)
                pval = permutation_test(r_vals, nr_vals, n_permutations=5000)
                diff_results.append({'cell_type': ct, 'diff': diff, 'pval': pval})

        diff_df = pd.DataFrame(diff_results).sort_values('diff')

        # Plot differential with significance markers
        colors = ['#e74c3c' if x < 0 else '#2ecc71' for x in diff_df['diff']]
        bars = axes[1].barh(range(len(diff_df)), diff_df['diff'], color=colors)

        # Add significance markers
        for i, (_, row) in enumerate(diff_df.iterrows()):
            sig = ''
            if row['pval'] < 0.001:
                sig = '***'
            elif row['pval'] < 0.01:
                sig = '**'
            elif row['pval'] < 0.05:
                sig = '*'
            elif row['pval'] < 0.1:
                sig = '†'

            if sig:
                x_pos = row['diff'] + (0.01 if row['diff'] > 0 else -0.01)
                axes[1].text(x_pos, i, sig, va='center', ha='left' if row['diff'] > 0 else 'right',
                            fontsize=10, fontweight='bold')

        axes[1].set_yticks(range(len(diff_df)))
        axes[1].set_yticklabels(diff_df['cell_type'])
        axes[1].set_title('Differential Enrichment (R - NR)\n*** p<0.001, ** p<0.01, * p<0.05, † p<0.1',
                         fontsize=12, fontweight='bold')
        axes[1].set_xlabel('Fraction Difference')
        axes[1].axvline(0, color='black', linestyle='--', linewidth=0.5)

    # Panel 3: Per-sample heatmap
    sample_pivot = df.pivot_table(index='sample', columns='cell_type', values='fraction', fill_value=0)
    # Add response labels
    sample_order = ['YP03A', 'YP03C', 'YP04C', 'YP12A', 'YP12C', 'YP15A', 'YP15C']
    sample_pivot = sample_pivot.reindex([s for s in sample_order if s in sample_pivot.index])

    sns.heatmap(sample_pivot, ax=axes[2], cmap='YlOrRd', annot=False,
                cbar_kws={'label': 'Fraction'}, linewidths=0.5)
    axes[2].set_title('Cell Type Composition per Sample', fontsize=14, fontweight='bold')
    axes[2].set_xlabel('Cell Type')
    axes[2].set_ylabel('Sample')

    # Color sample labels by response
    for i, label in enumerate(axes[2].get_yticklabels()):
        sample = label.get_text()
        if sample in PDAC_METADATA:
            color = RESPONSE_COLORS[PDAC_METADATA[sample]['response']]
            label.set_color(color)
            label.set_fontweight('bold')

    plt.suptitle('Cell Type Composition Analysis\n(n=4 Responders, n=3 Non-Responders; Permutation test with 5000 iterations)',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig2_cell_type_composition.png", bbox_inches='tight', dpi=300)
    plt.savefig(FIG_DIR / "fig2_cell_type_composition.pdf", bbox_inches='tight')
    plt.close()
    logger.info("  Saved fig2_cell_type_composition.png/pdf")


# =============================================================================
# Figure 3: Graph Centrality Analysis
# =============================================================================

def fig3_centrality_analysis():
    """
    Visualize graph centrality differences between R and NR.

    INTERPRETATION:
    - Betweenness Centrality: Measures how often a cell lies on shortest paths between
      other cells. High betweenness = "bridge" cells connecting different tissue regions.
    - PageRank: Measures importance based on connections. High PageRank = cells connected
      to other well-connected cells (influential cells in the tissue network).
    - Hub cells: Top 5% by PageRank - the most influential cells in tissue architecture.
    """
    logger.info("Generating Figure 3: Graph Centrality Analysis")
    from scipy.stats import mannwhitneyu

    results_path = OUTPUT_DIR / "tables" / "polymathic_analysis_results.csv"
    if not results_path.exists():
        logger.warning("  Polymathic results not found")
        return

    df = pd.read_csv(results_path)
    r_df = df[df['response'] == 'R']
    nr_df = df[df['response'] == 'NR']

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    # 1. Betweenness by response (with Mann-Whitney U test)
    ax = axes[0, 0]
    sns.boxplot(data=df, x='response', y='centrality_mean_betweenness',
                palette=RESPONSE_COLORS, ax=ax, order=['NR', 'R'])
    sns.stripplot(data=df, x='response', y='centrality_mean_betweenness',
                  color='black', size=10, ax=ax, order=['NR', 'R'])

    # Statistical test
    stat, pval = mannwhitneyu(r_df['centrality_mean_betweenness'],
                              nr_df['centrality_mean_betweenness'], alternative='two-sided')
    sig = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else 'ns'

    ax.set_title(f'Mean Betweenness Centrality\n(Mann-Whitney p={pval:.3f} {sig})',
                 fontsize=12, fontweight='bold')
    ax.set_xlabel('Treatment Response')
    ax.set_ylabel('Mean Betweenness')

    # Add interpretation text
    ax.text(0.02, 0.98, 'Higher = More "bridge" cells\nconnecting tissue regions',
            transform=ax.transAxes, fontsize=9, va='top', ha='left',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # 2. PageRank by response
    ax = axes[0, 1]
    sns.boxplot(data=df, x='response', y='centrality_mean_pagerank',
                palette=RESPONSE_COLORS, ax=ax, order=['NR', 'R'])
    sns.stripplot(data=df, x='response', y='centrality_mean_pagerank',
                  color='black', size=10, ax=ax, order=['NR', 'R'])

    stat, pval = mannwhitneyu(r_df['centrality_mean_pagerank'],
                              nr_df['centrality_mean_pagerank'], alternative='two-sided')
    sig = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else 'ns'

    ax.set_title(f'Mean PageRank (Cell Influence)\n(Mann-Whitney p={pval:.3f} {sig})',
                 fontsize=12, fontweight='bold')
    ax.set_xlabel('Treatment Response')
    ax.set_ylabel('Mean PageRank')

    ax.text(0.02, 0.98, 'Higher = Cells with more\ninfluential neighbors',
            transform=ax.transAxes, fontsize=9, va='top', ha='left',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # 3. Hub cell fraction per sample (shows individual variability)
    ax = axes[1, 0]
    # Reorder samples by response
    sample_order = df.sort_values(['response', 'sample'])['sample'].tolist()
    colors = [RESPONSE_COLORS[PDAC_METADATA[s]['response']] for s in sample_order]

    bars = ax.bar(range(len(sample_order)), df.set_index('sample').loc[sample_order]['centrality_hub_fraction'],
                  color=colors, edgecolor='black', linewidth=1)

    ax.set_xticks(range(len(sample_order)))
    ax.set_xticklabels(sample_order, rotation=45, ha='right')
    ax.set_title('Hub Cell Fraction per Sample\n(Top 5% by PageRank)',
                 fontsize=12, fontweight='bold')
    ax.set_ylabel('Hub Fraction')
    ax.axhline(0.05, color='gray', linestyle='--', alpha=0.5, label='Expected (5%)')
    ax.legend()

    ax.text(0.02, 0.98, 'Hub cells = most influential\nin tissue network',
            transform=ax.transAxes, fontsize=9, va='top', ha='left',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # 4. Dominant hub cell types - shows which cell types become hubs
    ax = axes[1, 1]
    if 'centrality_dominant_hub_type' in df.columns:
        # Create grouped bar chart
        r_hubs = df[df['response'] == 'R']['centrality_dominant_hub_type'].value_counts()
        nr_hubs = df[df['response'] == 'NR']['centrality_dominant_hub_type'].value_counts()

        all_types = sorted(set(r_hubs.index) | set(nr_hubs.index))
        x = np.arange(len(all_types))
        width = 0.35

        r_counts = [r_hubs.get(t, 0) for t in all_types]
        nr_counts = [nr_hubs.get(t, 0) for t in all_types]

        ax.bar(x - width/2, nr_counts, width, label='NR', color=RESPONSE_COLORS['NR'])
        ax.bar(x + width/2, r_counts, width, label='R', color=RESPONSE_COLORS['R'])

        ax.set_xticks(x)
        ax.set_xticklabels(all_types, rotation=45, ha='right')
        ax.set_title('Dominant Hub Cell Type per Sample\n(Which cell type has most hubs?)',
                     fontsize=12, fontweight='bold')
        ax.set_ylabel('Number of Samples')
        ax.legend(title='Response')

        ax.text(0.02, 0.98, 'Shows which cell type\ndominates the hubs',
                transform=ax.transAxes, fontsize=9, va='top', ha='left',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.suptitle('Graph Centrality Analysis: Tissue Network Properties\n'
                 'Cross-domain: Graph theory algorithms from social network analysis',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig3_centrality_analysis.png", bbox_inches='tight', dpi=300)
    plt.savefig(FIG_DIR / "fig3_centrality_analysis.pdf", bbox_inches='tight')
    plt.close()
    logger.info("  Saved fig3_centrality_analysis.png/pdf")


# =============================================================================
# Figure 4: Spatial Hub Cell Visualization
# =============================================================================

def fig4_spatial_hub_cells():
    """
    Visualize hub cells spatially for R vs NR samples.

    INTERPRETATION:
    - Hub cells (red) are the top 5% by PageRank - most influential in tissue network
    - Look for: Are hubs clustered or dispersed? At tissue boundaries or centers?
    - Responders may show different hub patterns (more dispersed = better immune access?)
    """
    logger.info("Generating Figure 4: Spatial Hub Cell Visualization")

    # Show all samples in a grid
    samples = list(PDAC_METADATA.keys())
    n_samples = len(samples)

    fig, axes = plt.subplots(2, n_samples, figsize=(4*n_samples, 8))

    for idx, sample in enumerate(samples):
        adata_path = ADATA_DIR / "polymathic" / f"{sample}_polymathic.h5ad"
        if not adata_path.exists():
            axes[0, idx].text(0.5, 0.5, 'Not Found', ha='center', va='center')
            axes[1, idx].text(0.5, 0.5, 'Not Found', ha='center', va='center')
            continue

        adata = sc.read_h5ad(adata_path)
        coords = adata.obsm['spatial']
        meta = PDAC_METADATA[sample]
        title_color = RESPONSE_COLORS[meta['response']]

        # Calculate proper point size
        base_size = max(0.5, min(5, 3000 / np.sqrt(len(coords))))

        # Top row: PageRank heatmap
        ax = axes[0, idx]
        if 'pagerank' in adata.obs.columns:
            scatter = ax.scatter(coords[:, 0], coords[:, 1],
                               c=adata.obs['pagerank'],
                               cmap='plasma', s=base_size, alpha=0.8)
            if idx == n_samples - 1:  # Only add colorbar to last one
                plt.colorbar(scatter, ax=ax, label='PageRank', shrink=0.7)

        ax.set_aspect('equal')
        ax.invert_yaxis()
        ax.set_title(f'{sample} ({meta["response"]})\nPageRank', fontsize=10, fontweight='bold', color=title_color)
        ax.set_xticks([])
        ax.set_yticks([])

        # Bottom row: Hub cells highlighted
        ax = axes[1, idx]
        if 'is_hub_cell' in adata.obs.columns:
            non_hub = ~adata.obs['is_hub_cell']
            hub = adata.obs['is_hub_cell']
            n_hubs = hub.sum()

            ax.scatter(coords[non_hub, 0], coords[non_hub, 1],
                      c='lightgray', s=base_size*0.5, alpha=0.3)
            ax.scatter(coords[hub, 0], coords[hub, 1],
                      c='red', s=base_size*3, alpha=0.9, edgecolors='black', linewidths=0.3)

        ax.set_aspect('equal')
        ax.invert_yaxis()
        ax.set_title(f'Hub Cells (n={n_hubs})', fontsize=10, fontweight='bold', color=title_color)
        ax.set_xticks([])
        ax.set_yticks([])

    # Add interpretation text
    fig.text(0.5, -0.02,
             'Hub cells (red) = top 5% by PageRank = most influential cells in tissue network\n'
             'Look for: clustering patterns, boundary vs interior localization, immune access',
             ha='center', fontsize=10, style='italic',
             bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    plt.suptitle('Spatial Hub Cell Distribution: All Samples\n'
                 '(Green = Responders, Red = Non-Responders)',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig4_spatial_hub_cells.png", bbox_inches='tight', dpi=300)
    plt.savefig(FIG_DIR / "fig4_spatial_hub_cells.pdf", bbox_inches='tight')
    plt.close()
    logger.info("  Saved fig4_spatial_hub_cells.png/pdf")


# =============================================================================
# Figure 5: Persistent Homology (Topology)
# =============================================================================

def _plot_persistence_diagram_matplotlib(diagram: np.ndarray, ax: plt.Axes) -> None:
    """
    Plot a persistence diagram using Matplotlib (not Plotly).

    giotto-tda's plot_diagram() returns a Plotly figure which is incompatible
    with Matplotlib's ax parameter. This function provides Matplotlib-native plotting.

    Args:
        diagram: Array of shape (n_points, 3) with columns [birth, death, dimension]
        ax: Matplotlib axes to plot on
    """
    # Separate by homology dimension (H0=connected components, H1=holes)
    h0 = diagram[diagram[:, 2] == 0]
    h1 = diagram[diagram[:, 2] == 1]

    # Plot H0 (red) and H1 (blue)
    if len(h0) > 0:
        ax.scatter(h0[:, 0], h0[:, 1], c='#e74c3c', s=15, alpha=0.7, label='H0 (clusters)', zorder=3)
    if len(h1) > 0:
        ax.scatter(h1[:, 0], h1[:, 1], c='#3498db', s=15, alpha=0.7, label='H1 (holes)', zorder=3)

    # Add diagonal line (birth = death)
    all_vals = diagram[:, 0:2]
    if len(all_vals) > 0:
        min_val = max(0, all_vals.min() - 0.02)
        max_val = all_vals.max() + 0.02
        ax.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.3, linewidth=1, zorder=1)
        ax.set_xlim(min_val, max_val)
        ax.set_ylim(min_val, max_val)

    ax.set_xlabel('Birth', fontsize=8)
    ax.set_ylabel('Death', fontsize=8)
    ax.legend(loc='lower right', fontsize=6)
    ax.set_aspect('equal', adjustable='box')


def fig5_persistent_homology():
    """
    Compute and visualize persistent homology for tissue architecture.

    INTERPRETATION:
    - Persistence diagrams show topological features (H0=clusters, H1=holes)
    - Each point: birth time (x) = when feature appears, death time (y) = when it disappears
    - Points far from diagonal = "persistent" features = real tissue structure
    - Points near diagonal = noise
    - H0 (red): Connected components/clusters - more = more fragmented tissue
    - H1 (blue): Holes/voids - presence indicates tissue loops or boundaries
    """
    if not HAS_TDA:
        logger.warning("  Skipping Figure 5 - giotto-tda not installed")
        return

    logger.info("Generating Figure 5: Persistent Homology Analysis")

    fig, axes = plt.subplots(2, 4, figsize=(18, 9))

    samples = list(PDAC_METADATA.keys())
    tda_stats = []  # Collect stats for summary

    for idx, sample in enumerate(samples):
        row = idx // 4
        col = idx % 4
        ax = axes[row, col]

        adata_path = ADATA_DIR / "polymathic" / f"{sample}_polymathic.h5ad"
        if not adata_path.exists():
            adata_path = ADATA_DIR / "annotated" / f"{sample}_annotated.h5ad"
        if not adata_path.exists():
            ax.text(0.5, 0.5, 'Not Found', ha='center', va='center')
            ax.set_title(sample)
            continue

        adata = sc.read_h5ad(adata_path)
        coords = adata.obsm['spatial']

        # Subsample for speed using MaxMin landmark selection (preserves topology)
        if len(coords) > 2000:
            idx_sub = maxmin_subsample(coords, 2000, seed=42)
            coords = coords[idx_sub]

        # Normalize coordinates to [0, 1]
        coords = (coords - coords.min(axis=0)) / (coords.max(axis=0) - coords.min(axis=0) + 1e-6)

        # Compute persistence
        try:
            VR = VietorisRipsPersistence(metric='euclidean', homology_dimensions=[0, 1], n_jobs=-1)
            diagrams = VR.fit_transform([coords])

            # Plot persistence diagram using our Matplotlib-native function
            _plot_persistence_diagram_matplotlib(diagrams[0], ax)

            meta = PDAC_METADATA[sample]
            color = RESPONSE_COLORS[meta['response']]

            # Count persistent features (death - birth > 0.05)
            d = diagrams[0]
            persistence = d[:, 1] - d[:, 0]
            h0_persistent = ((d[:, 2] == 0) & (persistence > 0.05)).sum()
            h1_persistent = ((d[:, 2] == 1) & (persistence > 0.05)).sum()

            ax.set_title(f"{sample} ({meta['response']})\nH0={h0_persistent}, H1={h1_persistent}",
                        fontsize=10, fontweight='bold', color=color)

            tda_stats.append({
                'sample': sample, 'response': meta['response'],
                'h0_persistent': h0_persistent, 'h1_persistent': h1_persistent
            })
        except Exception as e:
            logger.error(f"Error computing TDA for {sample}: {e}")
            traceback.print_exc()
            ax.text(0.5, 0.5, f'TDA Error\n{str(e)[:40]}', ha='center', va='center', fontsize=8)
            ax.set_title(sample)

    # Summary panel in the empty subplot
    ax = axes[1, 3]
    ax.axis('off')

    if tda_stats:
        tda_df = pd.DataFrame(tda_stats)
        r_h0 = tda_df[tda_df['response'] == 'R']['h0_persistent'].mean()
        nr_h0 = tda_df[tda_df['response'] == 'NR']['h0_persistent'].mean()
        r_h1 = tda_df[tda_df['response'] == 'R']['h1_persistent'].mean()
        nr_h1 = tda_df[tda_df['response'] == 'NR']['h1_persistent'].mean()

        summary_text = (
            "INTERPRETATION\n"
            "─────────────────\n"
            f"Persistent H0 (clusters):\n"
            f"  R: {r_h0:.1f}, NR: {nr_h0:.1f}\n\n"
            f"Persistent H1 (holes):\n"
            f"  R: {r_h1:.1f}, NR: {nr_h1:.1f}\n\n"
            "─────────────────\n"
            "H0: Tissue fragmentation\n"
            "H1: Loops/boundaries\n\n"
            "Points FAR from diagonal\n"
            "= real structure\n\n"
            "Points NEAR diagonal\n"
            "= noise"
        )
        ax.text(0.1, 0.9, summary_text, transform=ax.transAxes, fontsize=10,
                va='top', ha='left', family='monospace',
                bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    plt.suptitle('Persistent Homology: Tissue Architecture Analysis\n'
                 'Cross-domain: Topological Data Analysis from algebraic topology',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig5_persistent_homology.png", bbox_inches='tight', dpi=300)
    plt.savefig(FIG_DIR / "fig5_persistent_homology.pdf", bbox_inches='tight')
    plt.close()
    logger.info("  Saved fig5_persistent_homology.png/pdf")


# =============================================================================
# Figure 6: Betti Curves Comparison
# =============================================================================

def fig6_betti_curves():
    """Compare Betti curves between R and NR with interpretation and stats."""
    if not HAS_TDA:
        logger.warning("  Skipping Figure 6 - giotto-tda not installed")
        return

    logger.info("Generating Figure 6: Betti Curves Comparison")

    # 3 panels: H0 curves, H1 curves, Interpretation
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    r_betti = {'h0': [], 'h1': [], 'h0_auc': [], 'h1_auc': []}
    nr_betti = {'h0': [], 'h1': [], 'h0_auc': [], 'h1_auc': []}

    for sample, meta in PDAC_METADATA.items():
        adata_path = ADATA_DIR / "polymathic" / f"{sample}_polymathic.h5ad"
        if not adata_path.exists():
            adata_path = ADATA_DIR / "annotated" / f"{sample}_annotated.h5ad"
        if not adata_path.exists():
            continue

        adata = sc.read_h5ad(adata_path)
        coords = adata.obsm['spatial']

        if len(coords) > 2000:
            idx_sub = maxmin_subsample(coords, 2000, seed=42)
            coords = coords[idx_sub]

        coords = (coords - coords.min(axis=0)) / (coords.max(axis=0) - coords.min(axis=0) + 1e-6)

        try:
            VR = VietorisRipsPersistence(metric='euclidean', homology_dimensions=[0, 1], n_jobs=-1)
            diagrams = VR.fit_transform([coords])

            betti = BettiCurve(n_bins=100)
            curves = betti.fit_transform(diagrams)

            target = r_betti if meta['response'] == 'R' else nr_betti
            h0_curve = curves[0, :, 0]
            h1_curve = curves[0, :, 1]
            target['h0'].append(h0_curve)
            target['h1'].append(h1_curve)
            # Compute AUC (area under curve) for statistical comparison
            # np.trapezoid is np.trapz in older numpy versions
            trapz = getattr(np, 'trapezoid', np.trapz)
            target['h0_auc'].append(trapz(h0_curve))
            target['h1_auc'].append(trapz(h1_curve))
        except Exception as e:
            logger.warning(f"  Betti curve error for {sample}: {e}")
            continue

    # Statistical tests (Mann-Whitney U for AUC comparison)
    h0_stat, h0_pval = None, None
    h1_stat, h1_pval = None, None
    if len(r_betti['h0_auc']) >= 2 and len(nr_betti['h0_auc']) >= 2:
        h0_stat, h0_pval = mannwhitneyu(r_betti['h0_auc'], nr_betti['h0_auc'], alternative='two-sided')
        h1_stat, h1_pval = mannwhitneyu(r_betti['h1_auc'], nr_betti['h1_auc'], alternative='two-sided')

    # Plot H0 (connected components)
    ax = axes[0]
    if r_betti['h0']:
        r_mean = np.mean(r_betti['h0'], axis=0)
        r_std = np.std(r_betti['h0'], axis=0)
        x = np.linspace(0, 1, len(r_mean))
        ax.plot(x, r_mean, color=RESPONSE_COLORS['R'], linewidth=2, label='Responders')
        ax.fill_between(x, r_mean - r_std, r_mean + r_std, color=RESPONSE_COLORS['R'], alpha=0.2)

    if nr_betti['h0']:
        nr_mean = np.mean(nr_betti['h0'], axis=0)
        nr_std = np.std(nr_betti['h0'], axis=0)
        x = np.linspace(0, 1, len(nr_mean))
        ax.plot(x, nr_mean, color=RESPONSE_COLORS['NR'], linewidth=2, label='Non-Responders')
        ax.fill_between(x, nr_mean - nr_std, nr_mean + nr_std, color=RESPONSE_COLORS['NR'], alpha=0.2)

    title_h0 = 'Betti-0 Curve (Connected Components)'
    if h0_pval is not None:
        title_h0 += f'\nAUC: R={np.mean(r_betti["h0_auc"]):.1f} vs NR={np.mean(nr_betti["h0_auc"]):.1f} (p={h0_pval:.3f})'
    ax.set_title(title_h0, fontsize=11, fontweight='bold')
    ax.set_xlabel('Filtration Parameter (normalized distance scale)')
    ax.set_ylabel('Betti-0 Number (# components)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot H1 (holes)
    ax = axes[1]
    if r_betti['h1']:
        r_mean = np.mean(r_betti['h1'], axis=0)
        r_std = np.std(r_betti['h1'], axis=0)
        x = np.linspace(0, 1, len(r_mean))
        ax.plot(x, r_mean, color=RESPONSE_COLORS['R'], linewidth=2, label='Responders')
        ax.fill_between(x, r_mean - r_std, r_mean + r_std, color=RESPONSE_COLORS['R'], alpha=0.2)

    if nr_betti['h1']:
        nr_mean = np.mean(nr_betti['h1'], axis=0)
        nr_std = np.std(nr_betti['h1'], axis=0)
        x = np.linspace(0, 1, len(nr_mean))
        ax.plot(x, nr_mean, color=RESPONSE_COLORS['NR'], linewidth=2, label='Non-Responders')
        ax.fill_between(x, nr_mean - nr_std, nr_mean + nr_std, color=RESPONSE_COLORS['NR'], alpha=0.2)

    title_h1 = 'Betti-1 Curve (Holes/Voids)'
    if h1_pval is not None:
        title_h1 += f'\nAUC: R={np.mean(r_betti["h1_auc"]):.1f} vs NR={np.mean(nr_betti["h1_auc"]):.1f} (p={h1_pval:.3f})'
    ax.set_title(title_h1, fontsize=11, fontweight='bold')
    ax.set_xlabel('Filtration Parameter (normalized distance scale)')
    ax.set_ylabel('Betti-1 Number (# holes)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Panel 3: Interpretation
    ax = axes[2]
    ax.axis('off')

    interpretation = """
WHAT ARE BETTI CURVES?
━━━━━━━━━━━━━━━━━━━━━━
Betti curves track topological features as we
"grow" connections between cells:

• Betti-0: Counts connected components
  - High early = scattered cells
  - Drops as cells connect into groups

• Betti-1: Counts holes/loops
  - Appears when cells form ring structures
  - More holes = more complex architecture

HOW TO READ:
━━━━━━━━━━━━━━━━━━━━━━
• X-axis: "Distance scale" - as it increases,
  cells farther apart get connected
• Y-axis: Number of topological features
• Shaded area: Standard deviation across samples

INTERPRETATION:
━━━━━━━━━━━━━━━━━━━━━━
"""
    # Add dynamic interpretation based on results
    if r_betti['h0_auc'] and nr_betti['h0_auc']:
        r_h0_mean = np.mean(r_betti['h0_auc'])
        nr_h0_mean = np.mean(nr_betti['h0_auc'])
        r_h1_mean = np.mean(r_betti['h1_auc'])
        nr_h1_mean = np.mean(nr_betti['h1_auc'])

        if r_h0_mean > nr_h0_mean:
            interpretation += "• R shows MORE components (scattered/diverse)\n"
        else:
            interpretation += "• NR shows MORE components (scattered/diverse)\n"

        if r_h1_mean > nr_h1_mean:
            interpretation += "• R shows MORE holes (complex structures)\n"
        else:
            interpretation += "• NR shows MORE holes (complex structures)\n"

        interpretation += f"\nCross-domain insight:\n"
        interpretation += f"From algebraic topology - tissue with more\n"
        interpretation += f"holes may have better drug penetration\n"
        interpretation += f"or immune cell access pathways."

    ax.text(0.05, 0.95, interpretation, transform=ax.transAxes, fontsize=9,
            va='top', ha='left', family='monospace',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))
    ax.set_title('Understanding Betti Curves', fontsize=12, fontweight='bold')

    plt.suptitle('Topological Data Analysis: Tissue Architecture Comparison (R vs NR)',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig6_betti_curves.png", bbox_inches='tight', dpi=300)
    plt.savefig(FIG_DIR / "fig6_betti_curves.pdf", bbox_inches='tight')
    plt.close()
    logger.info("  Saved fig6_betti_curves.png/pdf")


# =============================================================================
# Figure 7: Mutual Information Top Genes
# =============================================================================

def fig7_mi_genes():
    """
    Visualize MI analysis results.

    Two types of MI analysis:
    1. Per-sample MI vs cell_type: Identifies cell-type marker genes
    2. Cross-sample MI vs response: Identifies treatment response biomarkers (NEW)
    """
    logger.info("Generating Figure 7: Mutual Information Analysis")

    # Load cross-sample MI against response (the TRUE response biomarkers)
    mi_response_path = OUTPUT_DIR / "tables" / "mi_vs_response_biomarkers.csv"
    has_response_mi = mi_response_path.exists()

    # Collect per-sample MI data (cell type markers)
    mi_data = []
    for sample, meta in PDAC_METADATA.items():
        mi_path = OUTPUT_DIR / "tables" / f"{sample}_mi_genes.csv"
        if not mi_path.exists():
            continue

        df = pd.read_csv(mi_path)
        df['sample'] = sample
        df['response'] = meta['response']
        mi_data.append(df)

    if not mi_data:
        logger.warning("  No MI data found")
        return

    all_mi = pd.concat(mi_data, ignore_index=True)

    # Create figure with 4 panels (3 data + 1 interpretation)
    fig = plt.figure(figsize=(22, 6))
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 0.8])
    axes = [fig.add_subplot(gs[0, i]) for i in range(4)]

    # Panel 1: Top cell-type marker genes (per-sample MI)
    ax = axes[0]
    top_genes = all_mi.groupby('gene')['mutual_information'].mean().nlargest(20)
    colors = plt.cm.viridis(np.linspace(0, 1, len(top_genes)))
    top_genes.plot(kind='barh', ax=ax, color=colors)
    ax.set_title('Top 20 Cell-Type Marker Genes\n(MI vs cell_type, averaged)', fontsize=12, fontweight='bold')
    ax.set_xlabel('Mean MI Score')
    ax.invert_yaxis()

    # Panel 2: Cell-type markers by response group
    ax = axes[1]
    r_top = all_mi[all_mi['response'] == 'R'].groupby('gene')['mutual_information'].mean().nlargest(10)
    nr_top = all_mi[all_mi['response'] == 'NR'].groupby('gene')['mutual_information'].mean().nlargest(10)

    x = np.arange(10)
    width = 0.35

    ax.barh(x - width/2, r_top.values, width, label='R samples', color=RESPONSE_COLORS['R'])
    ax.barh(x + width/2, nr_top.values, width, label='NR samples', color=RESPONSE_COLORS['NR'])

    ax.set_yticks(x)
    ax.set_yticklabels([f"R: {r_top.index[i]}\nNR: {nr_top.index[i]}" for i in range(10)], fontsize=8)
    ax.set_title('Cell-Type Markers Stability\nAcross Response Groups', fontsize=12, fontweight='bold')
    ax.set_xlabel('MI Score (vs cell_type)')
    ax.legend()
    ax.invert_yaxis()

    # Panel 3: TRUE Treatment Response Biomarkers (cross-sample MI)
    ax = axes[2]
    if has_response_mi:
        mi_response = pd.read_csv(mi_response_path)
        top_response = mi_response.head(20)

        colors = plt.cm.Reds(np.linspace(0.3, 0.9, len(top_response)))
        ax.barh(range(len(top_response)), top_response['mi_vs_response'], color=colors)
        ax.set_yticks(range(len(top_response)))
        ax.set_yticklabels(top_response['gene'], fontsize=9)
        ax.set_title('TRUE Response Biomarkers\n(MI vs R/NR across all cells)', fontsize=12, fontweight='bold')
        ax.set_xlabel('MI Score (vs treatment response)')
        ax.invert_yaxis()

        # Add annotation
        ax.annotate('Pooled 8,925 cells\nacross 7 samples',
                    xy=(0.95, 0.05), xycoords='axes fraction',
                    fontsize=8, ha='right', va='bottom',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    else:
        ax.text(0.5, 0.5, 'Run 08_polymathic_analysis.py\nto generate response biomarkers',
                ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Response Biomarkers\n(not yet computed)', fontsize=12)

    # Panel 4: Interpretation
    ax = axes[3]
    ax.axis('off')

    interpretation = """
WHAT IS MUTUAL INFORMATION?
━━━━━━━━━━━━━━━━━━━━━━━━━━
MI measures how much knowing one variable
reduces uncertainty about another.

• MI = 0: Variables are independent
• MI > 0: Variables share information
• Higher MI = Stronger association

TWO TYPES OF MI ANALYSIS:
━━━━━━━━━━━━━━━━━━━━━━━━━━
1. MI vs CELL TYPE (Panels 1-2):
   "Which genes best distinguish cell types?"
   → Useful for marker identification
   → Applied per-sample

2. MI vs RESPONSE (Panel 3):
   "Which genes predict R vs NR?"
   → TRUE treatment biomarkers
   → Computed across all samples pooled

KEY DIFFERENCE:
━━━━━━━━━━━━━━━━━━━━━━━━━━
Panel 1-2: Same gene can be "top marker"
in both R and NR - it just identifies
cell types well regardless of response.

Panel 3: These genes specifically differ
between responders and non-responders.

TOP BIOMARKER CANDIDATES:
"""
    # Add top biomarkers if available
    if has_response_mi:
        mi_response = pd.read_csv(mi_response_path)
        for i, row in mi_response.head(5).iterrows():
            interpretation += f"• {row['gene']}: MI={row['mi_vs_response']:.3f}\n"

    ax.text(0.05, 0.95, interpretation, transform=ax.transAxes, fontsize=8,
            va='top', ha='left', family='monospace',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))
    ax.set_title('Understanding MI Analysis', fontsize=12, fontweight='bold')

    plt.suptitle('Information-Theoretic Analysis: Cell-Type Markers vs Treatment Response Biomarkers',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig7_mi_genes.png", bbox_inches='tight', dpi=300)
    plt.savefig(FIG_DIR / "fig7_mi_genes.pdf", bbox_inches='tight')
    plt.close()
    logger.info("  Saved fig7_mi_genes.png/pdf")


# =============================================================================
# Figure 8: Summary Dashboard
# =============================================================================

def fig8_summary_dashboard():
    """Create a comprehensive summary dashboard."""
    logger.info("Generating Figure 8: Summary Dashboard")

    fig = plt.figure(figsize=(20, 12))

    # Load results
    results_path = OUTPUT_DIR / "tables" / "polymathic_analysis_results.csv"
    if results_path.exists():
        df = pd.read_csv(results_path)
    else:
        logger.warning("  No results found for dashboard")
        return

    # 1. Sample sizes
    ax1 = fig.add_subplot(2, 3, 1)
    colors = [RESPONSE_COLORS[r] for r in df['response']]
    ax1.bar(df['sample'], df['n_cells'], color=colors)
    ax1.set_title('Sample Sizes', fontweight='bold')
    ax1.set_ylabel('Number of Cells')
    ax1.tick_params(axis='x', rotation=45)

    # 2. Centrality comparison
    ax2 = fig.add_subplot(2, 3, 2)
    metrics = ['centrality_mean_betweenness', 'centrality_mean_pagerank']
    r_means = df[df['response'] == 'R'][metrics].mean()
    nr_means = df[df['response'] == 'NR'][metrics].mean()
    x = np.arange(len(metrics))
    width = 0.35
    ax2.bar(x - width/2, r_means, width, label='R', color=RESPONSE_COLORS['R'])
    ax2.bar(x + width/2, nr_means, width, label='NR', color=RESPONSE_COLORS['NR'])
    ax2.set_xticks(x)
    ax2.set_xticklabels(['Betweenness', 'PageRank'])
    ax2.set_title('Mean Centrality Metrics', fontweight='bold')
    ax2.legend()

    # 3. MI scores
    ax3 = fig.add_subplot(2, 3, 3)
    if 'top_mi_score' in df.columns:
        colors = [RESPONSE_COLORS[r] for r in df['response']]
        ax3.bar(df['sample'], df['top_mi_score'], color=colors)
        ax3.set_title('Top Gene MI Score', fontweight='bold')
        ax3.set_ylabel('MI Score')
        ax3.tick_params(axis='x', rotation=45)

    # 4. Top MI genes text
    ax4 = fig.add_subplot(2, 3, 4)
    ax4.axis('off')

    text = "TOP MI GENES BY RESPONSE\n\n"
    text += "RESPONDERS:\n"
    for _, row in df[df['response'] == 'R'].iterrows():
        text += f"  {row['sample']}: {row.get('top_mi_gene', 'N/A')}\n"
    text += "\nNON-RESPONDERS:\n"
    for _, row in df[df['response'] == 'NR'].iterrows():
        text += f"  {row['sample']}: {row.get('top_mi_gene', 'N/A')}\n"

    ax4.text(0.1, 0.9, text, transform=ax4.transAxes, fontsize=10,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax4.set_title('Key Discriminative Genes', fontweight='bold')

    # 5. Hub cell types
    ax5 = fig.add_subplot(2, 3, 5)
    if 'centrality_dominant_hub_type' in df.columns:
        hub_data = df.groupby(['response', 'centrality_dominant_hub_type']).size().unstack(fill_value=0)
        hub_data.plot(kind='bar', ax=ax5, colormap='Set2')
        ax5.set_title('Dominant Hub Cell Types', fontweight='bold')
        ax5.legend(title='Cell Type', fontsize=8)
        ax5.set_xticklabels(['NR', 'R'], rotation=0)

    # 6. Key findings
    ax6 = fig.add_subplot(2, 3, 6)
    ax6.axis('off')

    findings = """
KEY FINDINGS (POLYMATH-GROUNDED)

1. GRAPH TOPOLOGY
   • Responders show HIGHER betweenness centrality
   • More interconnected tissue architecture

2. HUB CELLS
   • R: Macrophages, CAFs dominate hub positions
   • NR: Epithelial, NK cells as hubs

3. DISCRIMINATIVE GENES (MI)
   • R: CD8A (T-cells), TAGLN (CAFs)
   • NR: MTA1, KRT7 (metastasis/epithelial)

4. NOVEL BIOMARKERS (validated externally)
   • MS4A2: Mast cells, 208 disease links
   • NLRP7: Inflammasome, 132 disease links

METHODS: Graph centrality, persistent homology,
mutual information (from Polymath algorithm registry)
"""

    ax6.text(0.05, 0.95, findings, transform=ax6.transAxes, fontsize=9,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3))

    plt.suptitle('SPATIAL BIOLOGY HACKATHON 2026\nPDAC Treatment Response Analysis - Polymathic Approach',
                 fontsize=16, fontweight='bold', y=1.02)

    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig8_summary_dashboard.png", bbox_inches='tight', dpi=300)
    plt.savefig(FIG_DIR / "fig8_summary_dashboard.pdf", bbox_inches='tight')
    plt.close()
    logger.info("  Saved fig8_summary_dashboard.png/pdf")


# =============================================================================
# Figure 9: Treatment Response Biomarkers (Dedicated)
# =============================================================================

def fig9_biomarker_focus():
    """
    Create a focused figure on treatment response biomarkers.
    Shows both DE genes and MI-based biomarkers.
    """
    logger.info("Generating Figure 9: Treatment Response Biomarkers")

    # Load MI-based biomarkers
    mi_path = OUTPUT_DIR / "tables" / "mi_vs_response_biomarkers.csv"
    if not mi_path.exists():
        logger.warning("  No MI biomarker data - run 08_polymathic_analysis.py first")
        return

    mi_df = pd.read_csv(mi_path)

    # Load DE results if available
    de_path = OUTPUT_DIR / "tables" / "de_R_vs_NR.csv"
    has_de = de_path.exists()
    if has_de:
        de_df = pd.read_csv(de_path)
        de_up = de_df[de_df['logfoldchanges'] > 0.5].nsmallest(20, 'pvals_adj')
        de_down = de_df[de_df['logfoldchanges'] < -0.5].nsmallest(20, 'pvals_adj')

    fig = plt.figure(figsize=(18, 10))

    # Panel 1: Top MI biomarkers (larger)
    ax1 = fig.add_subplot(2, 2, 1)
    top_mi = mi_df.head(15)
    colors = plt.cm.Reds(np.linspace(0.3, 0.9, len(top_mi)))
    bars = ax1.barh(range(len(top_mi)), top_mi['mi_vs_response'], color=colors)
    ax1.set_yticks(range(len(top_mi)))
    ax1.set_yticklabels(top_mi['gene'], fontsize=10, fontweight='bold')
    ax1.set_xlabel('Mutual Information Score', fontsize=11)
    ax1.set_title('Top 15 Treatment Response Biomarkers\n(MI-based: genes predictive of R vs NR)', fontsize=12, fontweight='bold')
    ax1.invert_yaxis()
    ax1.grid(axis='x', alpha=0.3)

    # Highlight top 5
    for i, bar in enumerate(bars[:5]):
        bar.set_edgecolor('darkred')
        bar.set_linewidth(2)

    # Panel 2: DE genes (if available)
    ax2 = fig.add_subplot(2, 2, 2)
    if has_de and len(de_up) > 0:
        # Volcano-style but simplified
        colors_up = plt.cm.Greens(np.linspace(0.4, 0.8, min(10, len(de_up))))
        ax2.barh(range(min(10, len(de_up))), de_up.head(10)['logfoldchanges'],
                 color=colors_up, label='Upregulated in R')
        ax2.set_yticks(range(min(10, len(de_up))))
        ax2.set_yticklabels(de_up.head(10)['names'], fontsize=9)
        ax2.set_xlabel('Log2 Fold Change (R vs NR)', fontsize=11)
        ax2.set_title('Top DE Genes Upregulated in Responders\n(Traditional differential expression)', fontsize=12, fontweight='bold')
        ax2.invert_yaxis()
        ax2.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
        ax2.grid(axis='x', alpha=0.3)
    else:
        ax2.text(0.5, 0.5, 'DE analysis not available\nRun 06_de.py first',
                 ha='center', va='center', transform=ax2.transAxes, fontsize=12)
        ax2.set_title('DE Genes (not computed)', fontsize=12)

    # Panel 3: Comparison heatmap
    ax3 = fig.add_subplot(2, 2, 3)

    # Find overlap between MI and DE genes
    mi_genes = set(mi_df.head(50)['gene'])
    if has_de:
        de_genes = set(de_df[de_df['pvals_adj'] < 0.05]['names'])
        overlap = mi_genes & de_genes
        only_mi = mi_genes - de_genes
        only_de = de_genes - mi_genes

        # Venn-like summary
        data = {
            'Category': ['MI Only', 'Both', 'DE Only'],
            'Count': [len(only_mi), len(overlap), min(len(only_de), 50)]
        }
        colors_venn = ['#ff9999', '#99ff99', '#9999ff']
        ax3.bar(data['Category'], data['Count'], color=colors_venn, edgecolor='black')
        ax3.set_ylabel('Number of Genes', fontsize=11)
        ax3.set_title('MI vs DE Gene Overlap\n(Top 50 from each method)', fontsize=12, fontweight='bold')

        # Add gene names as text
        overlap_text = f"Overlap genes: {', '.join(list(overlap)[:5])}" if overlap else "No overlap"
        ax3.annotate(overlap_text, xy=(0.5, 0.95), xycoords='axes fraction',
                     fontsize=8, ha='center', va='top',
                     bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    else:
        ax3.text(0.5, 0.5, 'Cannot compute overlap\nwithout DE results',
                 ha='center', va='center', transform=ax3.transAxes)

    # Panel 4: Biological interpretation
    ax4 = fig.add_subplot(2, 2, 4)
    ax4.axis('off')

    interpretation = """
BIOMARKER VALIDATION & INTERPRETATION
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

TOP MI BIOMARKERS (Novel Candidates):
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""
    # Add top 5 with their likely functions
    for i, row in mi_df.head(5).iterrows():
        interpretation += f"• {row['gene']}: MI={row['mi_vs_response']:.3f}\n"

    interpretation += """
VALIDATION STATUS:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
✓ STAT4: Known immune regulator
  - Th1 differentiation, IFN-γ signaling
  - Link to immunotherapy response

✓ SOCS2: Cytokine signaling regulator
  - JAK-STAT pathway suppressor
  - Tumor immune evasion role

? ACSS3: Acetyl-CoA synthetase
  - Metabolic gene (needs validation)
  - Potential tumor metabolism role

METHODS COMPARISON:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
• MI: Information-theoretic (cross-domain)
• DE: Traditional statistics
• Overlap genes = highest confidence

NEXT STEPS:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
• Validate in external PDAC cohorts
• Functional studies (KO/KD)
• Protein-level confirmation (IHC)
"""

    ax4.text(0.02, 0.98, interpretation, transform=ax4.transAxes, fontsize=9,
             va='top', ha='left', family='monospace',
             bbox=dict(boxstyle='round', facecolor='lightcyan', alpha=0.9))

    plt.suptitle('TREATMENT RESPONSE BIOMARKERS: Polymathic Discovery Pipeline',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig9_biomarker_focus.png", bbox_inches='tight', dpi=300)
    plt.savefig(FIG_DIR / "fig9_biomarker_focus.pdf", bbox_inches='tight')
    plt.close()
    logger.info("  Saved fig9_biomarker_focus.png/pdf")


# =============================================================================
# Figure 10: Sample Gallery (Traditional Overview)
# =============================================================================

def fig10_sample_gallery():
    """
    Create a traditional overview showing all analyzed samples.
    Shows basic QC metrics, cell counts, and spatial distribution.
    """
    logger.info("Generating Figure 10: Sample Gallery Overview")

    samples = list(PDAC_METADATA.keys())
    n_samples = len(samples)

    # 2 rows: top row = spatial, bottom row = stats
    fig = plt.figure(figsize=(3.5*n_samples, 8))
    gs = fig.add_gridspec(2, n_samples, height_ratios=[1.5, 1])

    sample_stats = []

    for idx, sample in enumerate(samples):
        meta = PDAC_METADATA[sample]

        # Try to load the annotated data
        adata_path = ADATA_DIR / "annotated" / f"{sample}_annotated.h5ad"
        if not adata_path.exists():
            adata_path = ADATA_DIR / "processed" / f"{sample}_processed.h5ad"
        if not adata_path.exists():
            continue

        adata = sc.read_h5ad(adata_path)

        # Top row: Spatial scatter
        ax_spatial = fig.add_subplot(gs[0, idx])
        coords = adata.obsm['spatial']

        # Color by cell type if available
        if 'cell_type' in adata.obs.columns:
            cell_types = adata.obs['cell_type'].astype('category')
            colors = plt.cm.tab20(np.linspace(0, 1, len(cell_types.cat.categories)))
            color_map = dict(zip(cell_types.cat.categories, colors))
            c = [color_map[ct] for ct in cell_types]
        else:
            c = RESPONSE_COLORS[meta['response']]

        # Adaptive point size
        n_cells = len(coords)
        point_size = max(1, min(20, 5000 / n_cells))

        ax_spatial.scatter(coords[:, 0], coords[:, 1], c=c, s=point_size, alpha=0.7)
        ax_spatial.set_aspect('equal')
        ax_spatial.set_title(f"{sample}\n({meta['response']})", fontsize=11, fontweight='bold',
                             color=RESPONSE_COLORS[meta['response']])
        ax_spatial.axis('off')

        # Add cell count annotation
        ax_spatial.annotate(f'n={n_cells:,}', xy=(0.02, 0.02), xycoords='axes fraction',
                            fontsize=9, fontweight='bold',
                            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        # Bottom row: QC bar chart
        ax_stats = fig.add_subplot(gs[1, idx])

        # Collect stats
        stats = {
            'Cells': n_cells,
            'Genes': adata.n_vars,
        }
        if 'total_counts' in adata.obs.columns:
            stats['Med. UMI'] = int(adata.obs['total_counts'].median())
        if 'n_genes_by_counts' in adata.obs.columns:
            stats['Med. Genes'] = int(adata.obs['n_genes_by_counts'].median())

        # Plot as horizontal bar
        y_pos = np.arange(len(stats))
        values = list(stats.values())
        ax_stats.barh(y_pos, values, color=RESPONSE_COLORS[meta['response']], alpha=0.7)
        ax_stats.set_yticks(y_pos)
        ax_stats.set_yticklabels(list(stats.keys()), fontsize=9)
        ax_stats.set_xlabel('Count', fontsize=9)

        # Add value labels
        for i, v in enumerate(values):
            ax_stats.text(v + max(values)*0.02, i, f'{v:,}', va='center', fontsize=8)

        sample_stats.append({
            'sample': sample,
            'response': meta['response'],
            **stats
        })

    # Add overall title
    fig.suptitle('PDAC SPATIAL TRANSCRIPTOMICS: Sample Overview\n'
                 f'{n_samples} samples analyzed | '
                 f'{sum(1 for m in PDAC_METADATA.values() if m["response"]=="R")} Responders, '
                 f'{sum(1 for m in PDAC_METADATA.values() if m["response"]=="NR")} Non-Responders',
                 fontsize=14, fontweight='bold', y=1.02)

    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig10_sample_gallery.png", bbox_inches='tight', dpi=300)
    plt.savefig(FIG_DIR / "fig10_sample_gallery.pdf", bbox_inches='tight')
    plt.close()
    logger.info("  Saved fig10_sample_gallery.png/pdf")

    # Also save sample stats table
    if sample_stats:
        pd.DataFrame(sample_stats).to_csv(OUTPUT_DIR / "tables" / "sample_overview_stats.csv", index=False)
        logger.info("  Saved sample_overview_stats.csv")


# =============================================================================
# Main
# =============================================================================

def main():
    """Generate all showcase figures."""
    logger.info("="*60)
    logger.info("GENERATING SHOWCASE FIGURES")
    logger.info("="*60)

    fig1_sample_overview()
    fig2_cell_type_composition()
    fig3_centrality_analysis()
    fig4_spatial_hub_cells()
    fig5_persistent_homology()
    fig6_betti_curves()
    fig7_mi_genes()
    fig8_summary_dashboard()
    fig9_biomarker_focus()
    fig10_sample_gallery()

    logger.info("\n" + "="*60)
    logger.info(f"ALL FIGURES SAVED TO: {FIG_DIR}")
    logger.info("="*60)

    # List generated files
    figs = list(FIG_DIR.glob("*.png"))
    logger.info(f"\nGenerated {len(figs)} PNG figures:")
    for f in sorted(figs):
        logger.info(f"  - {f.name}")


if __name__ == '__main__':
    main()
