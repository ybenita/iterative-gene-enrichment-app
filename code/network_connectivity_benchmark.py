#!/usr/bin/env python3
"""
Network Connectivity Benchmarking for iGEA

This script computes network connectivity metrics from permutation results
to establish a null distribution for benchmarking real gene lists.

The hypothesis: Real gene lists will have better network connectivity than
random gene lists, as measured by various graph metrics.

Connectivity metrics computed:
1. Average connections per gene
2. Network density
3. Number of connected components
4. Gene centrality distribution
5. Term connectivity
6. Network clustering coefficient
"""

import sys
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Set, Optional
import pandas as pd
import numpy as np
from collections import defaultdict
import json

# Add code directory to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "code"))

from iter_enrichment import IterativeEnrichment

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


class NetworkConnectivityAnalyzer:
    """
    Analyzes network connectivity from iGEA results.
    
    Builds a bipartite graph where:
    - Nodes: genes and terms
    - Edges: gene-term connections (gene is in term)
    """
    
    def __init__(self):
        self.genes: Set[str] = set()
        self.terms: Set[str] = set()
        self.edges: List[Tuple[str, str]] = []  # (gene, term) pairs
        self.gene_to_terms: Dict[str, Set[str]] = defaultdict(set)
        self.term_to_genes: Dict[str, Set[str]] = defaultdict(set)
        self.term_to_library: Dict[str, str] = {}  # Map term to library name
    
    def add_igea_results(self, results: List[Dict]) -> None:
        """
        Add iGEA results to build the network.
        
        Args:
            results: List of iGEA result dictionaries
        """
        for record in results:
            term = record.get("Term", "")
            if not term:
                continue
            
            # Get library name (if available)
            library = record.get("Library", "")
            
            # Get genes removed for this iteration (genes in the term)
            genes_removed = record.get("Genes removed for next iteration", [])
            if isinstance(genes_removed, str):
                genes_removed = [g.strip() for g in genes_removed.split(",") if g.strip()]
            elif not isinstance(genes_removed, list):
                continue
            
            # Add term
            self.terms.add(term)
            
            # Track library for this term
            if library:
                self.term_to_library[term] = library
            
            # Add genes and edges
            for gene in genes_removed:
                if gene:
                    self.genes.add(gene)
                    self.edges.append((gene, term))
                    self.gene_to_terms[gene].add(term)
                    self.term_to_genes[term].add(gene)
    
    def compute_metrics(self, include_library_diversity: bool = True) -> Dict:
        """
        Compute network connectivity metrics.
        
        Args:
            include_library_diversity: If False, skip library diversity metrics
                                     (use for individual library analysis where diversity is always 0/1)
        
        Returns:
            Dictionary of connectivity metrics
        """
        n_genes = len(self.genes)
        n_terms = len(self.terms)
        n_edges = len(self.edges)
        
        if n_genes == 0 or n_terms == 0:
            base_metrics = {
                "n_genes": 0,
                "n_terms": 0,
                "n_edges": 0,
                "avg_connections_per_gene": 0.0,
                "avg_connections_per_term": 0.0,
                "n_connected_components": 0,
                "largest_component_size": 0,
                "gene_centrality_mean": 0.0,
                "gene_centrality_std": 0.0,
                "gene_centrality_max": 0,
                "hub_genes_count": 0,
                "term_centrality_mean": 0.0,
                "term_centrality_std": 0.0,
                "term_centrality_max": 0,
                "largest_cluster_genes": 0,
                "largest_cluster_terms": 0,
                "largest_cluster_edges": 0,
                "largest_cluster_avg_edges_per_gene": 0.0,
                "avg_cluster_size": 0.0,
                "median_cluster_size": 0.0,
                "cluster_size_std": 0.0,
                "avg_cluster_avg_edges_per_gene": 0.0,
                "weighted_avg_cluster_avg_edges_per_gene": 0.0,
                "fraction_in_largest_cluster": 0.0,
                "fraction_edges_in_largest_cluster": 0.0,
                "clustering_coefficient": 0.0,
            }
            if include_library_diversity:
                base_metrics["largest_cluster_libraries"] = 0
                base_metrics["avg_cluster_libraries"] = 0.0
            return base_metrics
        
        # Basic metrics
        avg_connections_per_gene = n_edges / n_genes if n_genes > 0 else 0.0
        avg_connections_per_term = n_edges / n_terms if n_terms > 0 else 0.0
        
        # Gene centrality (degree)
        gene_centralities = [len(self.gene_to_terms[gene]) for gene in self.genes]
        gene_centrality_mean = np.mean(gene_centralities) if gene_centralities else 0.0
        gene_centrality_std = np.std(gene_centralities) if gene_centralities else 0.0
        gene_centrality_max = max(gene_centralities) if gene_centralities else 0
        
        # Term centrality (degree)
        term_centralities = [len(self.term_to_genes[term]) for term in self.terms]
        term_centrality_mean = np.mean(term_centralities) if term_centralities else 0.0
        term_centrality_std = np.std(term_centralities) if term_centralities else 0.0
        term_centrality_max = max(term_centralities) if term_centralities else 0
        
        # Connected components (clusters) - returns full cluster information
        n_components, largest_component_size, clusters = self._compute_connected_components()
        
        # Cluster-based metrics
        if clusters:
            cluster_sizes = [c['size'] for c in clusters]
            cluster_densities = [c.get('avg_edges_per_gene', c.get('density', 0)) for c in clusters]
            cluster_libraries = [c['n_libraries'] for c in clusters]
            
            # Find largest cluster
            largest_cluster = max(clusters, key=lambda x: x['size'])
            
            # Cluster distribution metrics
            avg_cluster_size = np.mean(cluster_sizes) if cluster_sizes else 0.0
            median_cluster_size = np.median(cluster_sizes) if cluster_sizes else 0.0
            cluster_size_std = np.std(cluster_sizes) if cluster_sizes else 0.0
            
            # Average edges per gene metrics
            avg_cluster_avg_edges_per_gene = np.mean(cluster_densities) if cluster_densities else 0.0
            # Weighted average (by cluster size)
            total_cluster_size = sum(cluster_sizes)
            weighted_avg_cluster_avg_edges_per_gene = (
                sum(c.get('avg_edges_per_gene', c.get('density', 0)) * c['size'] for c in clusters) / total_cluster_size
                if total_cluster_size > 0 else 0.0
            )
            
            # Largest cluster metrics
            largest_cluster_genes = largest_cluster['n_genes']
            largest_cluster_terms = largest_cluster['n_terms']
            largest_cluster_edges = largest_cluster['n_edges']
            largest_cluster_avg_edges_per_gene = largest_cluster.get('avg_edges_per_gene', largest_cluster.get('density', 0))
            
            # Library diversity metrics (only if requested and libraries are tracked)
            if include_library_diversity and self.term_to_library:
                avg_cluster_libraries = np.mean(cluster_libraries) if cluster_libraries else 0.0
                largest_cluster_libraries = largest_cluster['n_libraries']
            else:
                avg_cluster_libraries = None
                largest_cluster_libraries = None
            
            # Fraction of network in largest cluster
            fraction_in_largest_cluster = largest_cluster['size'] / (n_genes + n_terms) if (n_genes + n_terms) > 0 else 0.0
            fraction_edges_in_largest_cluster = largest_cluster['n_edges'] / n_edges if n_edges > 0 else 0.0
            
            # Hub genes (genes connecting to ≥3 terms)
            hub_genes_count = sum(1 for deg in gene_centralities if deg >= 3)
        else:
            avg_cluster_size = 0.0
            median_cluster_size = 0.0
            cluster_size_std = 0.0
            avg_cluster_avg_edges_per_gene = 0.0
            weighted_avg_cluster_avg_edges_per_gene = 0.0
            largest_cluster_genes = 0
            largest_cluster_terms = 0
            largest_cluster_edges = 0
            largest_cluster_avg_edges_per_gene = 0.0
            fraction_in_largest_cluster = 0.0
            fraction_edges_in_largest_cluster = 0.0
            hub_genes_count = 0
            
            # Library diversity (only if requested)
            if include_library_diversity:
                avg_cluster_libraries = 0.0
                largest_cluster_libraries = 0
            else:
                avg_cluster_libraries = None
                largest_cluster_libraries = None
        
        # Clustering coefficient (for bipartite graphs)
        clustering_coefficient = self._compute_bipartite_clustering()
        
        return {
            # Basic counts
            "n_genes": n_genes,
            "n_terms": n_terms,
            "n_edges": n_edges,
            
            # Gene connectivity
            "avg_connections_per_gene": avg_connections_per_gene,
            "gene_centrality_mean": gene_centrality_mean,
            "gene_centrality_std": gene_centrality_std,
            "gene_centrality_max": gene_centrality_max,
            "hub_genes_count": hub_genes_count,
            
            # Term connectivity
            "avg_connections_per_term": avg_connections_per_term,
            "term_centrality_mean": term_centrality_mean,
            "term_centrality_std": term_centrality_std,
            "term_centrality_max": term_centrality_max,
            
            # Cluster metrics (primary)
            "n_connected_components": n_components,
            "largest_component_size": largest_component_size,
            "largest_cluster_genes": largest_cluster_genes,
            "largest_cluster_terms": largest_cluster_terms,
            "largest_cluster_edges": largest_cluster_edges,
            "largest_cluster_avg_edges_per_gene": largest_cluster_avg_edges_per_gene,
            "avg_cluster_size": avg_cluster_size,
            "median_cluster_size": median_cluster_size,
            "cluster_size_std": cluster_size_std,
            "avg_cluster_avg_edges_per_gene": avg_cluster_avg_edges_per_gene,
            "weighted_avg_cluster_avg_edges_per_gene": weighted_avg_cluster_avg_edges_per_gene,
            "fraction_in_largest_cluster": fraction_in_largest_cluster,
            "fraction_edges_in_largest_cluster": fraction_edges_in_largest_cluster,
            
            # Library diversity (only included if computed)
            "largest_cluster_libraries": largest_cluster_libraries,
            "avg_cluster_libraries": avg_cluster_libraries,
            
            # Global clustering
            "clustering_coefficient": clustering_coefficient,
        }
    
    def _compute_connected_components(self) -> Tuple[int, int, List[Dict]]:
        """
        Compute connected components (clusters) using BFS.
        
        Returns:
            (number of components, size of largest component, list of cluster info)
        """
        visited_genes = set()
        visited_terms = set()
        components = []
        
        def bfs(start_gene: str) -> Tuple[Set[str], Set[str]]:
            """BFS from a gene node. Returns (genes, terms) in component."""
            component_genes = set()
            component_terms = set()
            queue = [("gene", start_gene)]
            
            while queue:
                node_type, node = queue.pop(0)
                
                if node_type == "gene":
                    if node in visited_genes:
                        continue
                    visited_genes.add(node)
                    component_genes.add(node)
                    # Add all connected terms
                    for term in self.gene_to_terms[node]:
                        if term not in visited_terms:
                            queue.append(("term", term))
                else:  # term
                    if node in visited_terms:
                        continue
                    visited_terms.add(node)
                    component_terms.add(node)
                    # Add all connected genes
                    for gene in self.term_to_genes[node]:
                        if gene not in visited_genes:
                            queue.append(("gene", gene))
            
            return component_genes, component_terms
        
        # Find all components
        for gene in self.genes:
            if gene not in visited_genes:
                cluster_genes, cluster_terms = bfs(gene)
                if cluster_genes or cluster_terms:
                    # Count edges within this cluster
                    cluster_edges = sum(
                        len([t for t in self.gene_to_terms[g] if t in cluster_terms])
                        for g in cluster_genes
                    )
                    
                    # Count unique libraries in this cluster
                    cluster_libraries = set()
                    for term in cluster_terms:
                        if term in self.term_to_library:
                            cluster_libraries.add(self.term_to_library[term])
                    n_libraries = len(cluster_libraries)
                    
                    # Average edges per gene: Average gene connectivity
                    # Formula: edges / genes
                    # This represents: average number of connections per gene
                    # Example: If a cluster has 10 genes and 15 edges:
                    #   avg_edges_per_gene = 15/10 = 1.5 (each gene connects to 1.5 terms on average)
                    if len(cluster_genes) > 0:
                        avg_edges_per_gene = cluster_edges / len(cluster_genes)
                    else:
                        avg_edges_per_gene = 0.0
                    
                    components.append({
                        'genes': cluster_genes,
                        'terms': cluster_terms,
                        'size': len(cluster_genes) + len(cluster_terms),
                        'n_genes': len(cluster_genes),
                        'n_terms': len(cluster_terms),
                        'n_edges': cluster_edges,
                        'density': avg_edges_per_gene,  # Keep 'density' key for backward compatibility
                        'avg_edges_per_gene': avg_edges_per_gene,  # New key name
                        'n_libraries': n_libraries,
                    })
        
        n_components = len(components)
        largest_component_size = max(c['size'] for c in components) if components else 0
        
        return n_components, largest_component_size, components
    
    def get_clusters(self) -> List[Dict]:
        """
        Get all clusters with their statistics.
        
        Returns:
            List of cluster dictionaries, sorted by size (largest first)
        """
        _, _, clusters = self._compute_connected_components()
        
        # Sort by size (largest first)
        clusters_sorted = sorted(clusters, key=lambda x: x['size'], reverse=True)
        
        # Add cluster number
        for i, cluster in enumerate(clusters_sorted, start=1):
            cluster['cluster_number'] = i
        
        return clusters_sorted
    
    def _compute_bipartite_clustering(self) -> float:
        """
        Compute clustering coefficient for bipartite graph.
        
        For bipartite graphs, we use the bipartite clustering coefficient:
        Measures how many neighbors of a node are also connected to each other
        (through paths of length 2).
        
        Returns:
            Average bipartite clustering coefficient
        """
        if len(self.genes) == 0:
            return 0.0
        
        clustering_coeffs = []
        
        for gene in self.genes:
            neighbors = self.gene_to_terms[gene]
            if len(neighbors) < 2:
                clustering_coeffs.append(0.0)
                continue
            
            # Count pairs of neighbors (terms) that share at least one gene
            neighbor_pairs = 0
            connected_pairs = 0
            
            neighbor_list = list(neighbors)
            for i in range(len(neighbor_list)):
                for j in range(i + 1, len(neighbor_list)):
                    term1 = neighbor_list[i]
                    term2 = neighbor_list[j]
                    neighbor_pairs += 1
                    
                    # Check if terms share at least one gene (other than current gene)
                    genes1 = self.term_to_genes[term1] - {gene}
                    genes2 = self.term_to_genes[term2] - {gene}
                    if genes1 & genes2:  # Intersection
                        connected_pairs += 1
            
            if neighbor_pairs > 0:
                clustering_coeffs.append(connected_pairs / neighbor_pairs)
            else:
                clustering_coeffs.append(0.0)
        
        return np.mean(clustering_coeffs) if clustering_coeffs else 0.0
    
    def reset(self):
        """Reset the analyzer for a new network."""
        self.genes.clear()
        self.terms.clear()
        self.edges.clear()
        self.gene_to_terms.clear()
        self.term_to_genes.clear()
        self.term_to_library.clear()


def compute_connectivity_from_permutation_results(
    permutation_file: str,
    output_file: Optional[str] = None
) -> pd.DataFrame:
    """
    Compute connectivity metrics from permutation results.
    
    Args:
        permutation_file: Path to merged permutation results TSV
        output_file: Optional path to save results
        
    Returns:
        DataFrame with connectivity metrics per permutation
    """
    logger.info(f"Loading permutation results from {permutation_file}")
    df = pd.read_csv(permutation_file, sep='\t')
    logger.info(f"Loaded {len(df)} permutation result rows")
    
    # Group by permutation (filename) and gene list size
    analyzer = NetworkConnectivityAnalyzer()
    all_metrics = []
    
    for (filename, gene_list_size), group in df.groupby(['filename', 'gene_list_size']):
        # Reset analyzer for this permutation
        analyzer.reset()
        
        # Convert group to iGEA results format
        results = []
        for _, row in group.iterrows():
            # Parse genes removed
            genes_removed = row.get('Genes removed for next iteration', '')
            if isinstance(genes_removed, str):
                genes_list = [g.strip() for g in genes_removed.split(',') if g.strip()]
            else:
                genes_list = []
            
            results.append({
                'Term': row.get('Term', ''),
                'Iteration': row.get('Iteration', 1),
                'Library': row.get('Library', ''),  # Include library name
                'Genes removed for next iteration': genes_list,
            })
        
        # Add to analyzer
        analyzer.add_igea_results(results)
        
        # Compute metrics (for permutations, we don't need library diversity since they're combined)
        # But we need cluster metrics, so use include_library_diversity=False (it won't matter for combined)
        metrics = analyzer.compute_metrics(include_library_diversity=False)
        metrics['filename'] = filename
        metrics['gene_list_size'] = gene_list_size
        all_metrics.append(metrics)
        
        if len(all_metrics) % 100 == 0:
            logger.info(f"Processed {len(all_metrics)} permutations...")
    
    # Create DataFrame
    metrics_df = pd.DataFrame(all_metrics)
    
    if output_file:
        metrics_df.to_csv(output_file, sep='\t', index=False)
        logger.info(f"Saved connectivity metrics to {output_file}")
    
    return metrics_df


def build_null_distribution(
    metrics_df: pd.DataFrame,
    output_file: Optional[str] = None
) -> Dict:
    """
    Build null distribution of connectivity metrics stratified by gene list size.
    
    Args:
        metrics_df: DataFrame with connectivity metrics
        output_file: Optional path to save null distribution JSON
        
    Returns:
        Dictionary with null distributions per gene list size
    """
    logger.info("Building null distribution...")
    
    null_distributions = {}
    
    for gene_list_size in sorted(metrics_df['gene_list_size'].unique()):
        subset = metrics_df[metrics_df['gene_list_size'] == gene_list_size]
        
        if len(subset) == 0:
            continue
        
        # Compute statistics for each metric
        # Priority order: Cluster metrics (primary) > Gene connectivity > Global metrics
        # Note: Library diversity metrics are optional and only included if present
        metrics_to_analyze = [
            # TIER 1: Cluster metrics (PRIMARY - focus on largest cluster)
            'largest_cluster_genes',           # Number of genes in largest cluster
            'largest_cluster_terms',           # Number of terms in largest cluster
            'largest_cluster_edges',           # Number of edges in largest cluster
            'largest_cluster_avg_edges_per_gene',  # Average edges per gene in largest cluster
            'largest_component_size',          # Total nodes in largest cluster
            'fraction_in_largest_cluster',     # Coverage of largest cluster
            'fraction_edges_in_largest_cluster', # Edge coverage in largest cluster
            # TIER 2: Cluster distribution
            'n_connected_components',         # Number of clusters (lower is better)
            'avg_cluster_size',                # Average cluster size
            'weighted_avg_cluster_avg_edges_per_gene',  # Size-weighted average edges per gene
            'avg_cluster_avg_edges_per_gene',  # Average edges per gene across clusters
            # TIER 3: Gene connectivity (within clusters)
            'hub_genes_count',                 # Genes connecting ≥3 terms
            'avg_connections_per_gene',         # Average terms per gene
            'gene_centrality_max',             # Max terms per gene
            # TIER 4: Global clustering
            'clustering_coefficient',          # Bipartite clustering coefficient
            # TIER 5: Library diversity (optional - only for combined networks)
            'largest_cluster_libraries',       # Number of libraries in largest cluster
            'avg_cluster_libraries',           # Average number of libraries per cluster
        ]
        
        distribution = {
            'gene_list_size': int(gene_list_size),
            'n_permutations': len(subset),
        }
        
        for metric in metrics_to_analyze:
            # Skip if metric doesn't exist (e.g., library diversity for individual libraries)
            if metric not in subset.columns:
                continue
            values = subset[metric].values
            # Filter out None/NaN values for optional metrics
            values = values[~pd.isna(values)]
            if len(values) == 0:
                continue
            distribution[metric] = {
                'mean': float(np.mean(values)),
                'std': float(np.std(values)),
                'median': float(np.median(values)),
                'q25': float(np.percentile(values, 25)),
                'q75': float(np.percentile(values, 75)),
                'q5': float(np.percentile(values, 5)),
                'q95': float(np.percentile(values, 95)),
                'min': float(np.min(values)),
                'max': float(np.max(values)),
            }
        
        null_distributions[int(gene_list_size)] = distribution
    
    if output_file:
        with open(output_file, 'w') as f:
            json.dump(null_distributions, f, indent=2)
        logger.info(f"Saved null distribution to {output_file}")
    
    return null_distributions


def benchmark_real_results(
    igea_results: List[Dict],
    gene_list_size: int,
    null_distribution: Dict,
    use_interpolation: bool = False
) -> Dict:
    """
    Benchmark real iGEA results against null distribution.
    
    Args:
        igea_results: List of iGEA result dictionaries
        gene_list_size: Size of the input gene list
        null_distribution: Null distribution dictionary
        use_interpolation: If True, interpolate between nearest sizes
        
    Returns:
        Dictionary with benchmark results (z-scores, percentiles, etc.)
    """
    # Find nearest gene list size in null distribution
    available_sizes = sorted([int(k) for k in null_distribution.keys()])
    if not available_sizes:
        return {}
    
    # Find nearest size(s)
    nearest_size = min(available_sizes, key=lambda x: abs(x - gene_list_size))
    size_diff = abs(gene_list_size - nearest_size)
    
    if use_interpolation and size_diff > 0 and len(available_sizes) > 1:
        # Find sizes above and below target for interpolation
        sizes_below = [s for s in available_sizes if s < gene_list_size]
        sizes_above = [s for s in available_sizes if s > gene_list_size]
        
        if sizes_below and sizes_above:
            # Use sizes on both sides for interpolation
            size_below = max(sizes_below)
            size_above = min(sizes_above)
            null_stats = _interpolate_null_distribution(
                gene_list_size,
                size_below,
                null_distribution[str(size_below)],
                size_above,
                null_distribution[str(size_above)]
            )
            logger.info(f"Interpolated null distribution for size {gene_list_size} between {size_below} and {size_above}")
        else:
            # Fall back to nearest size
            null_stats = null_distribution[str(nearest_size)]
            logger.warning(f"Could not interpolate (target size {gene_list_size} outside range), using nearest: {nearest_size}")
    else:
        null_stats = null_distribution[str(nearest_size)]
        if size_diff > 0:
            logger.info(f"Using null distribution for size {nearest_size} (difference: {size_diff} genes, {100*size_diff/gene_list_size:.1f}%) for gene list size {gene_list_size}")
    
    # Compute connectivity for real results
    analyzer = NetworkConnectivityAnalyzer()
    analyzer.add_igea_results(igea_results)
    real_metrics = analyzer.compute_metrics()
    
    # Compare against null distribution
    benchmark = {
        'gene_list_size': gene_list_size,
        'nearest_null_size': nearest_size,
        'size_difference': size_diff,
        'size_difference_pct': 100 * size_diff / gene_list_size if gene_list_size > 0 else 0,
        'real_metrics': real_metrics,
        'comparison': {},
    }
    
    # Metrics to compare, ordered by priority
    # Focus: Cluster metrics (primary) > Gene connectivity > Global metrics
    # Note: Library diversity metrics are optional and only compared if present
    metrics_to_compare = [
        # TIER 1: Cluster metrics (PRIMARY - focus on largest cluster)
        'largest_cluster_genes',           # Number of genes in largest cluster
        'largest_cluster_terms',           # Number of terms in largest cluster
        'largest_cluster_edges',           # Number of edges in largest cluster
        'largest_cluster_avg_edges_per_gene',  # Average edges per gene in largest cluster
        'largest_component_size',          # Total nodes in largest cluster
        'fraction_in_largest_cluster',     # Coverage of largest cluster
        'fraction_edges_in_largest_cluster', # Edge coverage in largest cluster
        # TIER 2: Cluster distribution
        'n_connected_components',         # Number of clusters (lower is better)
        'avg_cluster_size',                # Average cluster size
        'weighted_avg_cluster_avg_edges_per_gene',  # Size-weighted average edges per gene
        'avg_cluster_avg_edges_per_gene',  # Average edges per gene across clusters
        # TIER 3: Gene connectivity (within clusters)
        'hub_genes_count',                 # Genes connecting ≥3 terms
        'avg_connections_per_gene',         # Average terms per gene
        'gene_centrality_max',             # Max terms per gene
        # TIER 4: Global clustering
        'clustering_coefficient',          # Bipartite clustering coefficient
        # TIER 5: Library diversity (optional - only for combined networks)
        'largest_cluster_libraries',       # Number of libraries in largest cluster
        'avg_cluster_libraries',           # Average number of libraries per cluster
    ]
    
    for metric in metrics_to_compare:
        # Skip if metric not in real metrics or null distribution
        if metric not in real_metrics or metric not in null_stats:
            continue
        
        # Skip None values (for optional metrics like library diversity)
        if real_metrics.get(metric) is None:
            continue
        real_value = real_metrics.get(metric, 0)
        null_mean = null_stats.get(metric, {}).get('mean', 0)
        null_std = null_stats.get(metric, {}).get('std', 1)
        
        # Compute z-score
        if null_std > 0:
            z_score = (real_value - null_mean) / null_std
        else:
            z_score = 0.0
        
        # Compute percentile (assuming normal distribution)
        try:
            from scipy.stats import norm
            percentile = norm.cdf(z_score) * 100
        except ImportError:
            # Fallback: approximate percentile from z-score
            # Using standard normal approximation
            percentile = 50.0 + (z_score * 34.13)  # Rough approximation
        
        # Determine if higher or lower is better
        # Lower is better: n_connected_components (fewer clusters = better)
        # Higher is better: everything else
        if metric == 'n_connected_components':
            is_better = real_value < null_mean  # Lower is better
        else:
            is_better = real_value > null_mean  # Higher is better
        
        benchmark['comparison'][metric] = {
            'real_value': real_value,
            'null_mean': null_mean,
            'null_std': null_std,
            'z_score': z_score,
            'percentile': percentile,
            'is_better': is_better,
        }
    
    return benchmark


def _interpolate_null_distribution(
    target_size: int,
    size1: int,
    stats1: Dict,
    size2: int,
    stats2: Dict
) -> Dict:
    """
    Interpolate null distribution statistics between two sizes.
    
    Args:
        target_size: Target gene list size
        size1: First reference size
        stats1: Statistics for size1
        size2: Second reference size
        stats2: Statistics for size2
        
    Returns:
        Interpolated statistics dictionary
    """
    # Linear interpolation weight
    if size2 == size1:
        weight = 0.5
    else:
        weight = (target_size - size1) / (size2 - size1)
    
    # Interpolate each metric
    interpolated = {
        'gene_list_size': target_size,
        'n_permutations': int((1 - weight) * stats1.get('n_permutations', 0) + 
                              weight * stats2.get('n_permutations', 0)),
    }
    
    # Interpolate metric statistics (only metrics that exist in both)
    metrics = [
        'avg_connections_per_gene', 'gene_centrality_max', 'hub_genes_count',
        'n_connected_components', 'largest_component_size',
        'largest_cluster_genes', 'largest_cluster_terms', 'largest_cluster_edges',
        'largest_cluster_avg_edges_per_gene',
        'avg_cluster_size', 'avg_cluster_avg_edges_per_gene',
        'weighted_avg_cluster_avg_edges_per_gene',
        'fraction_in_largest_cluster',
        'fraction_edges_in_largest_cluster', 'clustering_coefficient',
        # Optional library diversity metrics
        'largest_cluster_libraries', 'avg_cluster_libraries'
    ]
    
    for metric in metrics:
        # Skip if metric doesn't exist in both stats (e.g., optional library diversity)
        if metric not in stats1 or metric not in stats2:
            continue
        stat1 = stats1[metric]
        stat2 = stats2[metric]
        
        interpolated[metric] = {
            'mean': (1 - weight) * stat1['mean'] + weight * stat2['mean'],
            'std': (1 - weight) * stat1['std'] + weight * stat2['std'],
            'median': (1 - weight) * stat1['median'] + weight * stat2['median'],
            'q25': (1 - weight) * stat1['q25'] + weight * stat2['q25'],
            'q75': (1 - weight) * stat1['q75'] + weight * stat2['q75'],
            'q5': (1 - weight) * stat1['q5'] + weight * stat2['q5'],
            'q95': (1 - weight) * stat1['q95'] + weight * stat2['q95'],
            'min': min(stat1['min'], stat2['min']),
            'max': max(stat1['max'], stat2['max']),
        }
    
    return interpolated


def benchmark_cluster(
    cluster: Dict,
    null_distribution: Dict,
    gene_list_size: int,
    use_interpolation: bool = False
) -> Dict:
    """
    Benchmark a single cluster against null distribution.
    
    Args:
        cluster: Cluster dictionary with metrics (from get_clusters())
        null_distribution: Null distribution dictionary
        gene_list_size: Size of the input gene list
        use_interpolation: If True, interpolate between nearest sizes
        
    Returns:
        Dictionary with benchmark statistics for the cluster
    """
    # Find nearest gene list size in null distribution
    available_sizes = sorted([int(k) for k in null_distribution.keys()])
    if not available_sizes:
        return {}
    
    nearest_size = min(available_sizes, key=lambda x: abs(x - gene_list_size))
    size_diff = abs(gene_list_size - nearest_size)
    
    if use_interpolation and size_diff > 0 and len(available_sizes) > 1:
        sizes_below = [s for s in available_sizes if s < gene_list_size]
        sizes_above = [s for s in available_sizes if s > gene_list_size]
        
        if sizes_below and sizes_above:
            size_below = max(sizes_below)
            size_above = min(sizes_above)
            null_stats = _interpolate_null_distribution(
                gene_list_size,
                size_below,
                null_distribution[str(size_below)],
                size_above,
                null_distribution[str(size_above)]
            )
        else:
            null_stats = null_distribution[str(nearest_size)]
    else:
        null_stats = null_distribution[str(nearest_size)]
    
    # Metrics to benchmark for this cluster
    cluster_metrics = {
        'cluster_size': cluster['size'],
        'cluster_genes': cluster['n_genes'],
        'cluster_terms': cluster['n_terms'],
        'cluster_edges': cluster['n_edges'],
        'cluster_avg_edges_per_gene': cluster.get('avg_edges_per_gene', cluster.get('density', 0)),
        'cluster_libraries': cluster.get('n_libraries', 0),
    }
    
    # Compare each metric against null distribution
    benchmark_results = {}
    
    # Map cluster metrics to null distribution metric names
    # Note: We compare each cluster against the "largest cluster" metrics from null
    # This tells us how this cluster compares to typical largest clusters in random networks
    # Some metrics may not be available in older null distributions
    metric_mapping = {
        'cluster_size': 'largest_component_size',  # Always available
        'cluster_genes': 'largest_cluster_genes',  # May not be in old null distributions
        'cluster_terms': 'largest_cluster_terms',  # May not be in old null distributions
        'cluster_edges': 'largest_cluster_edges',  # May not be in old null distributions
        'cluster_avg_edges_per_gene': 'largest_cluster_avg_edges_per_gene',  # May not be in old null distributions
        'cluster_libraries': 'largest_cluster_libraries',  # May not be in old null distributions
    }
    
    # Alternative: if cluster-specific metrics not available, we can't benchmark them
    # but we can still show the cluster size comparison
    
    for cluster_metric, null_metric in metric_mapping.items():
        real_value = cluster_metrics[cluster_metric]
        
        # Skip if metric not in null distribution (e.g., library diversity might not be present)
        if null_metric not in null_stats:
            continue
        
        null_mean = null_stats[null_metric].get('mean', 0)
        null_std = null_stats[null_metric].get('std', 1)
        
        # Compute z-score
        if null_std > 0:
            z_score = (real_value - null_mean) / null_std
        else:
            z_score = 0.0
        
        # Compute percentile
        try:
            from scipy.stats import norm
            percentile = norm.cdf(z_score) * 100
        except ImportError:
            percentile = 50.0 + (z_score * 34.13)  # Rough approximation
        
        # Determine if better (higher is better for all cluster metrics)
        is_better = real_value > null_mean
        
        benchmark_results[cluster_metric] = {
            'real_value': real_value,
            'null_mean': null_mean,
            'null_std': null_std,
            'null_min': null_stats[null_metric].get('min', 0),
            'null_max': null_stats[null_metric].get('max', 0),
            'null_median': null_stats[null_metric].get('median', 0),
            'null_n': null_stats[null_metric].get('n', 0),
            'z_score': z_score,
            'percentile': percentile,
            'is_better': is_better,
        }
    
    return benchmark_results


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Compute network connectivity benchmarks")
    parser.add_argument(
        "--permutation-file",
        type=str,
        default=str(PROJECT_ROOT / "archive" / "permutation_analysis" / "results" / "permutation_results" / "merged_permutation_results.tsv"),
        help="Path to merged permutation results"
    )
    parser.add_argument(
        "--output-metrics",
        type=str,
        default=str(PROJECT_ROOT / "results" / "connectivity_metrics.tsv"),
        help="Path to save connectivity metrics"
    )
    parser.add_argument(
        "--output-null",
        type=str,
        default=str(PROJECT_ROOT / "results" / "connectivity_null_distribution.json"),
        help="Path to save null distribution"
    )
    
    args = parser.parse_args()
    
    # Compute connectivity metrics from permutations
    metrics_df = compute_connectivity_from_permutation_results(
        args.permutation_file,
        args.output_metrics
    )
    
    # Build null distribution
    null_dist = build_null_distribution(
        metrics_df,
        args.output_null
    )
    
    print(f"\n✓ Computed connectivity metrics for {len(metrics_df)} permutations")
    print(f"✓ Built null distribution for {len(null_dist)} gene list sizes")
    print(f"\nGene list sizes in null distribution: {sorted(null_dist.keys())}")
