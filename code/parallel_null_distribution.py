#!/usr/bin/env python3
"""
Parallel null distribution computation from Parquet cluster statistics.

This module provides functions to:
1. Check which libraries have permutation data available
2. Compute null distribution from Parquet files in parallel with iGEA
3. Filter by user-selected libraries
"""

import pandas as pd
import numpy as np
import logging
from pathlib import Path
from typing import List, Dict, Set, Optional
from concurrent.futures import ThreadPoolExecutor
import threading

logger = logging.getLogger(__name__)


def get_available_libraries_from_parquet(parquet_dir: Path) -> Set[str]:
    """
    Get set of libraries that have permutation data available.
    
    Args:
        parquet_dir: Directory containing Parquet cluster statistics files
        
    Returns:
        Set of library names that have permutation data
    """
    if not parquet_dir.exists():
        logger.warning(f"Parquet directory not found: {parquet_dir}")
        return set()
    
    # Find any Parquet file to check library columns
    parquet_files = list(parquet_dir.glob("cluster_stats_size_*.parquet"))
    if not parquet_files:
        logger.warning(f"No Parquet files found in {parquet_dir}")
        return set()
    
    # Load one file to check available libraries
    sample_file = parquet_files[0]
    df = pd.read_parquet(sample_file)
    
    # Find boolean columns (has_*)
    bool_cols = [c for c in df.columns if c.startswith('has_')]
    
    # Also check the 'libraries' column to get actual library names
    all_libs_from_data = set()
    if 'libraries' in df.columns:
        for libs_str in df['libraries'].dropna():
            if isinstance(libs_str, str):
                libs = [lib.strip() for lib in libs_str.split(',')]
                all_libs_from_data.update(libs)
    
    # Convert column names back to library names
    available_libraries = set()
    for col in bool_cols:
        # has_go_bp -> GO BP
        # has_reactome -> Reactome
        # has_kegg -> KEGG
        lib_name = col.replace('has_', '').replace('_', ' ').title()
        # Handle special cases
        if 'go bp' in lib_name.lower():
            lib_name = 'GO BP'
        elif 'go cc' in lib_name.lower():
            lib_name = 'GO CC'
        elif 'go mf' in lib_name.lower():
            lib_name = 'GO MF'
        elif 'go ' in lib_name.lower():
            lib_name = lib_name.replace('Go ', 'GO ')
        elif 'kegg legacy' in lib_name.lower() or 'kegg_legacy' in lib_name.lower():
            lib_name = 'KEGG Legacy'
        elif 'kegg medicus' in lib_name.lower() or 'kegg_medicus' in lib_name.lower():
            lib_name = 'KEGG MEDICUS'
        elif 'kegg' in lib_name.lower():
            # Default to Legacy if just "KEGG" (most common case, and the one with permutation data)
            lib_name = 'KEGG Legacy'
        available_libraries.add(lib_name)
    
    # If we found libraries from the data, use those (more accurate)
    if all_libs_from_data:
        available_libraries = all_libs_from_data
    
    logger.info(f"Found {len(available_libraries)} libraries with permutation data: {sorted(available_libraries)}")
    return available_libraries


def normalize_library_name(lib_name: str) -> str:
    """
    Normalize library name for matching.
    Handles display names from alias.json like "H: Hallmark Gene Sets", "C5: Gene Ontology: Biological Process",
    and parquet library names like "C2: CP: Reactome Pathways".
    
    Since library names are now always from alias.json, this function normalizes those display names
    to match with parquet library names.
    
    Args:
        lib_name: Library name (should be from alias.json "name" field or parquet library name)
        
    Returns:
        Normalized library name for matching
    """
    lib_lower = lib_name.lower().strip()
    
    # Map common patterns to standard names (for display names from alias.json and parquet names)
    if 'hallmark' in lib_lower or lib_lower.startswith('h:'):
        return 'Hallmark'
    elif 'reactome' in lib_lower:
        return 'Reactome'
    elif 'kegg legacy' in lib_lower or 'kegg_legacy' in lib_lower:
        return 'KEGG Legacy'
    elif 'kegg medicus' in lib_lower or 'kegg_medicus' in lib_lower:
        return 'KEGG MEDICUS'
    elif 'kegg' in lib_lower:
        # Default to Legacy if just "KEGG" (most common case, and the one with permutation data)
        return 'KEGG Legacy'
    elif 'go bp' in lib_lower or 'biological process' in lib_lower or (lib_lower.startswith('c5:') and 'bp' in lib_lower):
        return 'GO BP'
    elif 'go cc' in lib_lower or 'cellular component' in lib_lower:
        return 'GO CC'
    elif 'go mf' in lib_lower or 'molecular function' in lib_lower:
        return 'GO MF'
    elif 'go ' in lib_lower and 'bp' in lib_lower:
        return 'GO BP'
    elif 'go ' in lib_lower and 'cc' in lib_lower:
        return 'GO CC'
    elif 'go ' in lib_lower and 'mf' in lib_lower:
        return 'GO MF'
    elif 'biocarta' in lib_lower:
        return 'BioCarta'
    elif 'wikipathways' in lib_lower:
        return 'WikiPathways'
    elif 'pathway interaction database' in lib_lower or 'pid' in lib_lower:
        return 'Pathway Interaction Database'
    elif 'canonical pathways' in lib_lower or (lib_lower.startswith('c2:') and 'canonical' in lib_lower):
        return 'Canonical pathways'
    
    # Remove common prefixes
    lib_clean = lib_lower
    for prefix in ['c2:', 'c5:', 'h:', 'c3:', 'c4:', 'c6:', 'c7:', 'c8:']:
        if lib_clean.startswith(prefix):
            lib_clean = lib_clean[len(prefix):].strip()
    
    # Also remove "cp:" prefix (common in parquet library names like "C2: CP: Reactome Pathways")
    if lib_clean.startswith('cp:'):
        lib_clean = lib_clean[len('cp:'):].strip()
    
    # Remove common suffixes
    for suffix in [' gene sets', ' pathways', ' pathway']:
        if lib_clean.endswith(suffix):
            lib_clean = lib_clean[:-len(suffix)].strip()
    
    # Title case
    return lib_clean.title()


def find_intersection_libraries(
    user_selected_libraries: List[str],
    parquet_dir: Path,
    use_all_available: bool = False
) -> tuple[List[str], List[str]]:
    """
    Find intersection between user-selected libraries and libraries with permutation data.
    
    Args:
        user_selected_libraries: List of library names selected by user
        parquet_dir: Directory containing Parquet cluster statistics files
        use_all_available: If True, use all available libraries from parquet files instead of filtering by user selection
        
    Returns:
        Tuple of (libraries_with_data, libraries_without_data)
    """
    available_libraries = get_available_libraries_from_parquet(parquet_dir)
    
    # If use_all_available is True, return all available libraries
    if use_all_available:
        libraries_with_data = sorted(list(available_libraries))
        libraries_without_data = []
        logger.info(f"Using all {len(libraries_with_data)} libraries available in parquet files: {libraries_with_data}")
        return libraries_with_data, libraries_without_data
    
    # Normalize available libraries
    available_normalized = {normalize_library_name(lib): lib for lib in available_libraries}
    
    # Find intersection
    libraries_with_data = []
    libraries_without_data = []
    
    for lib in user_selected_libraries:
        normalized = normalize_library_name(lib)
        
        # Try to find match
        matched = False
        for norm_key, orig_lib in available_normalized.items():
            # Check if normalized names match
            if normalized == norm_key:
                libraries_with_data.append(orig_lib)
                matched = True
                break
            # Also try direct substring match
            elif normalized in norm_key.lower() or norm_key.lower() in normalized.lower():
                libraries_with_data.append(orig_lib)
                matched = True
                break
        
        if not matched:
            libraries_without_data.append(lib)
    
    logger.info(f"Libraries with permutation data: {libraries_with_data}")
    if libraries_without_data:
        logger.warning(f"Libraries without permutation data (excluded from statistics): {libraries_without_data}")
    
    return libraries_with_data, libraries_without_data


def load_cluster_stats_for_size(
    parquet_dir: Path,
    gene_list_size: int,
    available_sizes: List[int]
) -> pd.DataFrame:
    """Load cluster statistics for the nearest gene list size."""
    # Round up to next increment of 50, cap at 1000
    original_size = gene_list_size
    if gene_list_size not in available_sizes:
        if gene_list_size > 1000:
            gene_list_size = 1000
            logger.warning(
                f"Gene list size {original_size} exceeds maximum permutation data size (1000). "
                f"Using size 1000 for null distribution comparison. "
                f"Note: Permutation data is only available up to size 1000."
            )
        else:
            # Round up to next multiple of 50
            gene_list_size = ((gene_list_size + 49) // 50) * 50
            logger.debug(
                f"Gene list size {original_size} not found, rounding up to nearest increment of 50: {gene_list_size}"
            )
        
        # Verify the rounded size exists
        if gene_list_size not in available_sizes:
            # Fallback: find nearest available size
            fallback_size = min(available_sizes, key=lambda x: abs(x - gene_list_size))
            logger.warning(
                f"Rounded size {gene_list_size} not found in Parquet data. "
                f"Using nearest available size: {fallback_size}"
            )
            gene_list_size = fallback_size
    
    parquet_file = parquet_dir / f"cluster_stats_size_{gene_list_size}.parquet"
    if not parquet_file.exists():
        raise FileNotFoundError(f"Parquet file not found: {parquet_file}")
    
    df = pd.read_parquet(parquet_file)
    logger.debug(f"Loaded {len(df):,} clusters from size {gene_list_size} (target: {original_size})")
    return df


def filter_by_libraries(
    df: pd.DataFrame,
    selected_libraries: List[str]
) -> pd.DataFrame:
    """
    Filter clusters by selected libraries.
    
    Includes clusters that contain at least one of the selected libraries.
    """
    if not selected_libraries:
        return df
    
    # Create filter mask: cluster must have at least one selected library
    mask = pd.Series([False] * len(df))
    
    for lib in selected_libraries:
        # Normalize library name for column matching
        col_name = f"has_{lib.lower().replace(' ', '_').replace(':', '_')}"
        if col_name in df.columns:
            mask = mask | df[col_name]
        else:
            # Fallback: check libraries string column
            mask = mask | df['libraries'].str.contains(lib, case=False, na=False)
    
    filtered_df = df[mask].copy()
    return filtered_df


def compute_null_distribution_from_raw_permutations(
    merged_permutation_file: Path,
    gene_list_size: int,
    selected_libraries: List[str],
    user_p_threshold: float,
    user_max_iterations: Optional[int] = None,
    metrics: List[str] = None
) -> Dict[str, Dict[str, float]]:
    """
    Compute null distribution from raw permutation results, filtering by user's p-value threshold and iteration count.
    
    Args:
        merged_permutation_file: Path to merged permutation results TSV file
        gene_list_size: Target gene list size
        selected_libraries: Libraries to filter by
        user_p_threshold: User's p-value threshold (must be <= 0.05)
        user_max_iterations: User's max iterations (permutation data will be filtered to this, capped at 30)
        metrics: List of metric names to compute
        
    Returns:
        Dictionary: {gene_list_size: {metric_name: {'mean': float, 'std': float, ...}}}
    """
    from network_connectivity_benchmark import NetworkConnectivityAnalyzer
    
    if user_p_threshold > 0.05:
        raise ValueError(f"User p-value threshold ({user_p_threshold}) exceeds 0.05. "
                        f"Statistical benchmarking is only available for p-value thresholds <= 0.05.")
    
    if metrics is None:
        metrics = [
            'largest_cluster_genes',
            'largest_cluster_terms',
            'largest_cluster_edges',
            'largest_cluster_avg_edges_per_gene',
            'largest_cluster_libraries',  # Added for library diversity benchmarking
            'largest_component_size',
            'n_connected_components',
            'avg_cluster_size',
            'avg_cluster_avg_edges_per_gene',
            'avg_cluster_libraries',  # Added for library diversity benchmarking
            'fraction_in_largest_cluster',
        ]
    
    logger.info(f"Loading raw permutation results from {merged_permutation_file}")
    df = pd.read_csv(merged_permutation_file, sep='\t')
    logger.info(f"Loaded {len(df):,} permutation result rows")
    
    # Find available gene list sizes
    available_sizes = sorted(df['gene_list_size'].unique())
    logger.info(f"Available gene list sizes in permutation data: {available_sizes}")
    
    # Find matching size: round UP to next increment of 50, cap at 1000
    original_size = gene_list_size
    if gene_list_size not in available_sizes:
        # Round up to next increment of 50
        if gene_list_size > 1000:
            # Cap at 1000 and warn
            gene_list_size = 1000
            logger.warning(
                f"Gene list size {original_size} exceeds maximum permutation data size (1000). "
                f"Using size 1000 for null distribution comparison. "
                f"Note: Permutation data is only available up to size 1000."
            )
        else:
            # Round up to next multiple of 50
            gene_list_size = ((gene_list_size + 49) // 50) * 50
            size_diff = gene_list_size - original_size
            logger.info(
                f"Gene list size {original_size} not found, rounding up to nearest increment of 50: "
                f"{gene_list_size} (difference: {size_diff})"
            )
        
        # Verify the rounded size exists in available sizes
        if gene_list_size not in available_sizes:
            # Fallback: find nearest available size
            nearest_size = min(available_sizes, key=lambda x: abs(x - gene_list_size))
            logger.warning(
                f"Rounded size {gene_list_size} not found in permutation data. "
                f"Using nearest available size: {nearest_size}"
            )
            gene_list_size = nearest_size
    
    # Filter by gene list size (now using nearest)
    df = df[df['gene_list_size'] == gene_list_size].copy()
    logger.info(f"After size filter: {len(df):,} rows for gene list size {gene_list_size}")
    
    if len(df) == 0:
        raise ValueError(f"No permutation data found for gene list size {gene_list_size}")
    
    # Filter by p-value threshold
    # Handle different possible column names for p-value
    p_value_col = None
    for col in ['iteration p-value', 'p-value', 'p_value', 'P-value', 'P_value', 'Full list p-value']:
        if col in df.columns:
            p_value_col = col
            break
    
    if p_value_col is None:
        raise ValueError("Could not find p-value column in permutation results. "
                        f"Available columns: {list(df.columns)}")
    
    # Convert p-value to float and filter
    df[p_value_col] = pd.to_numeric(df[p_value_col], errors='coerce')
    rows_before_pval = len(df)
    df = df[df[p_value_col] <= user_p_threshold].copy()
    logger.info(f"After p-value filter (<= {user_p_threshold}): {len(df):,} rows (was {rows_before_pval:,})")
    
    if len(df) == 0:
        raise ValueError(f"No permutation data found after p-value filtering (threshold: {user_p_threshold}). "
                        f"Permutation data was generated with p-value threshold 0.05, so very few results may pass stricter thresholds.")
    
    # Filter by selected libraries
    if selected_libraries:
        rows_before_lib = len(df)
        df = df[df['Library'].isin(selected_libraries)].copy()
        logger.info(f"After library filter: {len(df):,} rows (was {rows_before_lib:,}) for libraries: {selected_libraries}")
        
        if len(df) == 0:
            raise ValueError(f"No permutation data found after library filtering. "
                            f"Selected libraries: {selected_libraries}")
    
    # Filter by iteration count
    # Permutation data has max 30 iterations, cap user's filter at 30
    # Always apply iteration filter to match user's selection (or cap at 30 if user has more)
    if 'Iteration' in df.columns:
        rows_before_iter = len(df)
        df['Iteration'] = pd.to_numeric(df['Iteration'], errors='coerce')
        
        # Check iteration range in data
        if len(df) > 0:
            iter_min = df['Iteration'].min()
            iter_max = df['Iteration'].max()
            logger.info(f"Iteration range in data: {iter_min} to {iter_max}")
        
        if user_max_iterations is not None:
            # Cap at 30 (permutation data maximum)
            max_iter_filter = min(user_max_iterations, 30)
            df = df[df['Iteration'] <= max_iter_filter].copy()
            logger.info(f"After iteration filter (<= {max_iter_filter}): {len(df):,} rows (was {rows_before_iter:,}, user requested {user_max_iterations}, capped at 30)")
        else:
            # If user didn't specify, still cap at 30 (permutation data maximum)
            max_iter_filter = 30
            df = df[df['Iteration'] <= max_iter_filter].copy()
            logger.info(f"After iteration filter (<= {max_iter_filter}): {len(df):,} rows (was {rows_before_iter:,}, permutation data maximum)")
        
        if len(df) == 0:
            raise ValueError(f"No permutation data found after iteration filtering (max_iterations: {user_max_iterations or 'None'}, capped at 30). "
                            f"Data had iterations in range {iter_min} to {iter_max}.")
    else:
        logger.warning("Iteration column not found in permutation data, skipping iteration filter")
    
    # Get unique permutations
    unique_permutations = df.groupby('filename').size().reset_index()
    n_permutations = len(unique_permutations)
    logger.info(f"Processing {n_permutations} permutations")
    
    if n_permutations == 0:
        error_msg = (
            f"No permutation data available after all filtering steps. "
            f"Filters applied: size={gene_list_size}, p<={user_p_threshold}, "
            f"max_iterations={user_max_iterations if user_max_iterations is not None else 'None (capped at 30)'}, "
            f"libraries={selected_libraries}. "
            f"This may occur if: (1) p-value threshold is too strict (> 0.05), "
            f"(2) iteration filter is too restrictive, or (3) no matching libraries found."
        )
        raise ValueError(error_msg)
    
    # Initialize analyzer
    analyzer = NetworkConnectivityAnalyzer()
    
    # Process each permutation and compute cluster statistics
    all_cluster_stats = []
    for idx, (_, perm_info) in enumerate(unique_permutations.iterrows(), 1):
        filename = perm_info['filename']
        
        # Get data for this permutation
        perm_data = df[df['filename'] == filename]
        
        # Reset analyzer
        analyzer.reset()
        
        # Convert to iGEA results format
        results = []
        for _, row in perm_data.iterrows():
            # Try multiple possible column names for genes removed
            genes_removed = None
            for col in ['Genes removed for next iteration', 'Genes', 'genes_removed']:
                if col in row.index:
                    genes_removed = row.get(col, '')
                    break
            
            if genes_removed is None:
                genes_list = []
            elif isinstance(genes_removed, str):
                genes_list = [g.strip() for g in genes_removed.split(',') if g.strip()]
            elif isinstance(genes_removed, list):
                genes_list = genes_removed
            else:
                genes_list = []
            
            results.append({
                'Term': row.get('Term', ''),
                'Iteration': row.get('Iteration', 1),
                'Library': row.get('Library', ''),
                'Genes removed for next iteration': genes_list,
            })
        
        # Add to analyzer
        analyzer.add_igea_results(results)
        
        # Compute network-wide metrics (needed for fraction_in_largest_cluster)
        # This must be done before getting clusters to ensure metrics are available
        network_metrics = analyzer.compute_metrics(include_library_diversity=True)
        
        # Get clusters
        clusters = analyzer.get_clusters()
        
        # Store cluster statistics for this permutation
        if clusters:
            largest_cluster = clusters[0]  # Already sorted by size
            # Compute average cluster libraries
            cluster_libraries = [c.get('n_libraries', 0) for c in clusters]
            avg_cluster_libraries = np.mean(cluster_libraries) if cluster_libraries else 0.0
            
            all_cluster_stats.append({
                'filename': filename,
                'largest_cluster_genes': largest_cluster['n_genes'],
                'largest_cluster_terms': largest_cluster['n_terms'],
                'largest_cluster_edges': largest_cluster['n_edges'],
                'largest_cluster_avg_edges_per_gene': largest_cluster.get('avg_edges_per_gene', largest_cluster.get('density', 0)),
                'largest_cluster_libraries': largest_cluster.get('n_libraries', 0),
                'largest_component_size': largest_cluster['size'],
                'n_connected_components': len(clusters),
                'avg_cluster_size': np.mean([c['size'] for c in clusters]) if clusters else 0,
                'avg_cluster_avg_edges_per_gene': np.mean([c.get('avg_edges_per_gene', c.get('density', 0)) for c in clusters]) if clusters else 0,
                'avg_cluster_libraries': avg_cluster_libraries,
                'fraction_in_largest_cluster': network_metrics.get('fraction_in_largest_cluster', 0.0),
            })
        else:
            # No clusters - all metrics are 0
            all_cluster_stats.append({
                'filename': filename,
                'largest_cluster_genes': 0,
                'largest_cluster_terms': 0,
                'largest_cluster_edges': 0,
                'largest_cluster_avg_edges_per_gene': 0.0,
                'largest_cluster_libraries': 0,
                'largest_component_size': 0,
                'n_connected_components': 0,
                'avg_cluster_size': 0,
                'avg_cluster_avg_edges_per_gene': 0.0,
                'avg_cluster_libraries': 0.0,
                'fraction_in_largest_cluster': 0.0,
            })
        
        if idx % 100 == 0:
            logger.info(f"Processed {idx}/{n_permutations} permutations...")
    
    # Convert to DataFrame for easier statistics computation
    stats_df = pd.DataFrame(all_cluster_stats)
    
    # Compute statistics for each metric
    null_stats = {}
    for metric in metrics:
        if metric in stats_df.columns:
            values = stats_df[metric].dropna().values
            if len(values) > 0:
                null_stats[metric] = {
                    'mean': float(np.mean(values)),
                    'std': float(np.std(values)),
                    'n': int(len(values)),
                    'min': float(np.min(values)),
                    'max': float(np.max(values)),
                    'median': float(np.median(values))
                }
            else:
                logger.warning(f"No data for metric '{metric}'")
        else:
            logger.warning(f"Metric '{metric}' not found in cluster statistics")
    
    # The actual size used (may have been rounded)
    # Use the rounded gene_list_size, not the original
    actual_size = gene_list_size  # This is the rounded size after processing
    return null_stats, actual_size  # Return both stats and the actual size used


def compute_null_distribution_from_parquet(
    parquet_dir: Path,
    gene_list_size: int,
    selected_libraries: List[str],
    metrics: List[str] = None,
    user_p_threshold: float = None,
    user_max_iterations: Optional[int] = None,
    merged_permutation_file: Path = None
) -> Dict[str, Dict[str, float]]:
    """
    Compute null distribution statistics from Parquet cluster statistics.
    
    If user_p_threshold is provided and <= 0.05, filters raw permutation results by p-value.
    Otherwise, uses pre-computed Parquet cluster statistics (generated with p-value 0.05).
    
    Args:
        parquet_dir: Directory containing Parquet files
        gene_list_size: Target gene list size
        selected_libraries: Libraries to filter by
        metrics: List of metric names to compute (default: all cluster metrics)
        user_p_threshold: User's p-value threshold (if provided and <= 0.05, uses raw permutation filtering)
        merged_permutation_file: Path to merged permutation results TSV (required if user_p_threshold is provided)
        
    Returns:
        Dictionary: {gene_list_size: {metric_name: {'mean': float, 'std': float, ...}}}
    """
    # Check if user p-value threshold is provided and valid
    if user_p_threshold is not None:
        if user_p_threshold > 0.05:
            raise ValueError(
                f"User p-value threshold ({user_p_threshold}) exceeds 0.05. "
                f"Statistical benchmarking is only available for p-value thresholds <= 0.05. "
                f"Permutation data was generated with p-value threshold 0.05."
            )
        
        # If user threshold is exactly 0.05, use Parquet (faster, pre-computed)
        if user_p_threshold == 0.05:
            logger.info("User p-value threshold is 0.05, using pre-computed Parquet cluster statistics")
            # Continue to Parquet processing below
        else:
            # User threshold < 0.05, try to filter raw permutation results
            if merged_permutation_file is None:
                # Try to find merged permutation file
                project_root = parquet_dir.parent.parent
                possible_paths = [
                    project_root / "results" / "permutation_results" / "merged_permutation_results.tsv",
                    project_root / "archive" / "permutation_analysis" / "results" / "permutation_results" / "merged_permutation_results.tsv",
                    project_root / "merged_permutation_results.tsv",
                    parquet_dir.parent / "merged_permutation_results.tsv",
                ]
                
                merged_permutation_file = None
                for path in possible_paths:
                    if path.exists():
                        merged_permutation_file = path
                        logger.info(f"Found merged permutation file: {merged_permutation_file}")
                        break
            
            if merged_permutation_file is not None and merged_permutation_file.exists():
                logger.info(f"Using raw permutation results with p-value filtering (threshold: {user_p_threshold})")
                null_stats, actual_size = compute_null_distribution_from_raw_permutations(
                    merged_permutation_file,
                    gene_list_size,
                    selected_libraries,
                    user_p_threshold,
                    user_max_iterations,
                    metrics
                )
                # Return the actual size used (may have been rounded)
                return null_stats, actual_size
            else:
                # Merged TSV not available — fall back to pre-computed parquet stats
                logger.warning(
                    f"Merged permutation results file not found. "
                    f"Falling back to pre-computed Parquet cluster statistics (generated at p=0.05). "
                    f"Note: User threshold ({user_p_threshold}) cannot be applied without the merged TSV."
                )
                # Continue to Parquet processing below
    
    # Use pre-computed Parquet cluster statistics (generated with p-value 0.05)
    if metrics is None:
        metrics = [
            'largest_cluster_genes',
            'largest_cluster_terms',
            'largest_cluster_edges',
            'largest_cluster_avg_edges_per_gene',
            'largest_cluster_libraries',  # Added for library diversity benchmarking
            'largest_component_size',
            'n_connected_components',
            'avg_cluster_size',
            'avg_cluster_avg_edges_per_gene',
            'avg_cluster_libraries',  # Added for library diversity benchmarking
            'fraction_in_largest_cluster',
        ]
    
    # Find available sizes
    available_sizes = sorted([
        int(f.stem.split('_')[-1])
        for f in parquet_dir.glob("cluster_stats_size_*.parquet")
    ])
    
    if not available_sizes:
        raise FileNotFoundError(f"No Parquet files found in {parquet_dir}")
    
    # Load cluster statistics (this will round up to nearest increment of 50)
    original_size = gene_list_size
    df = load_cluster_stats_for_size(parquet_dir, gene_list_size, available_sizes)
    
    # Get the actual size used (may have been rounded)
    # The load_cluster_stats_for_size function rounds internally, so we need to determine
    # what size was actually loaded by checking what parquet file exists
    actual_size = original_size
    if original_size not in available_sizes:
        if original_size > 1000:
            actual_size = 1000
        else:
            actual_size = ((original_size + 49) // 50) * 50
        if actual_size not in available_sizes:
            actual_size = min(available_sizes, key=lambda x: abs(x - actual_size))
        if actual_size != original_size:
            logger.info(f"Using null distribution from size {actual_size} (rounded from {original_size})")
    
    # Filter by libraries
    if selected_libraries:
        df = filter_by_libraries(df, selected_libraries)
        logger.info(f"Filtered to {len(df):,} clusters containing selected libraries")
    
    # Compute statistics for each metric
    # All metrics need to be aggregated per permutation (one value per permutation)
    null_stats = {}
    
    # Helper function to get largest cluster metric per permutation
    def get_largest_cluster_metric(group, metric_col):
        """Get metric from largest cluster (cluster_number=1) for a permutation."""
        if len(group) == 0:
            return 0 if metric_col != 'density' else 0.0
        # cluster_number=1 is the largest cluster
        largest_clusters = group[group['cluster_number'] == 1]
        if len(largest_clusters) > 0:
            return largest_clusters.iloc[0][metric_col]
        # Fallback: get first cluster if no cluster_number=1
        return group.iloc[0][metric_col]
    
    for metric in metrics:
        values = None
        
        if metric == 'largest_cluster_genes':
            largest_per_perm = df.groupby('filename', group_keys=False).apply(
                lambda x: get_largest_cluster_metric(x, 'n_genes')
            )
            values = largest_per_perm.values
        elif metric == 'largest_cluster_terms':
            largest_per_perm = df.groupby('filename', group_keys=False).apply(
                lambda x: get_largest_cluster_metric(x, 'n_terms')
            )
            values = largest_per_perm.values
        elif metric == 'largest_cluster_edges':
            largest_per_perm = df.groupby('filename', group_keys=False).apply(
                lambda x: get_largest_cluster_metric(x, 'n_edges')
            )
            values = largest_per_perm.values
        elif metric == 'largest_cluster_avg_edges_per_gene':
            # Try new column name first, fallback to 'density' for backward compatibility
            largest_per_perm = df.groupby('filename', group_keys=False).apply(
                lambda x: get_largest_cluster_metric(x, 'avg_edges_per_gene') if 'avg_edges_per_gene' in x.columns else get_largest_cluster_metric(x, 'density')
            )
            values = largest_per_perm.values
        elif metric == 'largest_component_size':
            largest_per_perm = df.groupby('filename', group_keys=False).apply(
                lambda x: get_largest_cluster_metric(x, 'cluster_size')
            )
            values = largest_per_perm.values
        elif metric == 'n_connected_components':
            # Count clusters per permutation
            clusters_per_perm = df.groupby('filename').size()
            values = clusters_per_perm.values
        elif metric == 'avg_cluster_size':
            avg_per_perm = df.groupby('filename')['cluster_size'].mean()
            values = avg_per_perm.values
        elif metric == 'avg_cluster_avg_edges_per_gene':
            # Try new column name first, fallback to 'density' for backward compatibility
            if 'avg_edges_per_gene' in df.columns:
                avg_per_perm = df.groupby('filename')['avg_edges_per_gene'].mean()
            else:
                avg_per_perm = df.groupby('filename')['density'].mean()
            values = avg_per_perm.values
        elif metric == 'weighted_avg_cluster_avg_edges_per_gene':
            # Weighted by cluster size
            def weighted_density(group):
                if len(group) == 0:
                    return 0.0
                total_size = group['cluster_size'].sum()
                if total_size == 0:
                    return 0.0
                return (group['density'] * group['cluster_size']).sum() / total_size
            weighted_per_perm = df.groupby('filename', group_keys=False).apply(weighted_density)
            values = weighted_per_perm.values
        elif metric == 'fraction_in_largest_cluster':
            def compute_fraction(group):
                if len(group) == 0:
                    return 0.0
                largest = group[group['cluster_number'] == 1]
                if len(largest) == 0:
                    largest = group.iloc[[0]]
                total_genes = group['n_genes'].sum()
                total_terms = group['n_terms'].sum()
                total_nodes = total_genes + total_terms
                if total_nodes == 0:
                    return 0.0
                return largest.iloc[0]['cluster_size'] / total_nodes
            fractions = df.groupby('filename', group_keys=False).apply(compute_fraction)
            values = fractions.values
        elif metric == 'fraction_edges_in_largest_cluster':
            def compute_edge_fraction(group):
                if len(group) == 0:
                    return 0.0
                largest = group[group['cluster_number'] == 1]
                if len(largest) == 0:
                    largest = group.iloc[[0]]
                total_edges = group['n_edges'].sum()
                if total_edges == 0:
                    return 0.0
                return largest.iloc[0]['n_edges'] / total_edges
            edge_fractions = df.groupby('filename', group_keys=False).apply(compute_edge_fraction)
            values = edge_fractions.values
        elif metric == 'largest_cluster_libraries':
            largest_per_perm = df.groupby('filename', group_keys=False).apply(
                lambda x: get_largest_cluster_metric(x, 'n_libraries')
            )
            values = largest_per_perm.values
        elif metric == 'avg_cluster_libraries':
            avg_per_perm = df.groupby('filename')['n_libraries'].mean()
            values = avg_per_perm.values
        elif metric in ['cluster_size', 'n_genes', 'n_terms', 'n_edges', 'density', 'n_libraries']:
            # These are cluster-level metrics - aggregate per permutation (use largest cluster)
            largest_per_perm = df.groupby('filename', group_keys=False).apply(
                lambda x: get_largest_cluster_metric(x, metric)
            )
            values = largest_per_perm.values
        else:
            logger.warning(f"Metric '{metric}' not found and no computation method available")
            continue
        
        if values is not None and len(values) > 0:
            null_stats[metric] = {
                'mean': float(np.mean(values)),
                'std': float(np.std(values)),
                'n': int(len(values)),
                'min': float(np.min(values)),
                'max': float(np.max(values)),
                'median': float(np.median(values))
            }
        else:
            logger.warning(f"No data for metric '{metric}'")
    
    return null_stats, actual_size  # Return both stats and the actual size used


def compute_null_distribution_parallel(
    parquet_dir: Path,
    gene_list_size: int,
    selected_libraries: List[str],
    result_dict: Dict,
    lock: threading.Lock,
    user_p_threshold: float = None,
    user_max_iterations: Optional[int] = None,
    merged_permutation_file: Path = None
) -> None:
    """
    Compute null distribution in a separate thread and store result in shared dictionary.
    
    Args:
        parquet_dir: Directory containing Parquet files
        gene_list_size: Target gene list size
        selected_libraries: Libraries to filter by
        result_dict: Shared dictionary to store results
        lock: Thread lock for safe dictionary access
        user_p_threshold: User's p-value threshold (if > 0.05, benchmarking unavailable)
        merged_permutation_file: Path to merged permutation results TSV (for p-value filtering)
    """
    try:
        logger.info(f"Computing null distribution for size {gene_list_size} with libraries: {selected_libraries}")
        
        # Require p-value threshold <= 0.01 for benchmarking
        # This ensures statistical rigor and matches the permutation data generation
        if user_p_threshold is not None and user_p_threshold > 0.01:
            with lock:
                result_dict['status'] = 'unavailable'
                result_dict['error'] = (
                    f"Benchmarking requires p-value threshold <= 0.01. "
                    f"Your threshold: {user_p_threshold}. "
                    f"Permutation data was generated with p-value threshold 0.05, "
                    f"but benchmarking requires <= 0.01 for statistical rigor."
                )
            logger.warning(result_dict['error'])
            return
        
        if user_p_threshold is not None:
            logger.info(f"User p-value threshold: {user_p_threshold} (filtering permutation data)")
        else:
            logger.info(f"Using pre-computed cluster statistics (permutation data generated with 0.05)")
        
        # Compute null distribution (may round up the size)
        original_size = gene_list_size
        null_stats, actual_size = compute_null_distribution_from_parquet(
            parquet_dir,
            gene_list_size,
            selected_libraries,
            user_p_threshold=user_p_threshold,
            user_max_iterations=user_max_iterations,
            merged_permutation_file=merged_permutation_file
        )
        
        with lock:
            # Store with the actual size used (may be rounded up)
            result_dict['null_distribution'] = {str(actual_size): null_stats}
            result_dict['original_gene_list_size'] = original_size
            result_dict['null_distribution_size'] = actual_size
            result_dict['status'] = 'completed'
            result_dict['libraries_used'] = selected_libraries
            result_dict['permutation_p_threshold'] = 0.05  # Permutations always generated with 0.05
            if user_p_threshold is not None:
                result_dict['user_p_threshold'] = user_p_threshold
        
        logger.info(f"✓ Null distribution computation completed")
    except ValueError as e:
        # Handle p-value threshold errors
        logger.error(f"Error: {e}")
        with lock:
            result_dict['status'] = 'unavailable'
            result_dict['error'] = str(e)
    except Exception as e:
        logger.error(f"Error computing null distribution: {e}", exc_info=True)
        import traceback
        error_traceback = traceback.format_exc()
        logger.error(f"Full traceback:\n{error_traceback}")
        with lock:
            result_dict['status'] = 'error'
            result_dict['error'] = f"{str(e)}\n\nTraceback:\n{error_traceback}"
