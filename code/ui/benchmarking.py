"""
Benchmarking functionality for Streamlit app.
Computes network connectivity benchmarks and displays statistical tables.
"""

import logging
import threading
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import pandas as pd

from network_connectivity_benchmark import (
    NetworkConnectivityAnalyzer,
    benchmark_cluster
)
from parallel_null_distribution import (
    find_intersection_libraries,
    compute_null_distribution_parallel,
    normalize_library_name
)

# Setup logging with file handler for debugging
logger = logging.getLogger(__name__)
# Add file handler if not already present
if not any(isinstance(h, logging.FileHandler) for h in logger.handlers):
    ROOT = Path(__file__).resolve().parent.parent.parent
    log_file = ROOT / "sandbox" / "benchmarking.log"
    log_file.parent.mkdir(exist_ok=True)  # Ensure sandbox directory exists
    file_handler = logging.FileHandler(str(log_file))
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
    logger.addHandler(file_handler)
    logger.setLevel(logging.DEBUG)


def format_status(z_score: float, percentile: float, is_available: bool = True) -> str:
    """
    Format status indicator based on z-score and percentile.
    
    Args:
        z_score: Z-score value
        percentile: Percentile value
        is_available: Whether the metric is available in null distribution
        
    Returns:
        Status string
    """
    if not is_available or (z_score == 0.0 and percentile == 50.0):
        return "Not available"
    
    if z_score > 2.0 and percentile > 95.0:
        return "✓✓ SIGNIFICANTLY BETTER"
    elif z_score > 1.0 and percentile > 84.0:
        return "✓ BETTER"
    elif z_score > -1.0 and percentile > 16.0:
        return "~ SIMILAR"
    elif z_score < -1.0 and percentile < 16.0:
        return "✗ WORSE"
    else:
        return "~ SIMILAR"


def compute_benchmark_for_streamlit(
    iter_enrich: Dict,
    gene_list_size: int,
    p_threshold: float,
    parquet_dir: Path,
    user_max_iterations: Optional[int] = None,
    merged_permutation_file: Optional[Path] = None
    ) -> Tuple[Optional[Dict], Optional[List[Dict]], List[str], List[str], int, Optional[NetworkConnectivityAnalyzer]]:
    """
    Compute benchmarking for Streamlit app.
    
    Args:
        iter_enrich: Dictionary of IterativeEnrichment objects by library name
        gene_list_size: Size of the input gene list
        p_threshold: P-value threshold used for enrichment
        parquet_dir: Directory containing Parquet files with permutation statistics
        merged_permutation_file: Optional path to merged permutation TSV file
        
    Returns:
        Tuple of (null_distribution, cluster_benchmarks, libraries_with_data, libraries_without_data, actual_size_used, analyzer)
        Returns None for null_distribution and cluster_benchmarks if benchmarking unavailable
        actual_size_used: The gene list size from permutation data that was actually used
        analyzer: NetworkConnectivityAnalyzer instance used for the analysis
    """
    # Get all library names from iter_enrich
    # NOTE: All user-selected libraries are shown in enrichment results.
    # However, for statistical benchmarking, we only use libraries that match
    # between user selection and available permutation data in parquet files.
    user_selected_libraries = list(iter_enrich.keys())
    
    # Find which user-selected libraries have matching permutation data
    # Statistics will only use libraries that match exactly between user selection and parquet files
    libraries_with_data, libraries_without_data = find_intersection_libraries(
        user_selected_libraries,
        parquet_dir,
        use_all_available=False  # Match user-selected libraries with available libraries
    )
    
    # Require minimum of 3 matching libraries for benchmarking
    # This ensures statistical robustness of the null distribution comparison
    MIN_MATCHING_LIBRARIES = 3
    if len(libraries_with_data) < MIN_MATCHING_LIBRARIES:
        logger.warning(
            f"Benchmarking requires at least {MIN_MATCHING_LIBRARIES} matching libraries. "
            f"Found {len(libraries_with_data)} matching libraries: {libraries_with_data}. "
            f"User selected: {user_selected_libraries}. "
            f"Skipping benchmarking."
        )
        return None, None, libraries_with_data, libraries_without_data, gene_list_size, None
    
    if not libraries_with_data:
        logger.warning("No libraries with permutation data available for benchmarking")
        return None, None, [], user_selected_libraries, gene_list_size, None
    
    # Require p-value threshold <= 0.01 for benchmarking
    # This ensures statistical rigor and matches the permutation data generation
    if p_threshold > 0.01:
        logger.warning(
            f"Benchmarking requires p-value threshold <= 0.01. "
            f"Your threshold: {p_threshold}. "
            f"Skipping benchmarking."
        )
        return None, None, libraries_with_data, libraries_without_data, gene_list_size, None
    
    # Compute null distribution in parallel
    null_dist_result = {
        'null_distribution': None,
        'status': 'running',
        'libraries_used': libraries_with_data,
        'error': None
    }
    null_dist_lock = threading.Lock()
    
    # Try to find merged permutation file
    if merged_permutation_file is None:
        # Try multiple possible locations
        possible_paths = [
            parquet_dir.parent / "merged_permutation_results.tsv",  # New location in permutations folder
            parquet_dir.parent.parent / "permutations" / "merged_permutation_results.tsv",
            parquet_dir.parent / "permutation_results" / "merged_permutation_results.tsv",
            parquet_dir.parent.parent / "archive" / "permutation_analysis" / "results" / "permutation_results" / "merged_permutation_results.tsv",
            parquet_dir.parent.parent / "merged_permutation_results.tsv",
            parquet_dir.parent.parent.parent / "merged_permutation_results.tsv",  # In case parquet_dir is nested deeper
        ]
        for path in possible_paths:
            if path.exists():
                merged_permutation_file = path
                logger.info(f"Found merged permutation file: {merged_permutation_file}")
                break
        else:
            logger.warning("Could not find merged permutation file, proceeding without p-value filtering")
    
    null_dist_thread = threading.Thread(
        target=compute_null_distribution_parallel,
        args=(
            parquet_dir,
            gene_list_size,
            libraries_with_data,
            null_dist_result,
            null_dist_lock,
            p_threshold,
            user_max_iterations,
            merged_permutation_file
        ),
        daemon=True
    )
    null_dist_thread.start()
    
    # Wait for null distribution computation
    logger.info(f"Waiting for null distribution computation (timeout: 120 seconds)...")
    null_dist_thread.join(timeout=120)  # Wait up to 2 minutes
    
    # Check if thread is still alive
    if null_dist_thread.is_alive():
        logger.warning(f"Thread is still alive after timeout - computation may still be running")
    else:
        logger.info(f"Thread completed (alive: {null_dist_thread.is_alive()})")
    
    # Check status
    final_status = null_dist_result.get('status', 'unknown')
    error_msg = null_dist_result.get('error', 'No error message')
    
    logger.debug(f"Null distribution computation status: {final_status}")
    logger.debug(f"Result dict keys: {list(null_dist_result.keys())}")
    logger.debug(f"Result dict: {null_dist_result}")
    
    if final_status == 'error':
        logger.error(f"Null distribution computation failed: {error_msg}")
        logger.error(f"Error details: {null_dist_result}")
        logger.error(f"Full result dict: {null_dist_result}")
        return None, None, libraries_with_data, libraries_without_data, gene_list_size, None
    elif final_status == 'unavailable':
        logger.warning(f"Statistical benchmarking unavailable: {error_msg}")
        logger.warning(f"Full result dict: {null_dist_result}")
        return None, None, libraries_with_data, libraries_without_data, gene_list_size, None
    elif final_status == 'running':
        logger.warning("Null distribution computation still running after timeout (120 seconds)")
        logger.warning(f"Result dict: {null_dist_result}")
        logger.warning(f"This may indicate the computation is taking longer than expected or is stuck")
        return None, None, libraries_with_data, libraries_without_data, gene_list_size, None
    elif final_status == 'completed':
        null_distribution = null_dist_result.get('null_distribution')
        actual_size_used = null_dist_result.get('null_distribution_size', gene_list_size)
        
        # null_distribution is stored as {str(actual_size): {metric_name: {mean, std, ...}}}
        if not null_distribution:
            logger.warning("Null distribution computation completed but result is empty")
            logger.warning(f"Result dict keys: {list(null_dist_result.keys())}")
            logger.warning(f"Null distribution value: {null_distribution}")
            return None, None, libraries_with_data, libraries_without_data, gene_list_size, None
        
        # Verify structure
        if not isinstance(null_distribution, dict):
            logger.warning(f"Unexpected null_distribution type: {type(null_distribution)}")
            return None, None, libraries_with_data, libraries_without_data, gene_list_size, None
        
        # Check if it has the expected structure {str(size): stats}
        # The structure should be: {"200": {metric_name: {mean, std, ...}}}
        size_keys = [k for k in null_distribution.keys() if str(k).isdigit()]
        if not size_keys:
            logger.warning(f"Null distribution does not have size keys. Keys: {list(null_distribution.keys())}")
            logger.warning(f"Expected structure: {{str(size): {{metric: stats}}}}")
            return None, None, libraries_with_data, libraries_without_data, gene_list_size, None
        
        # Verify the actual size is in the distribution
        actual_size_key = str(actual_size_used)
        if actual_size_key not in null_distribution:
            logger.warning(f"Size {actual_size_used} not found in null_distribution. Available sizes: {size_keys}")
            # Try to use the first available size
            if size_keys:
                actual_size_key = str(size_keys[0])
                actual_size_used = int(size_keys[0])
                logger.info(f"Using available size {actual_size_used} instead")
            else:
                return None, None, libraries_with_data, libraries_without_data, gene_list_size, None
        
        # Get stats for this size
        stats_for_size = null_distribution[actual_size_key]
        if not stats_for_size or len(stats_for_size) == 0:
            logger.warning(f"Null distribution stats for size {actual_size_used} is empty")
            return None, None, libraries_with_data, libraries_without_data, gene_list_size, None
        
        logger.info(f"Successfully computed null distribution for size {actual_size_used} with {len(stats_for_size)} metrics")
    else:
        logger.warning(f"Unexpected status: {final_status}")
        logger.warning(f"Result dict: {null_dist_result}")
        return None, None, libraries_with_data, libraries_without_data, gene_list_size, None
    
    # Build combined network from ONLY matching libraries (libraries_with_data)
    # This ensures fair comparison: real network uses same libraries as null distribution
    # NOTE: All user-selected libraries are shown in enrichment results (handled by UI layer).
    # Here we filter iGEA results to only use matching libraries for statistical benchmarking.
    # This ensures the comparison is valid: same libraries, same p-value threshold, same max iterations.
    # Need to match user library names (e.g., "Reactome") to parquet library names (e.g., "C2: CP: Reactome Pathways")
    combined_analyzer = NetworkConnectivityAnalyzer()
    
    # Create a mapping from user library names to parquet library names
    user_to_parquet_mapping = {}
    for user_lib in iter_enrich.keys():
        normalized_user = normalize_library_name(user_lib)
        for parquet_lib in libraries_with_data:
            normalized_parquet = normalize_library_name(parquet_lib)
            if (normalized_user == normalized_parquet or 
                normalized_user in normalized_parquet.lower() or 
                normalized_parquet in normalized_user.lower()):
                user_to_parquet_mapping[user_lib] = parquet_lib
                logger.info(f"Matched user library '{user_lib}' to parquet library '{parquet_lib}'")
                break
    
    if not user_to_parquet_mapping:
        logger.warning("No matching libraries found between analysis and permutation data")
        logger.warning(f"User libraries: {list(iter_enrich.keys())}")
        logger.warning(f"Parquet libraries: {libraries_with_data}")
        return None, None, libraries_with_data, libraries_without_data, gene_list_size, None
    
    libraries_to_use = list(user_to_parquet_mapping.keys())
    logger.info(f"Filtering iGEA results: Using {len(libraries_to_use)} libraries for benchmarking "
               f"(matched to {len(user_to_parquet_mapping)} parquet libraries):")
    for user_lib, parquet_lib in user_to_parquet_mapping.items():
        logger.info(f"  - {user_lib} → {parquet_lib}")
    
    # Determine iteration filter (cap at 30 to match permutation data)
    # Always cap at 30 iterations max (permutation data maximum)
    # If user specified max_iterations, use that (capped at 30)
    # If user didn't specify but has more than 30 iterations, cap at 30
    max_iter_filter = 30  # Always cap at 30 for fair comparison with permutation data
    if user_max_iterations is not None:
        max_iter_filter = min(user_max_iterations, 30)
        logger.info(f"Filtering iGEA results to iteration <= {max_iter_filter} (user requested {user_max_iterations}, capped at 30)")
    else:
        # Check if user's results have more than 30 iterations
        max_user_iter = 0
        for lib_name in libraries_to_use:
            iter_enrich_obj = iter_enrich.get(lib_name)
            if iter_enrich_obj and iter_enrich_obj.results:
                for record in iter_enrich_obj.results:
                    iteration = record.get("Iteration", 1)
                    if isinstance(iteration, (int, float)):
                        max_user_iter = max(max_user_iter, int(iteration))
        if max_user_iter > 30:
            logger.info(f"User's results have {max_user_iter} iterations, trimming to 30 for fair comparison with permutation data")
    
    for lib_name in libraries_to_use:
        iter_enrich_obj = iter_enrich.get(lib_name)
        if iter_enrich_obj is None:
            continue
        
        results = []
        for record in iter_enrich_obj.results:
            # Filter by iteration count (always cap at 30)
            iteration = record.get("Iteration", 1)
            # Convert to int if needed
            if isinstance(iteration, str):
                try:
                    iteration = int(iteration)
                except (ValueError, TypeError):
                    iteration = 1
            elif not isinstance(iteration, (int, float)):
                iteration = 1
            else:
                iteration = int(iteration)
            
            # Filter: iteration must be <= max_iter_filter (which is always <= 30)
            if iteration > max_iter_filter:
                continue
            
            genes_removed = record.get("Genes removed for next iteration", [])
            if isinstance(genes_removed, str):
                genes_removed = [g.strip() for g in genes_removed.split(",") if g.strip()]
            
            results.append({
                'Term': f"{lib_name}: {record.get('Term', '')}",
                'Iteration': iteration,
                'Library': lib_name,
                'Genes removed for next iteration': genes_removed,
            })
        
        combined_analyzer.add_igea_results(results)
    
    # Get clusters
    clusters = combined_analyzer.get_clusters()
    
    if not clusters:
        logger.info("No clusters found in network")
        # Get actual size used
        actual_size_used = gene_list_size
        if null_distribution:
            available_sizes = sorted([int(k) for k in null_distribution.keys()])
            if available_sizes:
                actual_size_used = min(available_sizes, key=lambda x: abs(x - gene_list_size))
        return null_distribution, [], libraries_with_data, libraries_without_data, actual_size_used, combined_analyzer
    
    # Benchmark only the largest cluster (null distribution is built from largest clusters only)
    # Clusters are already sorted by size (largest first)
    largest_cluster = clusters[0]
    
    logger.info(f"Benchmarking largest cluster: {largest_cluster['size']} nodes "
               f"({largest_cluster['n_genes']} genes, {largest_cluster['n_terms']} terms)")
    logger.info(f"Note: Only the largest cluster is benchmarked against null distribution "
               f"(which is built from largest clusters in each permutation)")
    
    cluster_benchmarks = []
    try:
        benchmark = benchmark_cluster(
            largest_cluster,
            null_distribution,
            gene_list_size,
            use_interpolation=True
        )
        
        if benchmark:
            cluster_benchmarks.append({
                'cluster': largest_cluster,
                'benchmark': benchmark
            })
        else:
            logger.warning(f"benchmark_cluster returned empty result for largest cluster "
                         f"with {largest_cluster.get('n_genes', 0)} genes, {largest_cluster.get('n_terms', 0)} terms")
    except Exception as e:
        logger.error(f"Error benchmarking largest cluster: {e}", exc_info=True)
    
    # Get the actual size used from null distribution
    actual_size_used = gene_list_size
    if null_distribution:
        available_sizes = sorted([int(k) for k in null_distribution.keys()])
        if available_sizes:
            actual_size_used = min(available_sizes, key=lambda x: abs(x - gene_list_size))
    
    return null_distribution, cluster_benchmarks, libraries_with_data, libraries_without_data, actual_size_used, combined_analyzer


def extract_benchmark_table_data(cluster_benchmarks: List[Dict]) -> List[Dict]:
    """
    Extract table data for "Statistical Benchmarks vs Random Gene Lists" table.
    
    Args:
        cluster_benchmarks: List of cluster benchmark dictionaries
        
    Returns:
        List of dictionaries with table row data (one per cluster, showing largest cluster only)
    """
    if not cluster_benchmarks:
        return []
    
    # For now, show only the largest cluster (first one, as they're sorted by size)
    # In the future, we could show all clusters or allow user to select
    largest_cluster_data = cluster_benchmarks[0]
    cluster = largest_cluster_data['cluster']
    benchmark = largest_cluster_data['benchmark']
    
    table_rows = []
    
    # Define metrics to display
    metrics = [
        ('Cluster Size', 'cluster_size', int),
        ('Number of Genes', 'cluster_genes', int),
        ('Number of Terms', 'cluster_terms', int),
        ('Number of Edges', 'cluster_edges', int),
        ('Average Edges per Gene', 'cluster_avg_edges_per_gene', float),
        ('Number of Libraries', 'cluster_libraries', int),
    ]
    
    for metric_name, metric_key, value_type in metrics:
        if metric_key not in benchmark:
            continue
        
        metric_data = benchmark[metric_key]
        real_value = metric_data['real_value']
        null_mean = metric_data.get('null_mean', 0.0)
        null_std = metric_data.get('null_std', 0.0)
        z_score = metric_data.get('z_score', 0.0)
        percentile = metric_data.get('percentile', 50.0)
        
        # Check if metric is available
        is_available = not (z_score == 0.0 and percentile == 50.0 and null_mean == 0.0)
        status = format_status(z_score, percentile, is_available)
        
        # Format value based on type
        if value_type == int:
            value_str = f"{int(real_value)}"
            mean_str = f"{null_mean:.2f}" if is_available and null_mean > 0 else "N/A"
            std_str = f"{null_std:.2f}" if is_available and null_std > 0 else "N/A"
        else:  # float
            value_str = f"{real_value:.4f}"
            mean_str = f"{null_mean:.4f}" if is_available and null_mean > 0 else "N/A"
            std_str = f"{null_std:.4f}" if is_available and null_std > 0 else "N/A"
        
        table_rows.append({
            'Metric': metric_name,
            'Value': value_str,
            'Null Mean': mean_str,
            'Null Std': std_str,
            'Z-score': f"{z_score:.2f}",
            'Percentile': f"{percentile:.1f}%",
            'Status': status
        })
    
    return table_rows


def generate_statistical_report_text(
    cluster_benchmarks: List[Dict],
    gene_list_name: str,
    libraries_with_data: List[str],
    libraries_without_data: List[str],
    analyzer: Optional[NetworkConnectivityAnalyzer] = None,
    gene_list_size: Optional[int] = None,
    actual_size_used: Optional[int] = None
) -> str:
    """
    Generate full statistical report text (same format as generate_cluster_statistical_report.py).
    
    Args:
        cluster_benchmarks: List of cluster benchmark dictionaries
        gene_list_name: Name of the gene list
        libraries_with_data: Libraries used for statistics
        libraries_without_data: Libraries excluded from statistics
        analyzer: Optional NetworkConnectivityAnalyzer to calculate term centralities
        gene_list_size: Original size of the input gene list
        actual_size_used: Gene list size from permutation data that was actually used for null distribution
        
    Returns:
        Full report text as string
    """
    lines = []
    
    # Header
    lines.append("=" * 100)
    lines.append("CLUSTER-BY-CLUSTER STATISTICAL REPORT")
    lines.append("=" * 100)
    lines.append("")
    lines.append(f"Gene List: {gene_list_name}")
    if gene_list_size is not None:
        lines.append(f"Gene List Size: {gene_list_size}")
    if actual_size_used is not None:
        if actual_size_used != gene_list_size:
            lines.append(f"Actual Size Used for Null Distribution: {actual_size_used} (rounded from {gene_list_size})")
        else:
            lines.append(f"Actual Size Used for Null Distribution: {actual_size_used}")
    lines.append(f"Total Clusters in Network: {len(cluster_benchmarks) if cluster_benchmarks else 0}")
    lines.append(f"Note: Only the largest cluster is benchmarked (null distribution built from largest clusters only)")
    lines.append("")
    
    # Library information
    lines.append("=" * 100)
    lines.append("IMPORTANT: Library Information for Statistical Analysis")
    lines.append("=" * 100)
    lines.append("")
    lines.append(f"Statistics were computed using {len(libraries_with_data)} libraries with permutation data:")
    lines.append(f"  {', '.join(libraries_with_data)}")
    lines.append("")
    if libraries_without_data:
        lines.append(f"Libraries included in enrichment but EXCLUDED from statistics:")
        lines.append(f"  {', '.join(libraries_without_data)}")
        lines.append(f"  (Permutation data not available for these libraries)")
        lines.append("")
    lines.append("The full network visualization includes all selected libraries.")
    lines.append("However, statistical benchmarks are only computed for libraries with available")
    lines.append("permutation data to ensure accurate comparison against the null distribution.")
    lines.append("")
    lines.append("=" * 100)
    lines.append("")
    lines.append("This report shows the largest cluster with benchmark statistics comparing")
    lines.append("against random gene lists of similar size.")
    lines.append("Note: Only the largest cluster is benchmarked, as the null distribution")
    lines.append("is built from the largest cluster in each permutation.")
    lines.append("")
    lines.append("Interpretation:")
    lines.append("  • Z-score > 2.0 and Percentile > 95%: Significantly better than random (top 5%)")
    lines.append("  • Z-score > 1.0 and Percentile > 84%: Better than random (top 16%)")
    lines.append("  • Z-score ~ 0.0 and Percentile ~ 50%: Similar to random")
    lines.append("  • Z-score < -1.0 and Percentile < 16%: Worse than random")
    lines.append("")
    lines.append("=" * 100)
    lines.append("")
    
    # Process each cluster
    for idx, cluster_data in enumerate(cluster_benchmarks, 1):
        cluster = cluster_data['cluster']
        benchmark = cluster_data['benchmark']
        
        lines.append("-" * 100)
        lines.append(f"LARGEST CLUSTER")
        lines.append("-" * 100)
        lines.append("")
        
        # Basic metrics
        lines.append("Basic Cluster Metrics:")
        lines.append(f"  Cluster Size:        {cluster['size']} nodes (genes + terms)")
        lines.append(f"  Number of Genes:     {cluster['n_genes']}")
        lines.append(f"  Number of Terms:     {cluster['n_terms']}")
        lines.append(f"  Number of Edges:     {cluster['n_edges']}")
        lines.append(f"  Average Edges per Gene: {cluster.get('avg_edges_per_gene', cluster.get('density', 0.0)):.4f}")
        lines.append(f"  Number of Libraries: {cluster.get('n_libraries', 0)}")
        if 'libraries' in cluster:
            lines.append(f"  Libraries:           {', '.join(cluster['libraries'])}")
        lines.append("")
        
        # Statistical benchmarks
        lines.append("Statistical Benchmarks vs Random Gene Lists:")
        lines.append("  Metric                    Value      Null Mean   Null Std    Z-score   Percentile   Status")
        lines.append("  " + "-" * 110)
        
        # Define metrics
        metrics = [
            ('Cluster Size', 'cluster_size', int),
            ('Number of Genes', 'cluster_genes', int),
            ('Number of Terms', 'cluster_terms', int),
            ('Number of Edges', 'cluster_edges', int),
            ('Average Edges per Gene', 'cluster_avg_edges_per_gene', float),
            ('Number of Libraries', 'cluster_libraries', int),
        ]
        
        for metric_name, metric_key, value_type in metrics:
            if metric_key not in benchmark:
                continue
            
            metric_data = benchmark[metric_key]
            real_value = metric_data['real_value']
            null_mean = metric_data.get('null_mean', 0.0)
            null_std = metric_data.get('null_std', 0.0)
            z_score = metric_data.get('z_score', 0.0)
            percentile = metric_data.get('percentile', 50.0)
            
            is_available = not (z_score == 0.0 and percentile == 50.0 and null_mean == 0.0)
            status = format_status(z_score, percentile, is_available)
            
            if value_type == int:
                if is_available and null_mean > 0:
                    lines.append(f"  {metric_name:<25} {int(real_value):>8}  {null_mean:>10.2f}  {null_std:>10.2f}  {z_score:>8.2f}  {percentile:>8.1f}%  {status}")
                else:
                    lines.append(f"  {metric_name:<25} {int(real_value):>8}  {'N/A':>10}  {'N/A':>10}  {z_score:>8.2f}  {percentile:>8.1f}%  {status}")
            else:  # float
                if is_available and null_mean > 0:
                    lines.append(f"  {metric_name:<25} {real_value:>8.4f}  {null_mean:>10.4f}  {null_std:>10.4f}  {z_score:>8.2f}  {percentile:>8.1f}%  {status}")
                else:
                    lines.append(f"  {metric_name:<25} {real_value:>8.4f}  {'N/A':>10}  {'N/A':>10}  {z_score:>8.2f}  {percentile:>8.1f}%  {status}")
        
        lines.append("")
        
        # Terms in cluster (ranked by centrality)
        if analyzer is not None and 'terms' in cluster:
            terms = list(cluster['terms'])
            if terms:
                # Calculate term centralities (number of genes connected to each term)
                term_data = []
                for term in terms:
                    # Get library for this term
                    library = analyzer.term_to_library.get(term, "Unknown")
                    # Extract term name (remove library prefix if present)
                    term_name = term.split(": ", 1)[-1] if ": " in term else term
                    # Calculate term centrality (degree = number of genes connected to this term)
                    centrality = len(analyzer.term_to_genes.get(term, set()))
                    term_data.append({
                        'term': f"{library}:{term_name}",
                        'centrality': centrality,
                        'original_term': term
                    })
                
                # Rank terms by centrality (highest first)
                term_data.sort(key=lambda x: x['centrality'], reverse=True)
                
                lines.append(f"Terms in Cluster ({len(term_data)} total, ranked by centrality):")
                lines.append("  Rank  Term                                                          Centrality (genes)")
                lines.append("  " + "-" * 90)
                for i, term_info in enumerate(term_data[:30], 1):  # Show top 30 terms
                    lines.append(f"  {i:2d}.   {term_info['term']:<70} {term_info['centrality']:>3d}")
                if len(term_data) > 30:
                    lines.append(f"  ... and {len(term_data) - 30} more terms")
                lines.append("")
    
    # Footer
    lines.append("=" * 100)
    lines.append("END OF REPORT")
    lines.append("=" * 100)
    
    return '\n'.join(lines)

