import json
import logging
import math
import re
from io import StringIO
from pathlib import Path
from typing import Dict, List, Set
from datetime import datetime

import streamlit as st
import pandas as pd
from PIL import Image
from streamlit import session_state as state

# Existing imports
from background_gene_set import BackgroundGeneSet
from enrichment import Enrichment
from gene_set import GeneSet
from gene_set_library import GeneSetLibrary
from iter_enrichment import IterativeEnrichment
from ui.dot_utils import merge_iterative_dot
from ui.helpers import input_example, update_text_widgets, convert_and_validate_gene_input, display_conversion_results
from ui.processing import collect_results
from ui.rendering import (
    render_iter_results,
    render_network,
    render_results,
    render_validation,
    generate_regular_enrichment_json_analysis,
    render_ora_igea_comparison,
)
from ui.utils import download_link, download_file_link, update_aliases
from ui.benchmarking import (
    compute_benchmark_for_streamlit,
    extract_benchmark_table_data,
    generate_statistical_report_text
)

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)
ROOT = Path(__file__).resolve().parent.parent

st.set_page_config(
    page_title="Enrichment Analysis", layout="wide", initial_sidebar_state="expanded"
)


def get_msigdb_version_info() -> tuple[str, str]:
    """
    Get MSigDB version and last update date.
    
    Returns:
        Tuple of (version, update_date)
    """
    libraries_dir = ROOT / "data" / "libraries"
    backup_dir = libraries_dir / "backup"
    
    # Get version from GMT files
    version = "Unknown"
    gmt_files = list(libraries_dir.glob("*.gmt"))
    if gmt_files:
        # Extract version from first GMT file name
        sample_file = gmt_files[0].name
        if "v2025.1" in sample_file:
            version = "2025.1"
        elif "v2023.2" in sample_file:
            version = "2023.2"
        elif "v2023.1" in sample_file:
            version = "2023.1"
        else:
            # Try to extract version pattern
            version_match = re.search(r'v(\d{4}\.\d+)', sample_file)
            if version_match:
                version = version_match.group(1)
    
    # Get last update date from backup directories
    update_date = "Unknown"
    if backup_dir.exists():
        backup_dirs = [d for d in backup_dir.iterdir() if d.is_dir() and d.name.startswith("backup_")]
        if backup_dirs:
            # Sort by creation time and get the most recent
            latest_backup = max(backup_dirs, key=lambda x: x.stat().st_mtime)
            backup_name = latest_backup.name
            
            # Extract date from backup name (format: backup_YYYYMMDD_HHMMSS)
            if backup_name.startswith("backup_"):
                date_str = backup_name[7:]  # Remove "backup_" prefix
                try:
                    # Parse the date string
                    date_obj = datetime.strptime(date_str, "%Y%m%d_%H%M%S")
                    update_date = date_obj.strftime("%B %d, %Y")
                except ValueError:
                    update_date = "Unknown"
    
    return version, update_date


def _ensure_base_state():
    if "enrich" not in state:
        state.enrich = {}
    if "iter_results" not in state:
        state.iter_results: Dict[str, List[dict]] = {}
    if "iter_graph_parts" not in state:
        state.iter_graph_parts: Dict[str, dict] = {}
    if "results_ready" not in state:
        state.results_ready = False
    if "iter_min_overlap" not in state:
        state.iter_min_overlap = 3
    if "min_term_size" not in state:
        state.min_term_size = 10
    if "max_term_size" not in state:
        state.max_term_size = 600
    if "iter_max_term_size" not in state:
        state.iter_max_term_size = 600
    if "iter_ready" not in state:
        state.iter_ready = False
    if "selected_dot_paths" not in state:
        state.selected_dot_paths = []
    if "network_generated" not in state:
        state.network_generated = False  # flag to prevent clearing after checkbox changes
    if "libraries" not in state:
        state.libraries = []
    if "background_set" not in state:
        state.background_set = None
    if "gene_set_input" not in state:
        state.gene_set_input = ""
    if "gene_set_name" not in state:
        state.gene_set_name = ""
    if "selected_file" not in state:
        state.selected_file = "Select ..."
    if "gene_input_format" not in state:
        state.gene_input_format = 'symbols'
    if "bg_input_format" not in state:
        state.bg_input_format = 'symbols'
    if "advanced_settings_changed" not in state:
        state.advanced_settings_changed = False
    if "bt_submit_disabled" not in state:
        state.bt_submit_disabled = True
    if "bt_iter_disabled" not in state:
        state.bt_iter_disabled = True
    if "p_threshold" not in state:
        state.p_threshold = 0.01
    if "adjusted_p_threshold" not in state:
        state.adjusted_p_threshold = 0.05
    if "min_overlap" not in state:
        state.min_overlap = 3
    if "iter_p_threshold" not in state:
        state.iter_p_threshold = 0.01
    if "iter_max_iter" not in state:
        state.iter_max_iter = 10
    if "iter_min_term_size" not in state:
        state.iter_min_term_size = 10
    if "benchmark_computed" not in state:
        state.benchmark_computed = False
    if "benchmark_data" not in state:
        state.benchmark_data = None
    if "libraries_with_data" not in state:
        state.libraries_with_data = []
    if "libraries_without_data" not in state:
        state.libraries_without_data = []
    if "actual_size_used" not in state:
        state.actual_size_used = None
    if "benchmark_report_text" not in state:
        state.benchmark_report_text = None
    if "p_val_method" not in state:
        state.p_val_method = "Fisher's Exact Test"
    if "select_all_libraries" not in state:
        state.select_all_libraries = False


def reset_app() -> None:
    """Reset the app to default values."""
    logger.info("Resetting app to default values")
    
    # Clear ALL state variables to ensure complete reset
    for key in list(state.keys()):
        del state[key]
    
    # Reinitialize with default values
    _ensure_base_state()
    
    # Set default values for all parameters
    state.gene_input_format = 'symbols'
    state.bg_input_format = 'symbols'
    state.advanced_settings_changed = False
    state.bt_submit_disabled = True
    state.bt_iter_disabled = True
    
    # Reset all parameters to defaults
    state.p_threshold = 0.01
    state.min_overlap = 3
    state.min_term_size = 10
    state.max_term_size = 600
    state.iter_p_threshold = 0.01
    state.iter_max_iter = 10
    state.iter_min_overlap = 3
    state.iter_min_term_size = 10
    state.iter_max_term_size = 600
    state.p_val_method = "Fisher's Exact Test"
    
    # Clear all selections
    state.libraries = []
    state.background_set = None
    state.gene_set_input = ""
    state.gene_set_name = ""
    # Reset file selection dropdown to "Select ..."
    state.selected_file = "Select ..."
    
    # Clear any results and network state
    state.results_ready = False
    state.enrich = {}
    state.iter_enrich = {}
    state.iter_results.clear()
    state.iter_dot = {}
    state.selected_dot_paths = []
    state.network_generated = False
    state.last_merged_dot = ""
    
    # Clear any network checkbox states
    keys_to_remove = [key for key in state.keys() if key.startswith("use_") and key.endswith("_in_network")]
    for key in keys_to_remove:
        del state[key]
    
    keys_to_remove = [key for key in state.keys() if key.startswith("network_select_")]
    for key in keys_to_remove:
        del state[key]
    
    st.success("✅ App reset to default values!")


def _build_iterative_tables_download(all_iter_results: Dict[str, List[dict]]) -> str:
    """
    Build combined iterative enrichment results as TSV string.
    Uses the same format as to_dataframe() to match Streamlit UI columns.
    """
    import pandas as pd
    
    all_rows = []
    for lib, records in all_iter_results.items():
        for rec in records:
            # Get Library from record, fallback to dict key
            library_name = rec.get("Library", lib)
            
            # Get iteration p-value
            iter_pval = rec.get("iteration p-value", rec.get("p-value", ""))
            if isinstance(iter_pval, str) and iter_pval != "":
                try:
                    iter_pval_float = float(iter_pval)
                except ValueError:
                    iter_pval_float = 1.0
            elif isinstance(iter_pval, (int, float)):
                iter_pval_float = iter_pval
            else:
                iter_pval_float = 1.0
            
            # Calculate -log(p-value)
            iter_log_pval = -math.log10(iter_pval_float) if iter_pval_float > 0 else 0
            
            # Get genes removed for next iteration
            genes_removed = rec.get("Genes removed for next iteration", [])
            if isinstance(genes_removed, list):
                genes_removed_str = ', '.join(genes_removed)
            elif isinstance(genes_removed, str):
                genes_removed_str = genes_removed
            else:
                genes_removed_str = ""
            
            # Build row with all columns matching to_dataframe() format
            row = {
                "Library": library_name,
                "Iteration": rec.get("Iteration", ""),
                "Term": rec.get("Term", rec.get("term", "")),
                "Description": rec.get("Description", rec.get("description", "")),
                "iteration overlapping genes": rec.get("iteration overlapping genes", rec.get("Overlap size", "0/0")),
                "iteration p-value": iter_pval_float,
                "iteration -log(p-value)": iter_log_pval,
                "Genes removed for next iteration": genes_removed_str,
            }
            
            # Add optional columns if they exist
            if "Full list overlapping genes" in rec:
                row["Full list overlapping genes"] = rec["Full list overlapping genes"]
            if "Full list p-value" in rec:
                row["Full list p-value"] = rec["Full list p-value"]
            if "Regular FDR" in rec:
                row["Regular FDR"] = rec["Regular FDR"]
            if "Full list overlapping genes (gene names)" in rec:
                full_list_genes = rec["Full list overlapping genes (gene names)"]
                if isinstance(full_list_genes, list):
                    row["Full list overlapping genes (gene names)"] = ', '.join(full_list_genes)
                else:
                    row["Full list overlapping genes (gene names)"] = str(full_list_genes) if full_list_genes else ""
            
            all_rows.append(row)
    
    # Create DataFrame with proper column order
    if not all_rows:
        return ""
    
    df = pd.DataFrame(all_rows)
    
    # Define column order matching to_dataframe()
    expected_columns = [
        "Library", "Iteration", "Term", "Description", 
        "iteration overlapping genes", "iteration p-value", "iteration -log(p-value)",
        "Genes removed for next iteration",
        "Full list overlapping genes", "Full list p-value", "Regular FDR",
        "Full list overlapping genes (gene names)"
    ]
    
    # Reorder columns to match expected order
    existing_columns = [col for col in expected_columns if col in df.columns]
    other_columns = [col for col in df.columns if col not in expected_columns]
    column_order = existing_columns + other_columns
    
    df = df[column_order]
    
    return df.to_csv(sep="\t", index=False)


def _extract_iteration_1_as_ora(iterative_enrichments: Dict[str, IterativeEnrichment]) -> Dict[str, List[dict]]:
    """
    Extract iteration 1 results from iterative enrichment and format as ORA.
    
    Args:
        iterative_enrichments: Dictionary mapping library names to IterativeEnrichment objects
        
    Returns:
        Dictionary mapping library names to list of ORA-formatted results (iteration 1 only)
    """
    ora_results = {}
    
    for lib_name, iter_enrich in iterative_enrichments.items():
        if iter_enrich and iter_enrich.results:
            # Get iteration 1 records
            iteration_1_records = [r for r in iter_enrich.results if r.get("Iteration") == 1]
            
            # Convert to ORA format (similar to regular enrichment format)
            # Use the _regular_enrichment object if available for accurate ORA data
            ora_formatted = []
            
            # Try to get from _regular_enrichment first (most accurate)
            if hasattr(iter_enrich, '_regular_enrichment') and iter_enrich._regular_enrichment:
                regular_enrich = iter_enrich._regular_enrichment
                for rank, result in enumerate(regular_enrich.results, start=1):
                    overlap_genes = result.get("overlap", [])
                    if isinstance(overlap_genes, str):
                        genes_list = [g.strip() for g in overlap_genes.split(",") if g.strip()]
                    elif isinstance(overlap_genes, list):
                        genes_list = overlap_genes
                    else:
                        genes_list = []
                    
                    ora_formatted.append({
                        "Library": lib_name,
                        "Rank": rank,
                        "Term": result.get("term", ""),
                        "Description": result.get("description", ""),
                        "Overlap size": result.get("overlap_size", "0/0"),
                        "Genes": ", ".join(genes_list) if genes_list else "",
                        "p-value": result.get("p-value", 1.0),
                        "-log(p-value)": -math.log10(result.get("p-value", 1.0)) if result.get("p-value", 1.0) > 0 else 0,
                        "FDR": result.get("fdr", result.get("p-value", 1.0)),
                    })
            else:
                # Fallback: extract from iteration 1 records
                for rank, record in enumerate(iteration_1_records, start=1):
                    # Prefer "Full list overlapping genes (gene names)" if available (from merged data)
                    genes_list = []
                    full_list_genes = record.get("Full list overlapping genes (gene names)", [])
                    if full_list_genes:
                        if isinstance(full_list_genes, list):
                            genes_list = full_list_genes
                        elif isinstance(full_list_genes, str):
                            genes_list = [g.strip() for g in full_list_genes.split(",") if g.strip()]
                    
                    # Fallback to "Genes removed for next iteration"
                    if not genes_list:
                        genes_removed = record.get("Genes removed for next iteration", [])
                        if isinstance(genes_removed, str):
                            genes_list = [g.strip() for g in genes_removed.split(",") if g.strip()]
                        elif isinstance(genes_removed, list):
                            genes_list = genes_removed
                    
                    # Get overlap size - prefer "Full list overlapping genes" if available
                    overlap_size = record.get("Full list overlapping genes", record.get("iteration overlapping genes", "0/0"))
                    
                    # Get p-value - prefer "Full list p-value" if available
                    pval = record.get("Full list p-value", record.get("iteration p-value", record.get("p-value", 1.0)))
                    if isinstance(pval, str) and pval != "":
                        try:
                            pval = float(pval)
                        except ValueError:
                            pval = 1.0
                    
                    # Get FDR from merged data
                    fdr = record.get("Regular FDR", record.get("fdr", pval))
                    if isinstance(fdr, str) and fdr != "":
                        try:
                            fdr = float(fdr)
                        except ValueError:
                            fdr = pval
                    
                    ora_formatted.append({
                        "Library": lib_name,
                        "Rank": rank,
                        "Term": record.get("Term", record.get("term", "")),
                        "Description": record.get("Description", record.get("description", "")),
                        "Overlap size": overlap_size,
                        "Genes": ", ".join(genes_list) if genes_list else "",
                        "p-value": pval,
                        "-log(p-value)": -math.log10(pval) if pval > 0 else 0,
                        "FDR": fdr,
                    })
            
            if ora_formatted:
                ora_results[lib_name] = ora_formatted
    
    return ora_results


def _build_ora_tables_download(ora_results: Dict[str, List[dict]]) -> str:
    """
    Build combined ORA (iteration 1) results as TSV string.
    
    Args:
        ora_results: Dictionary mapping library names to list of ORA-formatted results
        
    Returns:
        TSV string with combined ORA results
    """
    rows = ["Library\tRank\tTerm\tDescription\tOverlap size\tp-value\t-log(p-value)\tFDR\tGenes"]
    for lib, records in ora_results.items():
        for rec in records:
            rows.append(
                f"{rec['Library']}\t{rec['Rank']}\t{rec['Term']}\t{rec['Description']}\t"
                f"{rec['Overlap size']}\t{rec['p-value']}\t{rec['-log(p-value)']}\t"
                f"{rec['FDR']}\t{rec['Genes']}"
            )
    return "\n".join(rows)


def _create_combined_ora_archive(iterative_enrichments: Dict[str, IterativeEnrichment]) -> str:
    """
    Create a single tar.gz archive containing all individual ORA (iteration 1) files.
    
    Args:
        iterative_enrichments: Dictionary mapping library names to IterativeEnrichment objects
        
    Returns:
        Path to the combined archive file
    """
    import tarfile
    from datetime import datetime
    import pandas as pd
    
    # Find the most recent run directory
    results_dir = ROOT / "results"
    results_dir.mkdir(exist_ok=True)
    
    if not results_dir.exists():
        logger.warning("No results directory found")
        return ""
    
    # Get all run directories and find the most recent one
    run_dirs = [d for d in results_dir.iterdir() if d.is_dir() and d.name.startswith("run_")]
    if not run_dirs:
        logger.warning("No run directories found")
        return ""
    
    # Sort by creation time and get the most recent
    latest_run_dir = max(run_dirs, key=lambda d: d.stat().st_ctime)
    logger.info(f"Using run directory: {latest_run_dir}")
    
    # Extract ORA results (iteration 1)
    ora_results = _extract_iteration_1_as_ora(iterative_enrichments)
    
    if not ora_results:
        logger.warning("No iteration 1 results found for ORA")
        return ""
    
    # Create combined archive filename
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    combined_archive_name = f"all_ora_files_{timestamp}.tar.gz"
    combined_archive_path = results_dir / combined_archive_name
    
    # Create tar.gz archive with all individual ORA files
    with tarfile.open(combined_archive_path, "w:gz") as tar:
        for lib_name, records in ora_results.items():
            # Create DataFrame and save as TSV
            df = pd.DataFrame(records)
            # Reorder columns
            column_order = ["Library", "Rank", "Term", "Description", "Overlap size", "p-value", "-log(p-value)", "FDR", "Genes"]
            df = df[[col for col in column_order if col in df.columns]]
            
            # Create temporary file for this library's ORA results
            import tempfile
            with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as tmp_file:
                df.to_csv(tmp_file.name, sep='\t', index=False)
                # Sanitize library name for filename
                safe_lib_name = lib_name.replace(':', '_').replace(' ', '_').replace('-', '_').replace('/', '_')
                arcname = f"{safe_lib_name}_ora_results.tsv"
                tar.add(tmp_file.name, arcname=arcname)
                logger.info(f"Added to archive: {arcname}")
                # Clean up temp file
                Path(tmp_file.name).unlink()
    
    logger.info(f"Created combined ORA files archive: {combined_archive_path}")
    return str(combined_archive_path)


def _create_combined_iteration_archive(iter_archives: Dict[str, str]) -> str:
    """
    Create a single tar.gz archive containing all individual iteration files directly.
    
    Args:
        iter_archives: Dictionary mapping library names to their archive paths (not used in new approach)
        
    Returns:
        Path to the combined archive file
    """
    import tarfile
    from datetime import datetime
    
    # Find the most recent run directory (the one created in this session)
    results_dir = ROOT / "results"
    # Ensure results directory exists
    results_dir.mkdir(exist_ok=True)
    
    if not results_dir.exists():
        logger.warning("No results directory found")
        return ""
    
    # Get all run directories and find the most recent one
    run_dirs = [d for d in results_dir.iterdir() if d.is_dir() and d.name.startswith("run_")]
    if not run_dirs:
        logger.warning("No run directories found")
        return ""
    
    # Sort by creation time and get the most recent
    latest_run_dir = max(run_dirs, key=lambda d: d.stat().st_ctime)
    logger.info(f"Using run directory: {latest_run_dir}")
    
    # Create combined archive filename
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    combined_archive_name = f"all_iteration_files_{timestamp}.tar.gz"
    combined_archive_path = results_dir / combined_archive_name
    
    # Create tar.gz archive with all individual files from the latest run
    with tarfile.open(combined_archive_path, "w:gz") as tar:
        # Find all individual iteration TSV files in the latest run directory
        tsv_files = list(latest_run_dir.glob("*_iteration_*.tsv"))
        
        # Add all TSV files
        for tsv_file in tsv_files:
            tar.add(tsv_file, arcname=tsv_file.name)
            logger.info(f"Added to archive: {tsv_file.name}")
    
    logger.info(f"Created combined iteration files archive: {combined_archive_path}")
    return str(combined_archive_path)


def _create_combined_regular_archive(regular_enrichments: Dict[str, Enrichment]) -> str:
    """
    Create a single tar.gz archive containing all individual regular enrichment files.
    
    Args:
        regular_enrichments: Dictionary mapping library names to Enrichment objects
        
    Returns:
        Path to the combined archive file
    """
    import tarfile
    from datetime import datetime
    
    # Find the most recent run directory (the one created in this session)
    results_dir = ROOT / "results"
    # Ensure results directory exists
    results_dir.mkdir(exist_ok=True)
    
    if not results_dir.exists():
        logger.warning("No results directory found")
        return ""
    
    # Get all run directories and find the most recent one
    run_dirs = [d for d in results_dir.iterdir() if d.is_dir() and d.name.startswith("run_")]
    if not run_dirs:
        logger.warning("No run directories found")
        return ""
    
    # Sort by creation time and get the most recent
    latest_run_dir = max(run_dirs, key=lambda d: d.stat().st_ctime)
    logger.info(f"Using run directory: {latest_run_dir}")
    
    # Create combined archive filename
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    combined_archive_name = f"all_regular_enrichment_files_{timestamp}.tar.gz"
    combined_archive_path = results_dir / combined_archive_name
    
    # Create tar.gz archive with all individual files from the latest run
    with tarfile.open(combined_archive_path, "w:gz") as tar:
        # Find all individual regular enrichment TSV files in the latest run directory
        # Regular enrichment files are named like: {library_name}_regular_results.tsv
        tsv_files = list(latest_run_dir.glob("*_regular_results.tsv"))
        
        # Also add JSON files if they exist
        json_files = list(latest_run_dir.glob("*.json"))
        
        # Add all TSV files
        for tsv_file in tsv_files:
            tar.add(tsv_file, arcname=tsv_file.name)
            logger.info(f"Added to archive: {tsv_file.name}")
        
        # Add JSON files (excluding iteration JSONs)
        for json_file in json_files:
            if "_iteration_" not in json_file.name:
                tar.add(json_file, arcname=json_file.name)
                logger.info(f"Added to archive: {json_file.name}")
    
    logger.info(f"Created combined regular enrichment files archive: {combined_archive_path}")
    return str(combined_archive_path)


def main() -> None:
    logger.info("Starting the Streamlit app")
    logo_path = ROOT / "code" / "static" / "logo.png"
    if logo_path.exists():
        st.sidebar.image(
            Image.open(logo_path),
            caption="Iterative Enrichment Analysis",
        )
    else:
        logger.warning(f"Logo not found at {logo_path}")
    st.sidebar.title("Enrichment analysis")
    st.sidebar.write(
        """This app tests the input gene list against selected pathway libraries. 
        
**Regular Mode**: Runs once and reports all results that pass the filters per library.

**Iterative Mode**: Iteratively removes genes from the top hit at each step until no term passes the p-value cutoff.

Results include ranked tables, bar charts, and network graphs."""
    )
    
    # Documentation button — opens rendered HTML on GitHub via htmlpreview.github.io
    import os
    github_owner = os.getenv("GITHUB_OWNER", "aion-labs")
    github_repo = os.getenv("GITHUB_REPO", "Gene-Enrichment-Analysis")
    github_branch = os.getenv("GITHUB_BRANCH", "streamlit-cloud")
    
    raw_url = f"https://raw.githubusercontent.com/{github_owner}/{github_repo}/{github_branch}/documentation/STREAMLIT_USER_GUIDE.html"
    file_url = f"https://htmlpreview.github.io/?{raw_url}"
    
    st.sidebar.markdown(
        f'''
        <div style="margin: 0.5rem 0;">
            <a href="{file_url}" target="_blank" style="text-decoration: none; display: block;">
                <button style="
                    background-color: #FF4B4B;
                    color: white;
                    padding: 0.75rem 1rem;
                    border: none;
                    border-radius: 0.5rem;
                    cursor: pointer;
                    width: 100%;
                    font-size: 1rem;
                    font-weight: 600;
                    transition: background-color 0.3s;
                    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                " onmouseover="this.style.backgroundColor='#FF6B6B'"
                   onmouseout="this.style.backgroundColor='#FF4B4B'">
                    📚 Documentation
                </button>
            </a>
        </div>
        ''',
        unsafe_allow_html=True
    )

    _ensure_base_state()

    mode = st.radio(
        "Mode",
        ["Regular", "Iterative"],
        index=1,
        horizontal=True,
        key="analysis_mode",
    )
    st.subheader(f"Enrichment analysis — {mode} mode")

    # Display MSigDB version and update date
    version, update_date = get_msigdb_version_info()
    if version != "Unknown" or update_date != "Unknown":
        st.caption(f"📊 **Gene Sets:** MSigDB v{version} • Last updated: {update_date}")

    state.lib_mapper = update_aliases("libraries")
    state.bg_mapper = update_aliases("backgrounds")
    state.advanced_settings_changed = False
    state.bt_submit_disabled = True

    analysis, advanced = st.tabs(["Analysis", "Advanced settings"])

    with analysis:
        col_input, col_settings = st.columns([5, 7])
        with col_input:
            # Gene input format selector
            if 'gene_input_format' not in state:
                state.gene_input_format = 'symbols'
            
            format_col1, format_col2 = st.columns([1, 3])
            with format_col1:
                state.gene_input_format = st.selectbox(
                    "Input Format",
                    ["symbols", "entrez_ids"],
                    index=0 if state.gene_input_format == 'symbols' else 1,
                    format_func=lambda x: "Gene Symbols" if x == "symbols" else "Entrez IDs"
                )
            
            with format_col2:
                if state.gene_input_format == 'symbols':
                    st.caption("💡 **Gene Symbols:** Enter official gene symbols (e.g., TP53, BRCA1)")
                else:
                    st.caption("💡 **Entrez IDs:** Enter numeric Entrez IDs (e.g., 7157, 672)")
            
            # Gene input area
            placeholder_text = (
                "Enter gene symbols (e.g., TP53, BRCA1) - one per line (max 800 genes)" 
                if state.gene_input_format == 'symbols' 
                else "Enter Entrez IDs (e.g., 7157, 672) - one per line (max 800 genes)"
            )
            
            st.text_area(
                "Input a set of genes",
                key="gene_set_input",
                height=400,
                placeholder=placeholder_text,
                label_visibility="collapsed",
            )
            st.caption("📝 **Note:** Maximum 800 genes allowed for optimal performance")
            
            # Auto-generate gene set name if not provided
            if not state.gene_set_name or state.gene_set_name.strip() == "":
                from datetime import datetime
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                state.gene_set_name = f"genelist_{timestamp}"
            
            st.text_input(
                "Input/edit gene list name:",
                key="gene_set_name",
                placeholder="Input a gene set name",
            )
            gene_files = [
                str(f).replace(f"{ROOT}/data/gene_lists/", "")
                for f in (ROOT / "data" / "gene_lists").rglob("*.txt")
            ]
            st.selectbox(
                "Or select a file from the `data` folder",
                ["Select ..."] + gene_files,
                index=0,
                on_change=update_text_widgets,
                key="selected_file",
            )
            
            # Add CSS for consistent button styling
            st.markdown("""
                <style>
                .stButton > button {
                    height: 40px;
                    margin-bottom: 8px;
                }
                </style>
                """, unsafe_allow_html=True)
            
            # Vertical button stack layout under the left column
            # Simplified activation: only require non-empty gene input
            ready_common = bool(state.gene_set_input and state.gene_set_input.strip())
            
            # Load Example button (top)
            st.button(
                "Load Example", 
                on_click=input_example,
                use_container_width=True,
                key="btn_load_example"
            )
            
            # Reset App button (middle)
            st.button(
                "Reset App", 
                on_click=reset_app, 
                type="secondary",
                use_container_width=True,
                key="btn_reset_app"
            )
            
            # Submit button (bottom)
            if mode == "Regular":
                state.bt_submit_disabled = not ready_common
                bt_submit = st.button(
                    "Submit", 
                    disabled=state.bt_submit_disabled, 
                    key="bt_reg",
                    use_container_width=True
                )
            else:
                state.bt_iter_disabled = not ready_common
                bt_iter = st.button(
                    "Submit",
                    disabled=state.bt_iter_disabled,
                    key="bt_iter",
                    use_container_width=True
                )
        with col_settings:
            # Always load default background
            if not hasattr(state, 'background_set') or not state.background_set:
                state.background_set = list(state.bg_mapper.keys())[0]  # Default to first background
            
            state.background_set = st.selectbox(
                "Background gene list", 
                state.bg_mapper.keys(),
                index=list(state.bg_mapper.keys()).index(state.background_set)
            )
            st.caption("Specifies the background list of genes...")
            
            # Initialize default libraries if not already set
            if not hasattr(state, 'libraries') or not state.libraries:
                default_libraries = [
                    "H: Hallmark Gene Sets",
                    "C2: Reactome Pathways", 
                    "C5: Gene Ontology: Biological Process"
                ]
                state.libraries = default_libraries
            
            # Simple library selection with 3 columns
            st.markdown("**Select libraries:**")
            
            # Get all available libraries
            all_libraries = list(state.lib_mapper.keys())
            
            # Create 3 columns for checkboxes with minimal spacing
            col1, col2, col3 = st.columns(3, gap="small")
            
            # Split libraries into 3 groups
            third = len(all_libraries) // 3
            libraries_col1 = all_libraries[:third]
            libraries_col2 = all_libraries[third:2*third]
            libraries_col3 = all_libraries[2*third:]
            
            # First column - compact layout
            with col1:
                for lib in libraries_col1:
                    is_selected = st.checkbox(
                        lib, 
                        value=lib in state.libraries,
                        key=f"lib_{lib.replace(':', '_').replace(' ', '_')}"
                    )
                    if is_selected and lib not in state.libraries:
                        state.libraries.append(lib)
                    elif not is_selected and lib in state.libraries:
                        state.libraries.remove(lib)
            
            # Second column - compact layout
            with col2:
                for lib in libraries_col2:
                    is_selected = st.checkbox(
                        lib, 
                        value=lib in state.libraries,
                        key=f"lib_{lib.replace(':', '_').replace(' ', '_')}"
                    )
                    if is_selected and lib not in state.libraries:
                        state.libraries.append(lib)
                    elif not is_selected and lib in state.libraries:
                        state.libraries.remove(lib)
            
            # Third column - compact layout
            with col3:
                for lib in libraries_col3:
                    is_selected = st.checkbox(
                        lib, 
                        value=lib in state.libraries,
                        key=f"lib_{lib.replace(':', '_').replace(' ', '_')}"
                    )
                    if is_selected and lib not in state.libraries:
                        state.libraries.append(lib)
                    elif not is_selected and lib in state.libraries:
                        state.libraries.remove(lib)
            
            # Load selected libraries
            if state.libraries:
                from ui.utils import get_library_name_from_alias
                state.gene_set_libraries = []
                for lib in state.libraries:
                    try:
                        lib_filename = state.lib_mapper.get(lib)
                        if not lib_filename:
                            logger.warning(f"Library '{lib}' not found in lib_mapper")
                            continue
                        lib_path = ROOT / "data" / "libraries" / lib_filename
                        if not lib_path.exists():
                            logger.error(f"Library file not found: {lib_path}")
                            st.error(f"❌ Library file not found: {lib_filename}")
                            continue
                        lib_name = get_library_name_from_alias(lib_path)
                        state.gene_set_libraries.append(
                            GeneSetLibrary(
                                str(lib_path),
                                name=lib_name,
                            )
                        )
                    except Exception as e:
                        logger.error(f"Error loading library '{lib}': {e}")
                        st.error(f"❌ Error loading library '{lib}': {e}")
            else:
                state.gene_set_libraries = []
            # Always load background gene set
            from ui.utils import get_background_info
            bg_file_path, bg_format = get_background_info(state.background_set)
            
            if bg_file_path:
                bg_path = Path(bg_file_path)
                if not bg_path.exists():
                    logger.error(f"Background file not found: {bg_path}")
                    st.error(f"❌ Background file not found: {bg_path.name}")
                    state.background_gene_set = None
                else:
                    try:
                        state.background_gene_set = BackgroundGeneSet(
                            bg_file_path,
                            input_format=bg_format  # Use format from alias file
                        )
                    except Exception as e:
                        logger.error(f"Error loading background: {e}")
                        st.error(f"❌ Error loading background: {e}")
                        state.background_gene_set = None
            else:
                st.error(f"❌ Could not load background: {state.background_set}")
                state.background_gene_set = None
            
            # Gene validation section (always runs)
            if state.gene_set_input:
                # Convert and validate gene input based on selected format
                converted_symbols, unrecognized_entrez, unrecognized_symbols, stats, conversions = convert_and_validate_gene_input(
                    state.gene_set_input, 
                    state.gene_input_format
                )
                
                # Display conversion results (always show validation results)
                display_conversion_results(converted_symbols, unrecognized_entrez, unrecognized_symbols, stats, state.gene_input_format, conversions)
                
                # Create gene set with converted symbols (if any valid genes found)
                if converted_symbols:
                    state.gene_set = GeneSet(
                        converted_symbols,
                        state.background_gene_set.genes,
                        state.gene_set_name,
                        hgcn=False,  # Skip background validation since genes are already validated
                        format=False,  # Skip formatting since genes are already formatted
                    )
                else:
                    state.gene_set = None
            if mode == "Regular":
                st.markdown("**Regular enrichment parameters**")
                st.caption("Configure filters for the enrichment analysis:")
                # P-value threshold filter for regular mode
                state.p_threshold = st.number_input(
                    "Raw p-value threshold",
                    min_value=1e-10,
                    max_value=0.5,
                    value=0.01,
                    step=0.001,
                    format="%.4f",
                    help="Maximum raw p-value for terms to be included in results."
                )
                # Adjusted p-value threshold filter for regular mode
                state.adjusted_p_threshold = st.number_input(
                    "Adjusted p-value threshold",
                    min_value=1e-10,
                    max_value=1.0,
                    value=0.05,
                    step=0.001,
                    format="%.4f",
                    help="Maximum adjusted p-value (FDR) for terms to be included in results."
                )
                # Minimum overlap filter for regular mode
                state.min_overlap = st.number_input(
                    "Minimum overlap with gene set",
                    min_value=1,
                    value=3,
                    step=1,
                    help="Minimum number of genes that must overlap between your input gene set and a gene set term."
                )
                # Term size filter for regular mode
                state.min_term_size, state.max_term_size = st.slider(
                    "Minimum and maximum term size",
                    min_value=1,
                    value=(10, 600),
                    step=10,
                    max_value=5000,
                    help="Filter gene sets by their size (number of genes)."
                )
            if mode == "Iterative":
                st.markdown("**Iterative parameters**")
                state.iter_p_threshold = st.number_input(
                    "P-value threshold",
                    min_value=1e-10,
                    max_value=0.5,
                    value=0.01,
                    step=0.001,
                    format="%.4f",
                )
                state.iter_max_iter = st.number_input(
                    "Max iterations (0 = no limit)",
                    min_value=0,
                    max_value=500,
                    value=10,
                    step=1,
                )
                state.iter_min_overlap = st.number_input(
                    "Minimum overlap with gene set",
                    min_value=1,
                    value=3,
                    step=1,
                )
                state.iter_min_term_size, state.iter_max_term_size = st.slider(
                    "Minimum and maximum term size",
                    min_value=1,
                    value=(10, 600),
                    step=10,
                    max_value=5000
                )


    # Filter gene sets based on mode-specific max_term_size
    if state.gene_set_libraries:
        if mode == "Regular":
            for gsl in state.gene_set_libraries:
                filtered_terms = [
                    t for t in gsl.library if t["size"] <= state.max_term_size
                ]
                gsl.library = filtered_terms
                gsl.num_terms = len(filtered_terms)
                gsl.unique_genes = gsl.compute_unique_genes()
                gsl.size = len(gsl.unique_genes)
        elif mode == "Iterative":
            for gsl in state.gene_set_libraries:
                filtered_terms = [
                    t for t in gsl.library if t["size"] <= state.iter_max_term_size
                ]
                gsl.library = filtered_terms
                gsl.num_terms = len(filtered_terms)
                gsl.unique_genes = gsl.compute_unique_genes()
                gsl.size = len(gsl.unique_genes)



    with advanced:
        if mode == "Regular":
            n_results = st.slider(
                "Number of results to display", 1, 100, 10, 1, key="n_res"
            )
        else:
            st.slider(
                "Number of results to display (regular only)",
                1,
                100,
                10,
                1,
                disabled=True,
            )
        # Use widget key to set session_state; do not assign to state directly
        st.selectbox(
            "P-value calculation method",
            ["Fisher's Exact Test", "Hypergeometric Test", "Chi-squared Test"],
            key="p_val_method",
        )
        if state.p_val_method != "Fisher's Exact Test":
            state.advanced_settings_changed = True
        # Background gene list format selector
        bg_format_col1, bg_format_col2 = st.columns([1, 3])
        with bg_format_col1:
            if 'bg_input_format' not in state:
                state.bg_input_format = 'symbols'
            
            state.bg_input_format = st.selectbox(
                "Background Format",
                ["symbols", "entrez_ids"],
                index=0 if state.bg_input_format == 'symbols' else 1,
                format_func=lambda x: "Gene Symbols" if x == "symbols" else "Entrez IDs"
            )
        
        with bg_format_col2:
            if state.bg_input_format == 'symbols':
                st.caption("💡 **Gene Symbols:** Upload file with gene symbols")
            else:
                st.caption("💡 **Entrez IDs:** Upload file with Entrez IDs (will be converted to symbols)")
        
        state.bg_custom = st.file_uploader(
            "Upload your background gene list", type=[".txt"]
        )
        if state.bg_custom:
            bgf = (ROOT / "data" / "backgrounds" / state.bg_custom.name).open("wb")
            bgf.write(state.bg_custom.getvalue())
            state.advanced_settings_changed = True
            
            # Update background aliases to refresh the menu
            state.bg_mapper = update_aliases("backgrounds")
            
            # Update the alias file to include format information for the uploaded file
            try:
                aliases_path = ROOT / "data" / "backgrounds" / "alias.json"
                with open(aliases_path, "r") as f:
                    aliases = json.load(f)
                
                # Find and update the entry for the uploaded file
                for entry in aliases:
                    if entry.get("file") == state.bg_custom.name:
                        entry["format"] = state.bg_input_format
                        break
                
                # Write back the updated aliases
                with open(aliases_path, "w") as f:
                    json.dump(aliases, f, indent=4)
                    
            except Exception as e:
                logger.warning(f"Could not update format in alias file: {e}")
            
            # Create background gene set with the uploaded file and selected format
            try:
                state.background_gene_set = BackgroundGeneSet(
                    str(ROOT / "data" / "backgrounds" / state.bg_custom.name),
                    name=state.bg_custom.name.replace(".txt", ""),
                    input_format=state.bg_input_format
                )
                st.success(f"✅ Background gene list loaded: {state.background_gene_set.size} genes")
                # Force page rerun to refresh the background menu
                st.rerun()
            except Exception as e:
                st.error(f"❌ Error loading background gene list: {str(e)}")
        # Gene set library upload disabled
        # state.libs_custom = st.file_uploader(
        #     "Upload gene set libraries",
        #     type=[".gmt"],
        #     accept_multiple_files=True,
        #     on_change=update_aliases,
        #     args=("libraries",),
        # )
        # if state.libs_custom:
        #     for libf in state.libs_custom:
        #         lf = (ROOT / "data" / "libraries" / libf.name).open("wb")
        #         lf.write(libf.getvalue())
        #         state.advanced_settings_changed = True
        # Apply settings button with improved layout
        col_apply, col_status = st.columns([1, 3])
        
        with col_apply:
            if state.advanced_settings_changed:
                if st.button("Apply Settings", use_container_width=True):
                    logger.info("Applied custom settings")
                    # Refresh aliases to ensure menus are updated
                    state.bg_mapper = update_aliases("backgrounds")
                    state.lib_mapper = update_aliases("libraries")
                    with col_status:
                        st.success("✅ Settings applied successfully!")
                    # Force page rerun to refresh the menus
                    st.rerun()
            else:
                st.button("Apply Settings", disabled=True, use_container_width=True)
        
        with col_status:
            if state.advanced_settings_changed:
                st.info("ℹ️ Settings have been modified. Click 'Apply Settings' to save changes.")
            else:
                st.success("✅ All settings are up to date.")

    # Regular execution
    if mode == "Regular" and "bt_submit" in locals() and bt_submit:
        logger.info("Running regular enrichment")
        # Clear previous results state
        state.results_ready = False
        state.enrich = {}
        render_validation()
        
        # Validation checks after submit
        validation_passed = True
        
        # Check 1: Gene set exists and meets size requirements
        if not hasattr(state, 'gene_set') or not state.gene_set:
            st.error("❌ **Gene validation failed!** No valid gene set created. Please check your gene input and try again.")
            validation_passed = False
        elif state.gene_set.size < 20:
            st.error(f"❌ **Gene validation failed!** Your input contains {state.gene_set.size} validated genes, but at least 20 genes are required. Please add more genes and try again.")
            validation_passed = False
        elif state.gene_set.size > 800:
            st.error(f"❌ **Gene validation failed!** Your input contains {state.gene_set.size} validated genes, but the maximum allowed is 800 genes. Please reduce your gene list size and try again.")
            validation_passed = False
        
        # Check 2: Background is valid
        if not state.background_gene_set:
            st.error("❌ **Background validation failed!** No valid background gene set loaded. Please try refreshing the page.")
            validation_passed = False
        
        # Check 3: At least 1 selected library
        if not state.gene_set_libraries or len(state.gene_set_libraries) == 0:
            st.error("❌ **Library validation failed!** No gene set libraries selected. Please select at least one library.")
            validation_passed = False
        
        # Only proceed if all validations pass
        if validation_passed and state.gene_set_input and ready_common:
            
            # Use validated gene count instead of raw input count
            n_genes = state.gene_set.size if state.gene_set else len(state.gene_set_input.split())
            if n_genes <= 100 or n_genes >= 5000:
                warn = "small" if n_genes <= 100 else "big"
                s = "s" if str(n_genes)[-1] != "1" else ""
                st.warning(f"You've entered {n_genes} validated gene{s}, which may be {warn}...")
            # Create progress container for regular enrichment
            progress_container = st.container()
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            total_libraries = len(state.gene_set_libraries)
            
            for i, gsl in enumerate(state.gene_set_libraries):
                # Update progress
                progress = (i + 1) / total_libraries
                progress_bar.progress(progress)
                status_text.text(f"Processing library {i+1}/{total_libraries}: {gsl.name} ({gsl.num_terms} terms)")
                
                logger.info(f"Starting regular enrichment for library: {gsl.name}")
                
                try:
                    enrich = Enrichment(
                        state.gene_set,
                        gsl,
                        state.background_gene_set,
                        min_term_size=state.min_term_size,
                        max_term_size=state.max_term_size,
                        p_value_method_name=state.p_val_method,
                    )
                    
                    # Filter results by p-value threshold, adjusted p-value threshold, and minimum overlap
                    filtered_results = [
                        result for result in enrich.results
                        if (result.get("p-value", 1.0) <= state.p_threshold and
                            result.get("fdr", 1.0) <= state.adjusted_p_threshold and
                            result.get("overlap_size", "").split("/")[0].isdigit() and 
                            int(result.get("overlap_size", "").split("/")[0]) >= state.min_overlap)
                    ]
                    
                    # Re-rank the filtered results (1-based indexing)
                    for i, result in enumerate(filtered_results):
                        result["rank"] = i + 1
                    
                    # Update the enrichment object with filtered results
                    enrich.results = filtered_results
                    state.enrich[gsl.name] = enrich
                    
                    # Ensure results directory exists
                    results_dir = ROOT / "results"
                    results_dir.mkdir(exist_ok=True)
                    with (results_dir / f"{enrich.name}.json").open("w") as js:
                        json.dump(enrich.to_snapshot(), js)
                    
                    logger.info(f"Completed regular enrichment for library: {gsl.name} ({len(enrich.results)} terms)")
                    
                except Exception as e:
                    logger.error(f"Error processing library {gsl.name}: {e}")
                    st.error(f"❌ Error processing library {gsl.name}: {str(e)}")
                    continue
            
            # Clear progress indicators
            progress_bar.empty()
            status_text.empty()
            
            state.results_ready = True
    if mode == "Regular" and state.results_ready:
        logger.info("Displaying regular results")
        
        # Download section with combined and all individual files
        st.markdown("**📥 Download Results:**")
        
        # 1. Combined regular enrichment TSV and JSON
        combined_tsv = collect_results(state.enrich)
        combined_json = generate_regular_enrichment_json_analysis(state.enrich)
        st.markdown(
            f"📊 **Combined Results:** {download_link(combined_tsv, 'regular_enrichment_results', 'tsv')}, "
            f"{download_link(combined_json, 'regular_enrichment_results', 'json')}",
            unsafe_allow_html=True,
        )
        
        # 2. All individual regular enrichment files as single archive
        combined_archive_path = _create_combined_regular_archive(state.enrich)
        if combined_archive_path:
            combined_archive_name = Path(combined_archive_path).name
            # Remove the .tar.gz extension since download_file_link will add it back
            base_filename = combined_archive_name.replace('.tar.gz', '')
            st.markdown(
                f"📁 **All Individual Files:** {download_file_link(combined_archive_path, base_filename, 'tar.gz')}",
                unsafe_allow_html=True,
            )
        
        # Display results for each library in the same order as the checkbox list
        def sort_libraries(libraries):
            hallmark = [lib for lib in libraries if lib.startswith("H:")]
            c2_libs = sorted([lib for lib in libraries if lib.startswith("C2:")])
            c5_libs = sorted([lib for lib in libraries if lib.startswith("C5:")])
            protein_interaction = [lib for lib in libraries if lib.startswith("Protein Interaction")]
            other = sorted([lib for lib in libraries if not any(lib.startswith(prefix) for prefix in ["H:", "C2:", "C5:", "Protein Interaction"])])
            
            return hallmark + c2_libs + c5_libs + protein_interaction + other
        
        available_libraries = sort_libraries(list(state.enrich.keys()))
        for lib in available_libraries:
            render_results(state.enrich[lib], lib, n_results)

    # Iterative execution
    if mode == "Iterative" and "bt_iter" in locals() and bt_iter:
        logger.info("Running iterative enrichment")
        # Clear previous results state
        state.iter_enrich = {}
        state.iter_results.clear()
        state.iter_dot = {}
        state.selected_dot_paths = []
        state.network_generated = False
        state.last_merged_dot = ""
        
        # Clear all network checkbox states from previous runs
        keys_to_remove = [key for key in state.keys() if key.startswith("use_") and key.endswith("_in_network")]
        for key in keys_to_remove:
            del state[key]
        
        render_validation()
        
        # Validation checks after submit
        validation_passed = True
        
        # Check 1: Gene set exists and meets size requirements
        if not hasattr(state, 'gene_set') or not state.gene_set:
            st.error("❌ **Gene input verfication failed!** No valid gene set created. Please check your gene input and try again.")
            validation_passed = False
        elif state.gene_set.size < 20:
            st.error(f"❌ **Gene input verfication failed!** Your input contains {state.gene_set.size} validated genes, but the minimum required is 20 genes. Please add more genes and try again.")
            validation_passed = False
        elif state.gene_set.size > 800:
            st.error(f"❌ **Gene input verfication failed!** Your input contains {state.gene_set.size} validated genes, but the maximum allowed is 800 genes. Please reduce your gene list size and try again.")
            validation_passed = False
        
        # Check 2: Background is valid
        if not state.background_gene_set:
            st.error("❌ **Background validation failed!** No valid background gene set loaded. Please try refreshing the page.")
            validation_passed = False
        
        # Check 3: At least 1 selected library
        if not state.gene_set_libraries or len(state.gene_set_libraries) == 0:
            st.error("❌ **Library validation failed!** No gene set libraries selected. Please select at least one library.")
            validation_passed = False
        
        # Only proceed if all validations pass
        if validation_passed and ready_common and state.gene_set_input:

            # Load background once and reuse for all libraries
            logger.info(f"Loading background gene set: {state.background_gene_set.name} ({state.background_gene_set.size} genes)")
            
            # Create progress container for iterative enrichment
            progress_container = st.container()
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            # Generate a shared run ID for all libraries in this submission
            from datetime import datetime
            shared_run_id = datetime.now().strftime("%Y%m%d_%H%M%S_%f")[:-3]  # Include milliseconds
            logger.info(f"Generated shared run ID for this submission: {shared_run_id}")
            
            total_libraries = len(state.gene_set_libraries)
            
            for i, gsl in enumerate(state.gene_set_libraries):
                # Update progress
                progress = (i + 1) / total_libraries
                progress_bar.progress(progress)
                status_text.text(f"Processing library {i+1}/{total_libraries}: {gsl.name} ({gsl.num_terms} terms)")
                
                logger.info(f"Starting iterative enrichment for library: {gsl.name}")
                
                try:
                    # Create a progress callback for this library
                    def progress_callback(message: str):
                        status_text.text(f"Processing library {i+1}/{total_libraries}: {gsl.name} - {message}")
                    
                    it = IterativeEnrichment(
                        gene_set=state.gene_set,
                        gene_set_library=gsl,
                        background_gene_set=state.background_gene_set,
                        min_term_size=state.min_term_size,
                        max_term_size=state.max_term_size,
                        p_value_method_name=state.p_val_method,
                        p_threshold=state.iter_p_threshold,
                        max_iterations=(
                            None if state.iter_max_iter == 0 else state.iter_max_iter
                        ),
                        min_overlap=state.iter_min_overlap,
                        progress_callback=progress_callback,
                        run_id=shared_run_id,
                    )
                    
                    # store enrichment object and results
                    state.iter_enrich[gsl.name] = it
                    state.iter_results[gsl.name] = it.results
                    
                    # Generate DOT from results
                    state.iter_dot[gsl.name] = it.to_dot()
                    
                    # Save main summary files to results folder (like regular mode)
                    if it.results:  # Only save if there are results
                        it.save_to_results_folder()
                    
                    logger.info(f"Completed iterative enrichment for library: {gsl.name} ({len(it.results)} iterations)")
                    
                except Exception as e:
                    logger.error(f"Error processing library {gsl.name}: {e}")
                    st.error(f"❌ Error processing library {gsl.name}: {str(e)}")
                    continue
            
            # Clear progress indicators
            progress_bar.empty()
            status_text.empty()
            
            # Reset benchmarking when new iGEA results are computed
            # This ensures benchmarking recomputes when iGEA is rerun with different parameters
            if 'benchmark_computed' in state:
                state.benchmark_computed = False
                state.benchmark_data = None
                state.benchmark_report_text = None
                state.libraries_with_data = []
                state.libraries_without_data = []
                if 'last_enrich_hash' in state:
                    del state.last_enrich_hash
            
            state.iter_ready = True

    # Iterative rendering
    if mode == "Iterative" and state.iter_ready:
        # Download section as a table
        st.markdown("**📥 Download Results:**")
        
        # Extract ORA results
        ora_results = _extract_iteration_1_as_ora(state.iter_enrich)
        combined_ora = _build_ora_tables_download(ora_results) if ora_results else None
        ora_archive_path = _create_combined_ora_archive(state.iter_enrich) if ora_results else None
        
        # Build iGEA results
        combined_igea = _build_iterative_tables_download(state.iter_results)
        combined_archive_path = _create_combined_iteration_archive({})
        
        # Create table data
        ora_combined_link = download_link(combined_ora, 'ora_enrichment_results', 'tsv') if combined_ora else "—"
        ora_individual_link = download_file_link(ora_archive_path, Path(ora_archive_path).name.replace('.tar.gz', ''), 'tar.gz') if ora_archive_path else "—"
        igea_combined_link = download_link(combined_igea, 'iterative_enrichment_results', 'tsv')
        igea_individual_link = download_file_link(combined_archive_path, Path(combined_archive_path).name.replace('.tar.gz', ''), 'tar.gz') if combined_archive_path else "—"
        
        # Display as table using markdown
        st.markdown("""
        | | ORA | iGEA |
        |---|---|---|
        | **Combined** | {ora_combined} | {igea_combined} |
        | **Individual library files** | {ora_individual} | {igea_individual} |
        """.format(
            ora_combined=ora_combined_link,
            ora_individual=ora_individual_link,
            igea_combined=igea_combined_link,
            igea_individual=igea_individual_link
        ), unsafe_allow_html=True)
        
        # Download validated gene symbols
        if hasattr(state, 'gene_set') and state.gene_set and state.gene_set.genes:
            validated_genes = sorted(list(state.gene_set.genes))
            validated_genes_text = '\n'.join(validated_genes)
            gene_set_name = state.gene_set_name if hasattr(state, 'gene_set_name') and state.gene_set_name else 'validated_genes'
            st.markdown(
                f"**📋 Validated Gene Symbols:** {download_link(validated_genes_text, gene_set_name, 'txt')}",
                unsafe_allow_html=True
            )

        # callback to keep checkbox state in session
        def toggle_library(lib_name):
            # Sanitize library name for widget key (replace special characters)
            safe_lib_name = lib_name.replace(':', '_').replace(' ', '_').replace('-', '_')
            if state[f"use_{safe_lib_name}_in_network"]:
                if lib_name not in state.selected_dot_paths:
                    state.selected_dot_paths.append(lib_name)
            else:
                if lib_name in state.selected_dot_paths:
                    state.selected_dot_paths.remove(lib_name)
            # Clear network when selection changes
            state.network_generated = False
            state.last_merged_dot = ""

        def toggle_network_selection(lib_name):
            # Sanitize library name for widget key (replace special characters)
            safe_lib_name = lib_name.replace(':', '_').replace(' ', '_').replace('-', '_')
            if state[f"network_select_{safe_lib_name}"]:
                if lib_name not in state.selected_dot_paths:
                    state.selected_dot_paths.append(lib_name)
            else:
                if lib_name in state.selected_dot_paths:
                    state.selected_dot_paths.remove(lib_name)
            # Clear network when selection changes
            state.network_generated = False
            state.last_merged_dot = ""

        # render each library's results with a persistent checkbox in the same order as the checkbox list
        def sort_libraries(libraries):
            hallmark = [lib for lib in libraries if lib.startswith("H:")]
            c2_libs = sorted([lib for lib in libraries if lib.startswith("C2:")])
            c5_libs = sorted([lib for lib in libraries if lib.startswith("C5:")])
            protein_interaction = [lib for lib in libraries if lib.startswith("Protein Interaction")]
            other = sorted([lib for lib in libraries if not any(lib.startswith(prefix) for prefix in ["H:", "C2:", "C5:", "Protein Interaction"])])
            
            return hallmark + c2_libs + c5_libs + protein_interaction + other
        
        available_libraries = sort_libraries(list(state.iter_enrich.keys()))
        for lib in available_libraries:
            it = state.iter_enrich[lib]
            render_iter_results(it, lib)
            # Sanitize library name for widget key (replace special characters)
            safe_lib_name = lib.replace(':', '_').replace(' ', '_').replace('-', '_')
            st.checkbox(
                "Use results in network",
                value=lib in state.selected_dot_paths,
                key=f"use_{safe_lib_name}_in_network",
                on_change=toggle_library,
                args=(lib,)
            )
        
        # Display comparison statistics (ORA = iteration 1, iGEA = all iterations)
        # No need for separate regular enrichment - iteration 1 IS the ORA
        logger.info(f"Checking comparison display: iter_enrich={bool(state.iter_enrich)}, len(iter_enrich)={len(state.iter_enrich) if state.iter_enrich else 0}")
        if state.iter_enrich and len(state.iter_enrich) > 0:
            # Check if we have at least iteration 1
            has_iteration_1 = any(
                iter_enrich and iter_enrich.results and 
                any(r.get("Iteration") == 1 for r in iter_enrich.results)
                for iter_enrich in state.iter_enrich.values()
            )
            
            if has_iteration_1:
                logger.info("Displaying ORA vs iGEA comparison statistics (ORA = iteration 1)")
                st.markdown("---")
                st.markdown("## 📊 ORA vs iGEA Comparison")
                st.caption("ORA = Iteration 1 | iGEA = All Iterations")
                try:
                    # Pass empty dict for regular_enrichments since we use iteration 1 as ORA
                    render_ora_igea_comparison(
                        {},  # No separate regular enrichment needed
                        state.iter_enrich,
                        p_threshold=state.iter_p_threshold,
                        fdr_threshold=state.adjusted_p_threshold,
                    )
                except Exception as e:
                    logger.error(f"Error displaying comparison: {e}")
                    st.error(f"Error displaying comparison statistics: {e}")
                    import traceback
                    st.code(traceback.format_exc())
            else:
                logger.info("Comparison not displayed - no iteration 1 found")
        else:
            logger.info("Comparison not displayed - no iterative results")

        # Network section
        st.markdown("---")
        st.header("Generate network for conversational AI systems")
        
        # Interactive library selection interface
        st.subheader("Library Selection for Network")
        
        # Get all available libraries in specific order: Hallmark, C2 (alphabetical), C5 (alphabetical), Protein Interaction
        def sort_libraries(libraries):
            hallmark = [lib for lib in libraries if lib.startswith("H:")]
            c2_libs = sorted([lib for lib in libraries if lib.startswith("C2:")])
            c5_libs = sorted([lib for lib in libraries if lib.startswith("C5:")])
            protein_interaction = [lib for lib in libraries if lib.startswith("Protein Interaction")]
            other = sorted([lib for lib in libraries if not any(lib.startswith(prefix) for prefix in ["H:", "C2:", "C5:", "Protein Interaction"])])
            
            return hallmark + c2_libs + c5_libs + protein_interaction + other
        
        available_libraries = sort_libraries(list(state.iter_enrich.keys()))
        
        # Create columns for better layout
        col1, col2 = st.columns([3, 1])
        
        with col1:
            # Interactive list with checkboxes for each library
            st.write("**Available Libraries:**")
            for lib in available_libraries:
                # Check if this library has any enrichment results
                has_results = not state.iter_enrich[lib].to_dataframe().empty
                
                # Sanitize library name for widget key (replace special characters)
                safe_lib_name = lib.replace(':', '_').replace(' ', '_').replace('-', '_')
                
                if has_results:
                    st.checkbox(
                        lib,
                        value=lib in state.selected_dot_paths,
                        key=f"network_select_{safe_lib_name}",
                        on_change=toggle_network_selection,
                        args=(lib,)
                    )
                else:
                    # Gray out libraries with no results
                    st.markdown(f"~~{lib}~~ *(no enrichment results)*")
        
        with col2:
            # Quick selection controls
            st.write("**Quick Actions:**")
            
            if st.button("Select All"):
                # Only select libraries that have results
                libraries_with_results = [lib for lib in available_libraries if not state.iter_enrich[lib].to_dataframe().empty]
                state.selected_dot_paths = libraries_with_results.copy()
                # Clear network when selection changes
                state.network_generated = False
                state.last_merged_dot = ""
                st.rerun()
            
            if st.button("Clear All"):
                state.selected_dot_paths = []
                # Clear network when selection changes
                state.network_generated = False
                state.last_merged_dot = ""
                st.rerun()

        # generate or re-display merged network
        if st.button("Generate Network"):
            if not state.selected_dot_paths:
                st.error("Please select at least one library for network generation.")
            else:
                state.network_generated = True
                available = set(state.iter_dot.keys())
                chosen = [lib for lib in state.selected_dot_paths if lib in available]
                if not chosen:
                    st.error("No valid libraries selected for network generation.")
                else:
                    selected_dots = {lib: state.iter_dot[lib] for lib in chosen}
                    state.last_merged_dot = merge_iterative_dot(selected_dots)
                    # Get gene count and library list for AI prompt
                    gene_count = state.gene_set.size if hasattr(state, 'gene_set') and state.gene_set else None
                    library_list = chosen if chosen else None
                    render_network(state.last_merged_dot, gene_count=gene_count, library_list=library_list)
        elif state.network_generated and state.selected_dot_paths:
            # Get gene count and library list for AI prompt
            gene_count = state.gene_set.size if hasattr(state, 'gene_set') and state.gene_set else None
            library_list = state.selected_dot_paths if state.selected_dot_paths else None
            render_network(state.last_merged_dot, gene_count=gene_count, library_list=library_list)
        
        # Statistical Benchmarking Section
        st.markdown("---")
        st.header("📊 Statistical Benchmarking")
        st.caption("Compare your network connectivity against random gene lists of similar size")
        
        # Check if we have results to benchmark
        if state.iter_enrich and len(state.iter_enrich) > 0:
            # Check if benchmarking needs to be computed or recomputed
            gene_list_size = state.gene_set.size if state.gene_set else 0
            p_threshold = state.iter_p_threshold
            
            # Get a hash of the current iter_enrich to detect if it changed
            # This ensures benchmarking recomputes when iGEA is rerun
            current_enrich_hash = hash(tuple(sorted(state.iter_enrich.keys()))) if state.iter_enrich else 0
            
            # Check if we need to recompute (e.g., if gene list size, p-threshold, or iGEA results changed)
            needs_computation = (
                'benchmark_computed' not in state or
                not state.benchmark_computed or
                state.get('last_benchmark_gene_list_size') != gene_list_size or
                state.get('last_benchmark_p_threshold') != p_threshold or
                state.get('last_enrich_hash') != current_enrich_hash
            )
            
            if needs_computation:
                with st.spinner("Computing statistical benchmark... This may take a minute."):
                    try:
                        # Set up paths
                        parquet_dir = ROOT / "permutations" / "permutation_cluster_statistics_parquet"
                        
                        # Get user's max_iterations (cap at 30 for permutation data)
                        user_max_iter = state.iter_max_iter if state.iter_max_iter > 0 else None
                        if user_max_iter is not None:
                            user_max_iter = min(user_max_iter, 30)  # Cap at 30 (permutation data max)
                        
                        # Compute benchmark
                        null_dist, cluster_benchmarks, libs_with_data, libs_without_data, actual_size_used, analyzer = compute_benchmark_for_streamlit(
                            state.iter_enrich,
                            gene_list_size,
                            p_threshold,
                            parquet_dir,
                            user_max_iter
                        )
                        
                        if null_dist and cluster_benchmarks:
                            state.benchmark_computed = True
                            state.benchmark_data = {
                                'null_distribution': null_dist,
                                'cluster_benchmarks': cluster_benchmarks
                            }
                            state.libraries_with_data = libs_with_data
                            state.libraries_without_data = libs_without_data
                            state.actual_size_used = actual_size_used
                            state.last_benchmark_gene_list_size = gene_list_size
                            state.last_benchmark_p_threshold = p_threshold
                            state.last_enrich_hash = current_enrich_hash
                            state.benchmark_analyzer = analyzer  # Store analyzer for report generation
                            
                            # Generate report text
                            gene_list_name = state.gene_set_name if hasattr(state, 'gene_set_name') and state.gene_set_name else "Gene List"
                            state.benchmark_report_text = generate_statistical_report_text(
                                cluster_benchmarks,
                                gene_list_name,
                                libs_with_data,
                                libs_without_data,
                                analyzer,
                                gene_list_size=gene_list_size,
                                actual_size_used=actual_size_used
                            )
                        else:
                            if not libs_with_data:
                                st.warning("⚠️ No libraries with permutation data available for benchmarking. Please ensure Parquet files exist in the results directory.")
                            elif null_dist is None:
                                st.warning("⚠️ Failed to compute null distribution. This may occur if:\n"
                                          "- P-value threshold is too strict (> 0.05)\n"
                                          "- No permutation data matches your gene list size\n"
                                          "- Permutation data files are missing or corrupted")
                            elif cluster_benchmarks is None or len(cluster_benchmarks) == 0:
                                st.warning("⚠️ No clusters found in your network. Benchmarking requires at least one cluster.")
                            else:
                                st.warning("⚠️ Failed to compute benchmark. Please check the logs for details.")
                    except Exception as e:
                        logger.error(f"Error computing benchmark: {e}", exc_info=True)
                        import traceback
                        error_details = traceback.format_exc()
                        logger.error(f"Full traceback: {error_details}")
                        st.error(f"❌ Error computing benchmark: {str(e)}\n\nCheck the terminal/console for detailed error logs.")
            
            # Display benchmark results if computed
            if state.benchmark_computed and state.benchmark_data:
                cluster_benchmarks = state.benchmark_data['cluster_benchmarks']
                
                if cluster_benchmarks:
                    # Gene list size information
                    actual_size = state.get('actual_size_used', gene_list_size)
                    size_info = f"**Gene list size:** {gene_list_size} genes"
                    if actual_size != gene_list_size:
                        size_info += f" (using permutation data for size {actual_size})"
                    else:
                        size_info += f" (exact match with permutation data)"
                    st.info(size_info)
                    
                    # Library information - improved message
                    if state.libraries_with_data:
                        libs_msg = f"**Libraries used for statistical benchmarking (based on permutation data availability):** {', '.join(state.libraries_with_data)} ({len(state.libraries_with_data)} librar"
                        if len(state.libraries_with_data) != 1:
                            libs_msg += "ies"
                        else:
                            libs_msg += "y"
                        libs_msg += ")"
                        st.info(libs_msg)
                    
                    if state.libraries_without_data:
                        # Check if excluded libraries actually have data in the cluster
                        excluded_msg = f"**Note:** The following libraries were included in your enrichment analysis but excluded from statistical benchmarking (no permutation data available): {', '.join(state.libraries_without_data)}"
                        st.warning(excluded_msg)
                    
                    # Add comprehensive note about filtering and parameters
                    filter_notes = []
                    if state.last_benchmark_p_threshold is not None:
                        filter_notes.append(f"p-value ≤ {state.last_benchmark_p_threshold}")
                    
                    user_max_iter = state.iter_max_iter if state.iter_max_iter > 0 else None
                    if user_max_iter is not None:
                        max_iter_display = min(user_max_iter, 30)  # Cap at 30
                        filter_notes.append(f"max iterations ≤ {max_iter_display}")
                    else:
                        filter_notes.append("max iterations ≤ 30 (permutation data maximum)")
                    
                    # Check if user parameters differ from permutation defaults
                    permutation_defaults = {
                        'min_overlap': 3,
                        'min_term_size': 10,
                        'max_term_size': 600
                    }
                    
                    user_params = {
                        'min_overlap': state.iter_min_overlap,
                        'min_term_size': state.iter_min_term_size,
                        'max_term_size': state.iter_max_term_size
                    }
                    
                    param_warnings = []
                    if user_params['min_overlap'] != permutation_defaults['min_overlap']:
                        param_warnings.append(f"Minimum overlap: {user_params['min_overlap']} (permutation default: {permutation_defaults['min_overlap']})")
                    if user_params['min_term_size'] != permutation_defaults['min_term_size']:
                        param_warnings.append(f"Minimum term size: {user_params['min_term_size']} (permutation default: {permutation_defaults['min_term_size']})")
                    if user_params['max_term_size'] != permutation_defaults['max_term_size']:
                        param_warnings.append(f"Maximum term size: {user_params['max_term_size']} (permutation default: {permutation_defaults['max_term_size']})")
                    
                    # Build comprehensive info message
                    info_parts = []
                    if filter_notes:
                        info_parts.append(f"**Statistical filtering:** Permutation data was filtered to match your analysis parameters: {', '.join(filter_notes)}.")
                    
                    if param_warnings:
                        info_parts.append("\n\n⚠️ **Parameter mismatch:** Your analysis uses different parameters than the permutation data:")
                        info_parts.append("\n".join(f"- {w}" for w in param_warnings))
                        info_parts.append(
                            f"\n\nFor more accurate statistical analysis, consider using the default parameters "
                            f"(min overlap: {permutation_defaults['min_overlap']}, "
                            f"min term size: {permutation_defaults['min_term_size']}, "
                            f"max term size: {permutation_defaults['max_term_size']}) "
                            "that were used to generate the permutation data."
                        )
                    
                    if info_parts:
                        info_msg = "".join(info_parts)
                        st.info(info_msg)
                    
                    # Extract table data for largest cluster
                    table_data = extract_benchmark_table_data(cluster_benchmarks)
                    
                    if table_data:
                        st.subheader("Statistical Benchmarks vs Random Gene Lists")
                        st.caption("Comparison of largest cluster metrics against null distribution from random gene lists")
                        
                        # Create DataFrame for display
                        df = pd.DataFrame(table_data)
                        
                        # Display table
                        st.dataframe(
                            df,
                            use_container_width=True,
                            hide_index=True
                        )
                        
                        # Interpretation guide
                        with st.expander("📖 How to interpret these results"):
                            st.markdown("""
                            **Z-score and Percentile:**
                            - **Z-score > 2.0 and Percentile > 95%**: ✓✓ SIGNIFICANTLY BETTER (top 5%)
                            - **Z-score > 1.0 and Percentile > 84%**: ✓ BETTER (top 16%)
                            - **Z-score ~ 0.0 and Percentile ~ 50%**: ~ SIMILAR to random
                            - **Z-score < -1.0 and Percentile < 16%**: ✗ WORSE than random
                            
                            **What this means:**
                            - Higher values (genes, terms, edges, density) indicate better network connectivity
                            - Your gene list is compared against thousands of random gene lists of similar size
                            - Significant results suggest biological coherence in your gene list
                            """)
                        
                        # Download link for full report
                        if state.benchmark_report_text:
                            st.markdown("---")
                            st.markdown("**📄 Download Full Statistical Report:**")
                            st.download_button(
                                label="Download Statistical Report (Text)",
                                data=state.benchmark_report_text,
                                file_name=f"{state.gene_set_name if hasattr(state, 'gene_set_name') and state.gene_set_name else 'statistical_report'}.txt",
                                mime="text/plain"
                            )
                else:
                    st.warning("No clusters found in the network. Benchmarking requires at least one cluster.")
            elif state.benchmark_computed:
                st.warning("Benchmark computation completed but no cluster data available.")
        else:
            st.info("Run iterative enrichment analysis first to enable statistical benchmarking.")

    logger.info("Finishing the Streamlit app")


if __name__ == "__main__":
    main()
