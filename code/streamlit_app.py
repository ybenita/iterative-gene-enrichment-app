import json
import logging
import math
import re
from io import StringIO
from pathlib import Path
from typing import Dict, List, Set
from datetime import datetime

import streamlit as st
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
)
from ui.utils import download_link, download_file_link, update_aliases

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
    rows = ["Library\tIteration\tTerm\tDescription\tOverlap size\tp-value\t-log(p-value)\tGenes"]
    for lib, records in all_iter_results.items():
        for rec in records:
            p = rec.get("p-value", float("nan"))
            overlap_size = rec.get("Overlap size", "0/0")
            description = rec.get("Description", "")
            genes = ', '.join(rec.get('Genes', []))
            log_p = -math.log10(p) if p and p > 0 else 0
            rows.append(
                f"{lib}\t{rec['Iteration']}\t{rec['Term']}\t{description}\t{overlap_size}\t{p}\t{log_p}\t{genes}"
            )
    return "\n".join(rows)


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


def main() -> None:
    logger.info("Starting the Streamlit app")
    st.sidebar.image(
        Image.open(ROOT / "code" / "static" / "logo.png"),
        caption="Iterative Enrichment Analysis",
    )
    st.sidebar.title("Enrichment analysis")
    st.sidebar.write(
        """This app tests the input gene list against selected pathway libraries. 
        
**Regular Mode**: Runs once and reports all results that pass the filters per library.

**Iterative Mode**: Iteratively removes genes from the top hit at each step until no term passes the p-value cutoff.

Results include ranked tables, bar charts, and network graphs."""
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
                "Enter gene symbols (e.g., TP53, BRCA1) - one per line (max 500 genes)" 
                if state.gene_input_format == 'symbols' 
                else "Enter Entrez IDs (e.g., 7157, 672) - one per line (max 500 genes)"
            )
            
            st.text_area(
                "Input a set of genes",
                key="gene_set_input",
                height=400,
                placeholder=placeholder_text,
                label_visibility="collapsed",
            )
            st.caption("📝 **Note:** Maximum 500 genes allowed for optimal performance")
            
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
                state.gene_set_libraries = [
                    GeneSetLibrary(
                        str(ROOT / "data" / "libraries" / state.lib_mapper[lib]),
                        name=lib,
                    )
                    for lib in state.libraries
                ]
            else:
                state.gene_set_libraries = []
            # Always load background gene set
            from ui.utils import get_background_info
            bg_file_path, bg_format = get_background_info(state.background_set)
            
            if bg_file_path:
                state.background_gene_set = BackgroundGeneSet(
                    bg_file_path,
                    input_format=bg_format  # Use format from alias file
                )
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
        elif state.gene_set.size > 500:
            st.error(f"❌ **Gene validation failed!** Your input contains {state.gene_set.size} validated genes, but the maximum allowed is 500 genes. Please reduce your gene list size and try again.")
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
        

        
        # Always render results section
        st.markdown(
            f"Download all results as {download_link(collect_results(state.enrich), 'regular_enrichment_results','tsv')}, "
            f"{download_link(generate_regular_enrichment_json_analysis(state.enrich), 'regular_enrichment_results','json')}",
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
        
        # Initialize network selection state if not exists
        if not hasattr(state, 'selected_regular_libraries'):
            state.selected_regular_libraries = []
        if not hasattr(state, 'regular_network_generated'):
            state.regular_network_generated = False
        if not hasattr(state, 'last_regular_network'):
            state.last_regular_network = ""
        
        # AI prompt construction section for regular enrichment
        st.markdown("---")
        st.header("Construction of AI Prompt")
        
        # callback to keep checkbox state in session
        def toggle_regular_library(lib_name):
            # Sanitize library name for widget key (replace special characters)
            safe_lib_name = lib_name.replace(':', '_').replace(' ', '_').replace('-', '_')
            if state[f"regular_network_select_{safe_lib_name}"]:
                if lib_name not in state.selected_regular_libraries:
                    state.selected_regular_libraries.append(lib_name)
            else:
                if lib_name in state.selected_regular_libraries:
                    state.selected_regular_libraries.remove(lib_name)
            # Clear network when selection changes
            state.regular_network_generated = False
            state.last_regular_network = ""
        
        # Interactive library selection interface
        st.subheader("Library Selection for AI Prompt")
        
        # Get all available libraries in specific order: Hallmark, C2 (alphabetical), C5 (alphabetical), Protein Interaction
        def sort_libraries(libraries):
            hallmark = [lib for lib in libraries if lib.startswith("H:")]
            c2_libs = sorted([lib for lib in libraries if lib.startswith("C2:")])
            c5_libs = sorted([lib for lib in libraries if lib.startswith("C5:")])
            protein_interaction = [lib for lib in libraries if lib.startswith("Protein Interaction")]
            other = sorted([lib for lib in libraries if not any(lib.startswith(prefix) for prefix in ["H:", "C2:", "C5:", "Protein Interaction"])])
            
            return hallmark + c2_libs + c5_libs + protein_interaction + other
        
        available_libraries = sort_libraries(list(state.enrich.keys()))
        
        # Create columns for better layout
        col1, col2 = st.columns([3, 1])
        
        with col1:
            # Interactive list with checkboxes for each library
            st.write("**Available Libraries:**")
            for lib in available_libraries:
                # Check if this library has any enrichment results
                has_results = not state.enrich[lib].to_dataframe().empty
                
                # Sanitize library name for widget key (replace special characters)
                safe_lib_name = lib.replace(':', '_').replace(' ', '_').replace('-', '_')
                
                if has_results:
                    st.checkbox(
                        lib,
                        value=lib in state.selected_regular_libraries,
                        key=f"regular_network_select_{safe_lib_name}",
                        on_change=toggle_regular_library,
                        args=(lib,)
                    )
                else:
                    # Gray out libraries with no results
                    st.markdown(f"~~{lib}~~ *(no enrichment results)*")
        
        with col2:
            # Quick selection controls
            st.write("**Quick Actions:**")
            
            if st.button("Select All", key="regular_select_all"):
                # Only select libraries that have results
                libraries_with_results = [lib for lib in available_libraries if not state.enrich[lib].to_dataframe().empty]
                state.selected_regular_libraries = libraries_with_results.copy()
                # Clear network when selection changes
                state.regular_network_generated = False
                state.last_regular_network = ""
                st.rerun()
            
            if st.button("Clear All", key="regular_clear_all"):
                state.selected_regular_libraries = []
                # Clear network when selection changes
                state.regular_network_generated = False
                state.last_regular_network = ""
                st.rerun()

        # Generate AI prompt
        if st.button("Generate AI Prompt", key="regular_generate_network"):
            if not state.selected_regular_libraries:
                st.error("Please select at least one library for AI prompt generation.")
            else:
                state.regular_network_generated = True
                # Generate JSON network from regular enrichment results
                network_json = generate_regular_network_json(state.enrich, state.selected_regular_libraries)
                state.last_regular_network = network_json
                render_regular_network_analysis(network_json, state.selected_regular_libraries)
        elif state.regular_network_generated and state.selected_regular_libraries:
            render_regular_network_analysis(state.last_regular_network, state.selected_regular_libraries)

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
        elif state.gene_set.size > 500:
            st.error(f"❌ **Gene input verfication failed!** Your input contains {state.gene_set.size} validated genes, but the maximum allowed is 500 genes. Please reduce your gene list size and try again.")
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
            

            
            state.iter_ready = True

    # Iterative rendering
    if mode == "Iterative" and state.iter_ready:
        # Download section with just two links
        st.markdown("**📥 Download Results:**")
        
        # 1. Combined iterative enrichment TSV
        combined = _build_iterative_tables_download(state.iter_results)
        st.markdown(
            f"📊 **Combined Results:** {download_link(combined, 'iterative_enrichment_results', 'tsv')}",
            unsafe_allow_html=True,
        )
        
        # 2. All individual iteration files as single archive
        combined_archive_path = _create_combined_iteration_archive({})
        if combined_archive_path:
            combined_archive_name = Path(combined_archive_path).name
            # Remove the .tar.gz extension since download_file_link will add it back
            base_filename = combined_archive_name.replace('.tar.gz', '')
            st.markdown(
                f"📁 **All Individual Files:** {download_file_link(combined_archive_path, base_filename, 'tar.gz')}",
                unsafe_allow_html=True,
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

        # Network section
        st.markdown("---")
        st.header("Network")
        
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
                    render_network(state.last_merged_dot)
        elif state.network_generated and state.selected_dot_paths:
            render_network(state.last_merged_dot)

    logger.info("Finishing the Streamlit app")


def generate_regular_network_json(enrich_results, selected_libraries):
    """
    Generate a lightweight JSON representation from regular enrichment results.
    This provides only essential data: Library, Term, Rank, Genes.
    
    :param enrich_results: Dictionary of enrichment results by library
    :param selected_libraries: List of selected library names
    :return: JSON string with lightweight format
    """
    import json
    
    # Extract data from enrichment results for selected libraries only
    results = []
    
    for lib_name in selected_libraries:
        if lib_name in enrich_results:
            enrich = enrich_results[lib_name]
            
            # Get the DataFrame with all results
            df = enrich.to_dataframe()
            
            # Extract only the required columns: Library, Term, Rank, Genes
            for _, row in df.iterrows():
                # Convert genes string back to list
                genes = [gene.strip() for gene in row['Genes'].split(',') if gene.strip()]
                
                results.append({
                    "library": row['Library'],
                    "term": row['Term'],
                    "rank": int(row['Rank']),
                    "genes": genes
                })
    
    # Create the lightweight JSON structure
    analysis_data = {
        "analysis_type": "over_representation_analysis",
        "results": results
    }
    
    return json.dumps(analysis_data, indent=2)


def render_regular_network_analysis(network_json, selected_libraries):
    """
    Render the network analysis for regular enrichment results.
    
    :param network_json: JSON network content
    :param selected_libraries: List of selected library names
    """
    st.subheader("AI Analysis")
    
    # Generate AI prompt for regular enrichment
    ai_prompt = generate_regular_ai_prompt(network_json, selected_libraries)
    
    # Display the prompt in an expander
    with st.expander("📋 AI Analysis Prompt", expanded=False):
        st.text_area(
            "Copy this prompt to use with AI analysis tools:",
            value=ai_prompt,
            height=400,
            key="regular_ai_prompt"
        )
    
    # Offer download
    st.markdown(
        f'Download {download_link(ai_prompt, "regular_ai_analysis_prompt", "txt")} for AI analysis',
        unsafe_allow_html=True,
    )


def generate_regular_ai_prompt(network_json, selected_libraries):
    """
    Generate AI analysis prompt for regular enrichment results.
    
    :param network_json: JSON enrichment data
    :param selected_libraries: List of selected library names
    :return: Formatted AI analysis prompt
    """
    import json
    
    # Parse the JSON data
    analysis_data = json.loads(network_json)
    
    prompt = f"""# OVER-REPRESENTATION ANALYSIS (ORA)
# Ranked List of Enriched Biological Terms

## ANALYSIS CONTEXT
You are a computational biologist analyzing over-representation analysis (ORA) results. This analysis identifies biological pathways, processes, and functions that are significantly enriched in a given gene set compared to a background gene set. The results are provided as a ranked list of enriched terms, where each term contains a list of genes that are over-represented in that biological process.

**DATA STRUCTURE:**
- **Ranked List**: Terms are ranked by statistical significance (lower rank = more significant)
- **Gene Lists**: Each term contains the specific genes that are enriched in that biological process
- **Library Sources**: Different databases provide complementary biological context
- **Table Format**: Results are organized as a table rather than a network structure

## ENRICHMENT DATA (JSON Format)
```json
{network_json}
```

## LIBRARY SOURCES AND THEIR BIOLOGICAL CONTEXT:
"""
    
    # Add library descriptions only for selected libraries
    library_descriptions = {
        "H: Hallmark Gene Sets": "Curated gene sets representing well-defined biological states or processes",
        "C2: BioCarta": "Canonical pathways from BioCarta database",
        "C2: KEGG MEDICUS": "Metabolic and signaling pathways from KEGG",
        "C2: Pathway Interaction Database": "Curated human signaling pathways",
        "C2: Reactome Pathways": "Expert-curated biological pathways",
        "C2: WikiPathways": "Community-curated biological pathways",
        "C5: Gene Ontology: Biological Process": "Biological processes from Gene Ontology",
        "C5: Gene Ontology: Cellular Component": "Cellular components from Gene Ontology",
        "C5: Gene Ontology: Molecular Function": "Molecular functions from Gene Ontology",
        "C5: Human Phenotype Ontology": "Human phenotypes and diseases",
        "Protein Interaction": "Protein-protein interaction networks"
    }
    
    for lib in selected_libraries:
        description = library_descriptions.get(lib, "Biological pathway/process database")
        prompt += f"- **{lib}**: {description}\n"
    
    prompt += """
## ANALYSIS REQUEST:

Please analyze this ranked list of enriched terms and provide:

1. **Key Biological Insights:**
   - What are the most significant biological processes/pathways identified?
   - Which genes appear in multiple terms and may be central to the biological response?

2. **Ranked List Analysis:**
   - How do the top-ranked terms relate to each other biologically?
   - Are there patterns in the gene composition across different terms?
   - What does the ranking suggest about the hierarchy of biological processes?

3. **Biological Hypothesis:**
   - Based on the ranked list structure, what biological hypothesis can you generate?
   - How do the different library sources support or complement each other?
   - Try to formulate one coherent hypothesis that explains as many terms as possible and generates a coherent biological story.
   
4. **Estimated Experimental Context**
   - Can you hypothesize on what was the experiment that generated this result?

## RESPONSE STRUCTURE:
- **Executive Summary** (2-3 sentences)
- **Key Biological Processes Identified** (by library source and rank)
- **Gene-Term Relationship Analysis**
- **Proposed Biological Hypothesis** (coherent story explaining multiple terms)
- **Estimated Experimental Context**

## IMPORTANT NOTES:
- Focus on biological interpretation, not just statistical significance
- Consider the functional relationships between genes within and across terms
- Look for genes that appear in multiple terms as potential key regulators
- Consider the broader biological context and literature
- Different library sources may provide complementary biological insights
- Build a coherent story that explains as many terms as possible rather than treating them in isolation
- The goal is to find a unifying biological hypothesis that connects the enriched terms into a meaningful narrative"""
    
    return prompt


if __name__ == "__main__":
    main()
