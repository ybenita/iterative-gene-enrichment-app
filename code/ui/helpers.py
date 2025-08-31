import logging
from pathlib import Path

import streamlit as st
from streamlit import session_state as state

from gene_converter import GeneConverter

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)
ROOT = Path(__file__).resolve().parent.parent.parent


def input_example() -> None:
    """Input an example gene set from the example_gene_list.txt file."""
    if not hasattr(state, 'gene_input_format'):
        state.gene_input_format = 'symbols'
    
    # Read from the example_gene_list.txt file
    example_file_path = ROOT / "data" / "gene_lists" / "example_gene_list.txt"
    
    try:
        with open(example_file_path, "r") as f:
            example_genes = f.read().strip()
        
        state.gene_set_input = example_genes
        state.gene_set_name = "Example gene set"
        
        logger.info(f"Loaded example gene set from {example_file_path}")
        
    except FileNotFoundError:
        logger.error(f"Example gene list file not found: {example_file_path}")
        # Fallback to hardcoded examples if file doesn't exist
        if state.gene_input_format == 'entrez_ids':
            example_genes = """7157
672
675
7159
7160
7158
7161
7162
7163
7164"""
            state.gene_set_name = "Example gene set (Entrez IDs)"
        else:
            example_genes = """ABHD11
ABHD14A
ABHD3
ACAA1A
ACBD4
ACO1
ADH5
ADHFE1
ADK
AFAP1L1"""
            state.gene_set_name = "Example gene set (Symbols)"
        
        state.gene_set_input = example_genes
        st.warning("⚠️ Example gene list file not found. Using fallback examples.")


def update_text_widgets() -> None:
    """Update text widgets when file selection changes."""
    if state.selected_file != "Select ...":
        file_path = ROOT / "data" / "gene_lists" / state.selected_file
        try:
            with open(file_path, "r") as f:
                content = f.read()
                state.gene_set_input = content
                state.gene_set_name = state.selected_file.replace(".txt", "")
        except FileNotFoundError:
            st.error(f"File {file_path} not found")


def convert_and_validate_gene_input(input_text: str, input_format: str) -> tuple[list[str], list[str], list[str], dict, list[str]]:
    """
    Convert gene input based on specified format to symbols and provide validation feedback.
    
    Args:
        input_text: Raw input text containing gene identifiers
        input_format: Either 'symbols' or 'entrez_ids'
        
    Returns:
        Tuple of (converted_symbols, unrecognized_entrez, unrecognized_symbols, stats, conversions)
    """
    if not hasattr(state, 'gene_converter'):
        state.gene_converter = GeneConverter()
    
    # Clear previous conversions
    state.gene_converter.clear_conversions()
    
    # Parse input lines
    lines = [line.strip() for line in input_text.split('\n') if line.strip()]
    
    converted_symbols = []
    unrecognized_entrez = []
    unrecognized_symbols = []
    
    if input_format == 'entrez_ids':
        # Convert Entrez IDs to symbols
        for line in lines:
            gene_id = line.strip()
            if not gene_id:
                continue
            
            symbol = state.gene_converter.get_symbol(gene_id)
            if symbol:
                converted_symbols.append(symbol)
            else:
                unrecognized_entrez.append(gene_id)
    else:
        # Validate symbols (input_format == 'symbols')
        for line in lines:
            gene_id = line.strip()
            if not gene_id:
                continue
            
            # Use validate_and_map_symbol to get the mapped symbol
            mapped_symbol = state.gene_converter.validate_and_map_symbol(gene_id)
            if mapped_symbol:
                converted_symbols.append(mapped_symbol)
            else:
                unrecognized_symbols.append(gene_id)
    
    # Get conversion statistics and conversions
    stats = state.gene_converter.get_stats()
    conversions = state.gene_converter.get_conversions()
    
    return converted_symbols, unrecognized_entrez, unrecognized_symbols, stats, conversions


def display_conversion_results(converted_symbols: list[str], 
                              unrecognized_entrez: list[str], 
                              unrecognized_symbols: list[str],
                              stats: dict,
                              input_format: str,
                              conversions: list[str] = None) -> None:
    """
    Display the results of gene input conversion and validation.
    
    Args:
        converted_symbols: List of successfully converted gene symbols
        unrecognized_entrez: List of unrecognized Entrez IDs
        unrecognized_symbols: List of unrecognized gene symbols
        stats: Dictionary with conversion statistics
        input_format: The input format that was used ('symbols' or 'entrez_ids')
        conversions: List of gene symbol conversions (e.g., ["P53→TP53", "P21→CDKN1A"])
    """
    total_input = len(converted_symbols) + len(unrecognized_entrez) + len(unrecognized_symbols)
    
    if total_input == 0:
        return
    
    # Display conversion summary
    with st.expander("Gene Input Validation Summary", expanded=True):
        col1, col2 = st.columns(2)
        
        with col1:
            st.metric("Valid Genes", len(converted_symbols))
        
        with col2:
            if input_format == 'entrez_ids':
                st.metric("Invalid Entrez IDs", len(unrecognized_entrez))
            else:
                st.metric("Invalid Symbols", len(unrecognized_symbols))
        
        # Show validation details
        if converted_symbols:
            st.success(f"✅ Successfully processed {len(converted_symbols)} genes")
        
        # Show conversions if any
        if conversions and len(conversions) > 0:
            st.info(f"🔄 {len(conversions)} gene symbols converted")
            st.caption(", ".join(conversions))
        
        if input_format == 'entrez_ids':
            if unrecognized_entrez:
                st.warning(f"⚠️ {len(unrecognized_entrez)} Entrez IDs not found in database")
                st.caption("Invalid Entrez IDs: " + ", ".join(unrecognized_entrez))
        else:
            if unrecognized_symbols:
                st.warning(f"⚠️ {len(unrecognized_symbols)} gene symbols not found in database")
                st.caption("Invalid symbols: " + ", ".join(unrecognized_symbols))
        

