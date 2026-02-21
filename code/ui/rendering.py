import logging
from math import log10
from typing import Dict, List, Set

import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st
from streamlit import session_state as state

from enrichment import Enrichment
from iter_enrichment import IterativeEnrichment
from ui.dot_utils import dot_to_plotly
from ui.utils import download_link

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def render_table(result: pd.DataFrame) -> None:
    """
    Render a styled DataFrame within the Streamlit app.

    This function takes a DataFrame containing results data and applies custom formatting before
    displaying it within the Streamlit app using the `st.dataframe` function. Custom formatting includes
    specific number formatting for 'p-value' and 'fdr' columns.

    :param result: The DataFrame containing results data to display.
    """

    logger.info("Rendering DataFrame in Streamlit app.")

    def custom_format(n):
        if n > 0.001:
            return f"{n:.3f}"
        else:
            return f"{n:.3e}"

    # Create a copy of the dataframe to avoid modifying the original
    display_df = result.copy()
    
    # Convert underscores to spaces in term names for better readability
    if 'Term' in display_df.columns and not display_df.empty:
        # Ensure term column contains strings and handle NaN/None values
        display_df['Term'] = display_df['Term'].fillna('').astype(str).str.replace('_', ' ', regex=False)

    st.dataframe(
        display_df.style.format({"p-value": custom_format, "fdr": custom_format}),
        use_container_width=True,
        column_config={
            "Rank": None,
            "Term": "Term",
            "Overlap size": "Overlap size",
            "Genes": "Overlap (click to expand)",
            "Description": None,
            "Library": None,
            "p-value": "P-value",
            "FDR": "FDR",
        },
    )


def render_barchart(result: pd.DataFrame, file_name: str = "") -> None:
    """
    Render a bar chart visualization of the results data within the Streamlit app.

    This function takes a DataFrame containing results data, specifically terms and p-values, and creates
    a bar chart using Plotly Express, which is then displayed in the Streamlit app using the `st.plotly_chart` function.

    :param result: The DataFrame containing results data to visualize.
    :param file_name: Optional file name to create unique chart keys.
    """
    logger.info("Rendering bar chart in Streamlit app.")
    bar = result[["Term", "p-value"]].copy()
    
    # Convert underscores to spaces in term names for better readability
    if not bar.empty:
        bar['Term'] = bar['Term'].fillna('').astype(str).str.replace('_', ' ', regex=False)
    
    bar.loc[:, "p-value"] = bar.loc[:, "p-value"].apply(lambda x: -1 * log10(x))
    bar = bar.sort_values(by=["p-value"])
    bar.columns = ["Term", "-log10(p-value)"]
    fig = px.bar(
        bar,
        x="-log10(p-value)",
        y="Term",
        orientation="h",
        labels={"-log10(p-value)": "−log₁₀(p‐value)", "Term": "Term"},
    )
    st.plotly_chart(fig, use_container_width=True, key=f"barchart_{file_name}")


def render_results(result: Enrichment, file_name: str, n_results: int = 10) -> None:
    """
    Render a results section within the Streamlit app.

    This function processes and visualizes the result data within the Streamlit app. It provides
    a table visualization and a bar chart visualization of the top results. Additionally, it
    allows the user to download the results in various formats.

    :param result: The DataFrame containing results data to display.
    :param file_name: The name of the file to be used for downloading results.
    :param n_results: Numbers of results to display
    """
    logger.info(f"Rendering results for file: {file_name}")
    result_df = result.to_dataframe().head(n_results)
    result_df = result_df.set_index("Rank")
    st.divider()
    st.subheader(file_name)
    
    # Calculate library-specific sizes for display
    # Get the library-specific background size from the enrichment object
    filtered_terms = [term for term in result.gene_set_library.library 
                     if result.min_term_size <= term["size"] <= result.max_term_size]
    
    # Calculate unique genes from filtered terms only
    filtered_unique_genes = set()
    for term in filtered_terms:
        filtered_unique_genes.update(term["genes"])
    
    # Calculate library-specific background size
    library_background_size = len(result.background_gene_set.genes & filtered_unique_genes)
    
    # Calculate library-specific input size (genes that are in both input and library background)
    library_input_size = len(result.gene_set.genes & filtered_unique_genes)
    
    # Display library-specific size information
    st.caption(f"{library_input_size}/{result.gene_set.size} genes, {library_background_size}/{result.background_gene_set.size} background")
    
    # Check if there are any results to display
    if result_df.empty:
        st.warning("⚠️ No enrichment results found.")
        return
    
    table, bar = st.tabs(["Results", "Bar chart"])
    with table:
        render_table(result_df)
    with bar:
        render_barchart(result_df, file_name)

    st.markdown(
        f'Download results as {download_link(result.to_tsv(), file_name, "tsv")}, {download_link(result.to_json(), file_name, "json")}',
        unsafe_allow_html=True,
    )


def render_validation() -> None:
    """
    Validate and render the input gene set information.

    This function checks the `gene_set` in the session state for duplicates and invalid entries.
    It then provides a feedback to the user in the Streamlit app on the validation results.
    """
    logger.info("Validating and rendering the input gene set information.")
    if "gene_set" in state and state.gene_set is not None:
        total = state.gene_set.size
        dups = len(state.gene_set.validation["duplicates"])
        non_gene = len(state.gene_set.validation["non_valid"])
        if dups:
            dups_st = f", ⚠️ {dups} duplicates"
        else:
            dups_st = ""
        if non_gene:
            non_gene_st = f", ⛔ {non_gene} non valid"
        else:
            non_gene_st = ""
        caption = f"{total} genes{dups_st}{non_gene_st}"
        with st.expander(caption):
            if dups:
                st.data_editor(
                    pd.json_normalize(state.gene_set.validation)["duplicates"],
                    column_config={
                        "duplicates": st.column_config.ListColumn(
                            dups_st[2:], width="large"
                        ),
                    },
                    hide_index=True,
                )
            if non_gene:
                st.data_editor(
                    pd.json_normalize(state.gene_set.validation)["non_valid"],
                    column_config={
                        "non_valid": st.column_config.ListColumn(
                            non_gene_st[2:], width="large"
                        ),
                    },
                    hide_index=True,
                )
            if not dups and not non_gene:
                st.success("✅ All genes are valid")


def render_iter_table(result: pd.DataFrame) -> None:
    """
    Render a styled DataFrame of iterative enrichment results.

    :param result: DataFrame indexed by 'iteration' with columns ['term', 'p-value', 'overlap_size', 'genes'].
    :type result: pandas.DataFrame
    """
    logger.info("Rendering iterative results table.")

    # Create a copy to avoid modifying the original
    df = result.copy()
    
    # Convert underscores to spaces in term names for better readability
    if 'Term' in df.columns and not df.empty:
        df['Term'] = df['Term'].fillna('').astype(str).str.replace('_', ' ', regex=False)

    # Apply custom formatting to p-value columns
    def custom_format(n):
        if n == "" or pd.isna(n):
            return ""
        try:
            n = float(n)
            if n > 0.001:
                return f"{n:.3f}"
            return f"{n:.2e}"
        except (ValueError, TypeError):
            return str(n)

    # Format all p-value related columns
    format_dict = {}
    if "iteration p-value" in df.columns:
        format_dict["iteration p-value"] = custom_format
    if "Full list p-value" in df.columns:
        format_dict["Full list p-value"] = custom_format
    if "Regular FDR" in df.columns:
        format_dict["Regular FDR"] = custom_format

    styled = df.style.format(format_dict)

    # Build column config
    column_config = {
        "Term": "Term",
        "iteration p-value": "Iteration P-value",
        "iteration overlapping genes": "Iteration Overlapping Genes",
        "iteration -log(p-value)": "Iteration -log(p-value)",
        "Genes removed for next iteration": "Genes Removed for Next Iteration",
        "Description": None,
        "Library": None,
    }
    
    # Add regular enrichment columns if they exist
    if "Full list overlapping genes" in df.columns:
        column_config["Full list overlapping genes"] = "Full List Overlapping Genes"
    if "Full list p-value" in df.columns:
        column_config["Full list p-value"] = "Full List P-value"
    if "Regular FDR" in df.columns:
        column_config["Regular FDR"] = "Full List FDR"
    if "Full list overlapping genes (gene names)" in df.columns:
        column_config["Full list overlapping genes (gene names)"] = "Full List Overlapping Genes (Gene Names)"

    st.dataframe(
        styled,
        use_container_width=True,
        column_config=column_config,
    )


def render_iter_barchart(result: pd.DataFrame, file_name: str = "") -> None:
    """
    Render a bar chart visualization of iterative enrichment results.

    :param result: DataFrame indexed by 'iteration' with at least a 'p-value' column.
    :type result: pandas.DataFrame
    :param file_name: Optional file name to create unique chart keys.
    :type file_name: str
    """
    logger.info("Rendering iterative bar chart.")

    # Prepare bar plot data - use iteration p-value if available, otherwise fall back to p-value
    p_value_col = "iteration p-value" if "iteration p-value" in result.columns else "p-value"
    columns_to_select = ["Iteration", "Term", p_value_col]
    bar = result.reset_index()[columns_to_select].copy()
    
    # Convert underscores to spaces in term names for better readability
    if not bar.empty:
        bar['Term'] = bar['Term'].fillna('').astype(str).str.replace('_', ' ', regex=False)
    
    bar["-log10(p-value)"] = bar[p_value_col].apply(
        lambda x: -log10(x) if x and x > 0 else None
    )
    bar = bar.sort_values(by=[p_value_col], ascending=False)

    hover_col = p_value_col if p_value_col in bar.columns else "p-value"
    fig = px.bar(
        bar,
        x="-log10(p-value)",
        y="Term",
        orientation="h",
        hover_data=[hover_col],
        labels={"Term": "Term", "-log10(p-value)": "-log10(p-value)"},
        title="Iterative Enrichment p-value per Iteration",
    )
    st.plotly_chart(fig, use_container_width=True, key=f"iter_barchart_{file_name}")


def render_iter_results(result: IterativeEnrichment, file_name: str) -> None:
    """
    Render a results section for iterative enrichment within the Streamlit app.

    :param result: IterativeEnrichment object containing iteration records.
    :type result: src.iter_enrichment.IterativeEnrichment
    :param file_name: Name of the library or section header.
    :type file_name: str
    """
    logger.info(f"Rendering iterative results for {file_name}.")

    st.divider()
    st.subheader(file_name)
    
    # Calculate and display library-specific size information for the initial gene set
    # Get the first iteration's enrichment object to show initial sizes
    if result._iteration_enrichments:
        initial_enrichment = result._iteration_enrichments[0]
        
        # Calculate library-specific sizes for display
        filtered_terms = [term for term in initial_enrichment.gene_set_library.library 
                         if initial_enrichment.min_term_size <= term["size"] <= initial_enrichment.max_term_size]
        
        # Calculate unique genes from filtered terms only
        filtered_unique_genes = set()
        for term in filtered_terms:
            filtered_unique_genes.update(term["genes"])
        
        # Calculate library-specific background size
        library_background_size = len(initial_enrichment.background_gene_set.genes & filtered_unique_genes)
        
        # Calculate library-specific input size (genes that are in both input and library background)
        library_input_size = len(initial_enrichment.gene_set.genes & filtered_unique_genes)
        
        # Display library-specific size information
        st.caption(f"Library specific input: {library_input_size}/{initial_enrichment.gene_set.size} genes; Library specific background: {library_background_size}/{initial_enrichment.background_gene_set.size}")
        
        # Show final remaining genes if different from initial
        if result.results:
            # Calculate final remaining genes by starting with library-specific input size
            # and subtracting genes removed in each iteration
            final_remaining = library_input_size - sum(len(record.get("genes", [])) for record in result.results)
            if final_remaining != library_input_size:
                st.caption(f"Final: {final_remaining} genes remaining after {len(result.results)} iterations")

    # Tabs for table and chart
    table_tab, bar_tab = st.tabs(["Iterations", "Bar chart"])
    df = result.to_dataframe()

    if df.empty or "Iteration" not in df.columns:
        logger.warning("No iterative enrichment results.")
        st.warning("No iterative enrichment results.")
    else:
        df = df.set_index("Iteration")
        with table_tab:
            render_iter_table(df)

        with bar_tab:
            render_iter_barchart(df, file_name)

        # Download links
        st.markdown(
            f'Download iterative results as {download_link(result.to_tsv(), file_name, "tsv")}, '
            f'{download_link(result.to_json(), file_name, "json")}',
            unsafe_allow_html=True,
        )


def calculate_ora_igea_comparison(
    regular_enrichments: Dict[str, Enrichment],
    iterative_enrichments: Dict[str, IterativeEnrichment],
    p_threshold: float = 0.01,
    fdr_threshold: float = 0.05,
) -> Dict:
    """
    Calculate comparison statistics between ORA (regular enrichment) and iGEA (iterative enrichment).
    ORA is split into two sets: one based on raw p-value threshold and one based on FDR threshold.
    
    Args:
        regular_enrichments: Dictionary mapping library names to Enrichment objects
        iterative_enrichments: Dictionary mapping library names to IterativeEnrichment objects
        p_threshold: P-value threshold for significance (default: 0.01)
        fdr_threshold: FDR threshold for significance (default: 0.05)
        
    Returns:
        Dictionary containing comparison statistics including:
        - genes_ora_raw: set of unique genes in ORA (raw p-value)
        - genes_ora_fdr: set of unique genes in ORA (FDR)
        - genes_igea: set of unique genes in iGEA
        - genes_overlap_raw: set of genes in both ORA (raw) and iGEA
        - genes_overlap_fdr: set of genes in both ORA (FDR) and iGEA
        - terms_ora_raw: set of significant terms in ORA (raw p-value)
        - terms_ora_fdr: set of significant terms in ORA (FDR)
        - terms_igea: set of significant terms in iGEA
        - terms_overlap_raw: set of terms in both ORA (raw) and iGEA
        - terms_overlap_fdr: set of terms in both ORA (FDR) and iGEA
        - stats: dictionary with counts and percentages
    """
    # Collect genes and terms from regular enrichment (ORA) - split by raw p-value and FDR
    # If regular enrichment is not available, use iteration 1 from iterative enrichment as ORA
    genes_ora_raw: Set[str] = set()
    terms_ora_raw: Set[str] = set()
    genes_ora_fdr: Set[str] = set()
    terms_ora_fdr: Set[str] = set()
    
    def _process_ora_result(result: dict, genes_raw: Set[str], terms_raw: Set[str], 
                           genes_fdr: Set[str], terms_fdr: Set[str]) -> None:
        """Helper function to process a single ORA result and add to appropriate sets."""
        # Get p-value
        pval = result.get("p-value", 1.0)
        if isinstance(pval, str) and pval != "":
            try:
                pval_float = float(pval)
            except ValueError:
                pval_float = 1.0
        elif isinstance(pval, (int, float)):
            pval_float = pval
        else:
            pval_float = 1.0
        
        # Get FDR
        fdr = result.get("fdr", 1.0)
        if isinstance(fdr, str) and fdr != "":
            try:
                fdr_float = float(fdr)
            except ValueError:
                fdr_float = 1.0
        elif isinstance(fdr, (int, float)):
            fdr_float = fdr
        else:
            fdr_float = 1.0
        
        # Check significance using raw p-value threshold
        is_significant_raw = pval_float <= p_threshold
        
        # Check significance using FDR threshold
        is_significant_fdr = fdr_float <= fdr_threshold
        
        term_name = result.get("term", "")
        if not term_name:
            return
        
        # Extract genes from overlap
        overlap_genes = result.get("overlap", [])
        if isinstance(overlap_genes, str):
            overlap_genes = [g.strip() for g in overlap_genes.split(",") if g.strip()]
        elif not isinstance(overlap_genes, list):
            overlap_genes = []
        
        # Add to raw p-value sets
        if is_significant_raw:
            terms_raw.add(term_name)
            genes_raw.update(overlap_genes)
        
        # Add to FDR sets
        if is_significant_fdr:
            terms_fdr.add(term_name)
            genes_fdr.update(overlap_genes)
    
    # First, try to get ORA from regular enrichment
    if regular_enrichments:
        for lib_name, enrich in regular_enrichments.items():
            if enrich and enrich.results:
                for result in enrich.results:
                    _process_ora_result(result, genes_ora_raw, terms_ora_raw, genes_ora_fdr, terms_ora_fdr)
    
    # If no regular enrichment results, use iteration 1's full enrichment results as ORA
    # ORA should include ALL significant terms from iteration 1, not just the top term
    if not genes_ora_raw and iterative_enrichments:
        for lib_name, iter_enrich in iterative_enrichments.items():
            # Use _regular_enrichment which contains ALL terms from iteration 1
            if iter_enrich and hasattr(iter_enrich, '_regular_enrichment') and iter_enrich._regular_enrichment:
                regular_enrich = iter_enrich._regular_enrichment
                if regular_enrich.results:
                    for result in regular_enrich.results:
                        _process_ora_result(result, genes_ora_raw, terms_ora_raw, genes_ora_fdr, terms_ora_fdr)
    
    # Collect genes from iterative enrichment (iGEA) - ALL iterations (including iteration 1)
    # iGEA represents the combined output of all iterations
    genes_igea: Set[str] = set()
    terms_igea: Set[str] = set()
    
    for lib_name, iter_enrich in iterative_enrichments.items():
        if iter_enrich and iter_enrich.results:
            for record in iter_enrich.results:
                # Include ALL iterations in iGEA (including iteration 1)
                # Check if term is significant using raw p-value threshold only
                iter_pval = record.get("iteration p-value", record.get("p-value", 1.0))
                if isinstance(iter_pval, str) and iter_pval != "":
                    try:
                        iter_pval_float = float(iter_pval)
                    except ValueError:
                        iter_pval_float = 1.0
                elif isinstance(iter_pval, (int, float)):
                    iter_pval_float = iter_pval
                else:
                    iter_pval_float = 1.0
                
                # Use only p-value threshold (not FDR) for significance
                is_significant = iter_pval_float <= p_threshold
                
                if is_significant:
                    # Use Term_raw for matching with ORA (raw term names), fallback to term if not available
                    term_name = record.get("Term_raw", record.get("term", ""))
                    if not term_name:
                        # If Term_raw is not available, try to get from Term (formatted) and convert back
                        formatted_term = record.get("Term", "")
                        if formatted_term:
                            # Convert formatted term back to raw (spaces to underscores)
                            # This is a fallback - Term_raw should always be present
                            term_name = formatted_term.replace(" ", "_")
                    if term_name:
                        terms_igea.add(term_name)
                    
                    # Extract genes removed for next iteration (these are the genes in the term)
                    genes_removed = record.get("Genes removed for next iteration", [])
                    if isinstance(genes_removed, list):
                        genes_igea.update(genes_removed)
                    elif isinstance(genes_removed, str):
                        # Handle comma-separated string
                        genes_igea.update([g.strip() for g in genes_removed.split(",") if g.strip()])
    
    # Calculate overlaps - separate for raw p-value and FDR
    genes_overlap_raw = genes_ora_raw & genes_igea
    genes_overlap_fdr = genes_ora_fdr & genes_igea
    genes_ora_raw_only = genes_ora_raw - genes_igea
    genes_ora_fdr_only = genes_ora_fdr - genes_igea
    genes_igea_only = genes_igea - (genes_ora_raw | genes_ora_fdr)
    
    terms_overlap_raw = terms_ora_raw & terms_igea
    terms_overlap_fdr = terms_ora_fdr & terms_igea
    terms_ora_raw_only = terms_ora_raw - terms_igea
    terms_ora_fdr_only = terms_ora_fdr - terms_igea
    terms_igea_only = terms_igea - (terms_ora_raw | terms_ora_fdr)
    
    # Debug logging to help diagnose counting issues
    logger.info(f"ORA vs iGEA comparison (p_threshold={p_threshold}, fdr_threshold={fdr_threshold}):")
    logger.info(f"  ORA (raw) terms: {len(terms_ora_raw)} (sample: {list(terms_ora_raw)[:3] if terms_ora_raw else []})")
    logger.info(f"  ORA (FDR) terms: {len(terms_ora_fdr)} (sample: {list(terms_ora_fdr)[:3] if terms_ora_fdr else []})")
    logger.info(f"  iGEA terms: {len(terms_igea)} (sample: {list(terms_igea)[:3] if terms_igea else []})")
    logger.info(f"  Overlapping terms (raw): {len(terms_overlap_raw)} (sample: {list(terms_overlap_raw)[:3] if terms_overlap_raw else []})")
    logger.info(f"  Overlapping terms (FDR): {len(terms_overlap_fdr)} (sample: {list(terms_overlap_fdr)[:3] if terms_overlap_fdr else []})")
    
    # Calculate statistics
    total_genes_ora_raw = len(genes_ora_raw)
    total_genes_ora_fdr = len(genes_ora_fdr)
    total_genes_igea = len(genes_igea)
    total_genes_overlap_raw = len(genes_overlap_raw)
    total_genes_overlap_fdr = len(genes_overlap_fdr)
    
    total_terms_ora_raw = len(terms_ora_raw)
    total_terms_ora_fdr = len(terms_ora_fdr)
    total_terms_igea = len(terms_igea)
    total_terms_overlap_raw = len(terms_overlap_raw)
    total_terms_overlap_fdr = len(terms_overlap_fdr)
    
    stats = {
        "genes": {
            "ora_raw_total": total_genes_ora_raw,
            "ora_fdr_total": total_genes_ora_fdr,
            "igea_total": total_genes_igea,
            "overlap_raw": total_genes_overlap_raw,
            "overlap_fdr": total_genes_overlap_fdr,
            "ora_raw_only": len(genes_ora_raw_only),
            "ora_fdr_only": len(genes_ora_fdr_only),
            "igea_only": len(genes_igea_only),
        },
        "terms": {
            "ora_raw_total": total_terms_ora_raw,
            "ora_fdr_total": total_terms_ora_fdr,
            "igea_total": total_terms_igea,
            "overlap_raw": total_terms_overlap_raw,
            "overlap_fdr": total_terms_overlap_fdr,
            "ora_raw_only": len(terms_ora_raw_only),
            "ora_fdr_only": len(terms_ora_fdr_only),
            "igea_only": len(terms_igea_only),
        }
    }
    
    return {
        "genes_ora_raw": genes_ora_raw,
        "genes_ora_fdr": genes_ora_fdr,
        "genes_igea": genes_igea,
        "genes_overlap_raw": genes_overlap_raw,
        "genes_overlap_fdr": genes_overlap_fdr,
        "genes_ora_raw_only": genes_ora_raw_only,
        "genes_ora_fdr_only": genes_ora_fdr_only,
        "genes_igea_only": genes_igea_only,
        "terms_ora_raw": terms_ora_raw,
        "terms_ora_fdr": terms_ora_fdr,
        "terms_igea": terms_igea,
        "terms_overlap_raw": terms_overlap_raw,
        "terms_overlap_fdr": terms_overlap_fdr,
        "terms_ora_raw_only": terms_ora_raw_only,
        "terms_ora_fdr_only": terms_ora_fdr_only,
        "terms_igea_only": terms_igea_only,
        "stats": stats,
    }


def render_ora_igea_comparison(
    regular_enrichments: Dict[str, Enrichment],
    iterative_enrichments: Dict[str, IterativeEnrichment],
    p_threshold: float = 0.01,
    fdr_threshold: float = 0.05,
) -> None:
    """
    Render comparison statistics between ORA and iGEA results.
    
    Args:
        regular_enrichments: Dictionary mapping library names to Enrichment objects
        iterative_enrichments: Dictionary mapping library names to IterativeEnrichment objects
        p_threshold: P-value threshold for significance
        fdr_threshold: FDR threshold for significance
    """
    # Check if we have iterative enrichment results (required)
    if not iterative_enrichments:
        st.info("ℹ️ Iterative enrichment results are required for comparison.")
        return
    
    # Check if we have actual iterative results
    has_iterative_results = any(
        iter_enrich and iter_enrich.results 
        for iter_enrich in iterative_enrichments.values()
    )
    
    if not has_iterative_results:
        st.info("ℹ️ Iterative enrichment results with data are required for comparison.")
        return
    
    # Check if we have ORA data (either from regular enrichment or iteration 1)
    has_regular_results = False
    if regular_enrichments:
        has_regular_results = any(
            enrich and enrich.results 
            for enrich in regular_enrichments.values()
        )
    
    # If no regular enrichment, check if we have iteration 1 to use as ORA
    if not has_regular_results:
        has_iteration_1 = any(
            iter_enrich and iter_enrich.results and 
            any(r.get("Iteration") == 1 for r in iter_enrich.results)
            for iter_enrich in iterative_enrichments.values()
        )
        if not has_iteration_1:
            st.info("ℹ️ Either regular enrichment results or iterative enrichment with at least iteration 1 is required for comparison.")
            return
    
    comparison = calculate_ora_igea_comparison(
        regular_enrichments, iterative_enrichments, p_threshold, fdr_threshold
    )
    
    stats = comparison["stats"]
    
    st.markdown("### 📊 ORA vs iGEA Comparison Statistics")
    
    # Create comparison table
    gene_stats = stats["genes"]
    term_stats = stats["terms"]
    
    # Display as expanded table with ORA split into raw p-value and FDR
    st.markdown("""
    | | ORA (raw p-value) | ORA (FDR) | iGEA | Overlap (raw) | Overlap (FDR) |
    |---|---|---|---|---|---|
    | **Genes** | {ora_raw_genes} | {ora_fdr_genes} | {igea_genes} | {overlap_raw_genes} | {overlap_fdr_genes} |
    | **Terms** | {ora_raw_terms} | {ora_fdr_terms} | {igea_terms} | {overlap_raw_terms} | {overlap_fdr_terms} |
    """.format(
        ora_raw_genes=gene_stats["ora_raw_total"],
        ora_fdr_genes=gene_stats["ora_fdr_total"],
        igea_genes=gene_stats["igea_total"],
        overlap_raw_genes=gene_stats["overlap_raw"],
        overlap_fdr_genes=gene_stats["overlap_fdr"],
        ora_raw_terms=term_stats["ora_raw_total"],
        ora_fdr_terms=term_stats["ora_fdr_total"],
        igea_terms=term_stats["igea_total"],
        overlap_raw_terms=term_stats["overlap_raw"],
        overlap_fdr_terms=term_stats["overlap_fdr"]
    ))
    
    # Venn diagram data for download
    st.markdown("---")
    st.markdown("#### 📥 Venn Diagram Data")
    
    # Create TSV with Venn diagram data - expanded to include raw and FDR sets
    venn_data_rows = [
        "Category\tSet\tCount\tItems",
        f"Genes\tORA (raw) Only\t{len(comparison['genes_ora_raw_only'])}\t{', '.join(sorted(comparison['genes_ora_raw_only']))}",
        f"Genes\tORA (FDR) Only\t{len(comparison['genes_ora_fdr_only'])}\t{', '.join(sorted(comparison['genes_ora_fdr_only']))}",
        f"Genes\tiGEA Only\t{len(comparison['genes_igea_only'])}\t{', '.join(sorted(comparison['genes_igea_only']))}",
        f"Genes\tOverlap (raw)\t{len(comparison['genes_overlap_raw'])}\t{', '.join(sorted(comparison['genes_overlap_raw']))}",
        f"Genes\tOverlap (FDR)\t{len(comparison['genes_overlap_fdr'])}\t{', '.join(sorted(comparison['genes_overlap_fdr']))}",
        f"Terms\tORA (raw) Only\t{len(comparison['terms_ora_raw_only'])}\t{', '.join(sorted(comparison['terms_ora_raw_only']))}",
        f"Terms\tORA (FDR) Only\t{len(comparison['terms_ora_fdr_only'])}\t{', '.join(sorted(comparison['terms_ora_fdr_only']))}",
        f"Terms\tiGEA Only\t{len(comparison['terms_igea_only'])}\t{', '.join(sorted(comparison['terms_igea_only']))}",
        f"Terms\tOverlap (raw)\t{len(comparison['terms_overlap_raw'])}\t{', '.join(sorted(comparison['terms_overlap_raw']))}",
        f"Terms\tOverlap (FDR)\t{len(comparison['terms_overlap_fdr'])}\t{', '.join(sorted(comparison['terms_overlap_fdr']))}",
    ]
    venn_data_tsv = "\n".join(venn_data_rows)
    
    st.markdown(
        f"Download Venn diagram data: {download_link(venn_data_tsv, 'ora_igea_venn_data', 'tsv')}",
        unsafe_allow_html=True,
    )


def render_network(
    dot: str, 
    title: str = "Iterative Enrichment Network",
    gene_count: int = None,
    library_list: List[str] = None
) -> None:
    """
    Render a bipartite network graph of iterative enrichment across libraries.

    :param dot: Graphviz DOT-format string.
    :type dot: str
    :param title: Header title for the network graph.
    :type title: str
    :param gene_count: Number of genes in the original input gene list.
    :type gene_count: int, optional
    :param library_list: List of library names that were tested.
    :type library_list: list, optional
    """
    logger.info("Rendering iterative network graph.")

    st.divider()
    st.subheader(title)
    
    # Count edges to check network size
    from ui.dot_utils import count_edges_in_dot
    edge_count = count_edges_in_dot(dot)
    
    # Check if network is too large for Plotly rendering
    if edge_count <= 500:
        # Render Plotly network as usual
        st.plotly_chart(dot_to_plotly(dot), use_container_width=True, key="network_chart")
    else:
        # Show warning and skip Plotly rendering
        st.warning(
            f"⚠️ Network too large ({edge_count} edges). "
            "Plotly visualization skipped to prevent crashes. "
            "Download the DOT file to view in external tools like Graphviz."
        )
    
    # Generate AI prompt - this is the source of truth for both clipboard and download
    ai_prompt = generate_ai_analysis_prompt(dot, gene_count=gene_count, library_list=library_list)
    
    # Display AI prompt section with platform buttons
    st.markdown("---")
    st.subheader("🤖 AI Analysis Prompt")
    st.caption("Use this prompt with conversational AI systems to analyze your enrichment network")
    
    # Display prompt in text area (selectable for easy copying)
    # Note: The text_area displays the same ai_prompt that will be downloaded
    prompt_key = f"ai_prompt_{hash(dot) % 10000}"  # Unique key based on dot content
    st.text_area(
        "AI Analysis Prompt:",
        value=ai_prompt,  # Same text as downloadable file
        height=400,
        key=prompt_key,
        help="Select all text (Ctrl+A / Cmd+A) and copy (Ctrl+C / Cmd+C) to paste into your preferred AI platform"
    )
    
    # AI Platform buttons
    st.markdown("**Open AI Platform:**")
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.link_button(
            "💬 ChatGPT",
            "https://chat.openai.com/",
            use_container_width=True,
            help="Opens ChatGPT in a new tab. Copy the prompt from the text box above and paste it into the chat."
        )
    
    with col2:
        st.link_button(
            "🧠 Claude",
            "https://claude.ai/",
            use_container_width=True,
            help="Opens Claude in a new tab. Copy the prompt from the text box above and paste it into the chat."
        )
    
    with col3:
        st.link_button(
            "✨ Gemini",
            "https://gemini.google.com/",
            use_container_width=True,
            help="Opens Google Gemini in a new tab. Copy the prompt from the text box above and paste it into the chat."
        )
    
    # Instructions
    st.info("💡 **How to use:** 1) Select all text in the prompt box above (Ctrl+A / Cmd+A) 2) Copy it (Ctrl+C / Cmd+C) 3) Click your preferred AI platform button 4) Paste the prompt (Ctrl+V / Cmd+V) into the chat")
    
    # Always offer downloads with clear labels
    st.markdown(
        f'**Downloads:** Download network in DOT format {download_link(dot, "iterative_network", "dot")} | Download AI prompt as text file {download_link(ai_prompt, "ai_analysis_prompt", "txt")}',
        unsafe_allow_html=True,
    )


def generate_ai_analysis_prompt(
    dot_content: str, 
    gene_count: int = None, 
    library_list: List[str] = None
) -> str:
    """
    Generate an AI analysis prompt based on the DOT network content with library source annotations.
    
    :param dot_content: The DOT network content
    :param gene_count: Number of genes in the original input gene list
    :param library_list: List of library names that were tested
    :return: Formatted AI analysis prompt with library information
    """
    import re
    from ui.dot_utils import load_library_colors
    
    # Load library color mapping
    color_to_library = load_library_colors()
    
    # Parse the DOT content to add library annotations
    lines = dot_content.split('\n')
    annotated_lines = []
    
    for line in lines:
        line = line.strip()
        
        # Check if this is a term node line
        if 'type="term"' in line and 'fillcolor=' in line:
            # Extract the fillcolor
            color_match = re.search(r'fillcolor="([^"]+)"', line)
            if color_match:
                color = color_match.group(1)
                library_name = color_to_library.get(color, "Unknown Library")
                
                # Add library information to the node attributes
                if 'library=' not in line:
                    # Insert library attribute before the closing bracket
                    line = line.rstrip('];') + f', library="{library_name}"];'
        
        annotated_lines.append(line)
    
    # Create the annotated DOT content
    annotated_dot = '\n'.join(annotated_lines)
    
    prompt = """You are a computational biologist analyzing an iterative gene set enrichment network. Iterative means that for each library source, the top enriched gene set is found, saved and the genes are removed from the initial gene list. The remaining genes are then tested for enrichment again. Thus, each gene appears only once per library source tested but can be linked to multiple terms that originate from different libraries. The results are provided as a DOT network represents the relationships between genes and gene sets(=terms).

"""
    
    # Add gene count and library list information if provided
    if gene_count is not None or library_list is not None:
        info_parts = []
        if gene_count is not None:
            info_parts.append(f"The original gene list included {gene_count} genes")
        if library_list is not None and len(library_list) > 0:
            library_str = ", ".join(library_list)
            if gene_count is not None:
                info_parts.append(f"and was tested across the following libraries: {library_str}")
            else:
                info_parts.append(f"The analysis was tested across the following libraries: {library_str}")
        
        if info_parts:
            prompt += f"{' '.join(info_parts)}.\n\n"
    
    prompt += """**NETWORK STRUCTURE:**
"""
    
    # Add the annotated DOT content
    prompt += annotated_dot + "\n\n"
    
    # Add the rest of the prompt
    prompt += """**NETWORK INTERPRETATION GUIDE:**

**Node Types:**
- **Gene nodes** (type="gene"): Individual genes from the input gene set
- **Term nodes** (type="term"): Enriched biological pathways/processes/gene sets
- **Library Source**: Each term node includes a "library" attribute indicating the source database

**Library Sources and Their Biological Context:**
- **H: Hallmark gene sets**: Curated gene sets representing well-defined biological states or processes
- **C2: BioCarta**: Canonical pathways from BioCarta database
- **C2: KEGG MEDICUS**: Metabolic and signaling pathways from KEGG
- **C2: Pathway Interaction Database**: Curated human signaling pathways
- **C2: Reactome Pathways**: Expert-curated biological pathways
- **C2: WikiPathways**: Community-curated biological pathways
- **C5: Gene Ontology: Biological Process**: Biological processes from Gene Ontology
- **C5: Gene Ontology: Cellular Component**: Cellular components from Gene Ontology
- **C5: Gene Ontology: Molecular Function**: Molecular functions from Gene Ontology
- **C5: Human Phenotype Ontology**: Human phenotypes and diseases
- **Protein Interaction**: Protein-protein interaction networks

**Edge Connections:**
- Each edge (gene -- term) represents a gene's membership in that biological term
- More connections indicate genes that are part of multiple enriched processes
- The network topology shows which genes are central to the biological response

**Iteration Information:**
- Term nodes are labeled with iteration numbers (term_1_, term_2_, etc.)
- Earlier iterations (lower numbers) represent the most statistically significant enrichments
- Ignore iteration number for biological interpretation, instead focus on network hubs densely connected for interesting biological insights

**ANALYSIS REQUEST:**

Please analyze this network and provide:

1. **Key Biological Insights:**
   - What are the most significant biological processes/pathways identified?
   - Which genes appear to be central to the biological response?

2. **Network Topology Analysis:**
   - Which genes are "hub" genes (connected to multiple terms)?
   - Are there distinct clusters or modules in the network?
   - What does the connectivity pattern suggest about biological organization?

3. **Biological Hypothesis:**
   - Based on the network structure, what biological hypothesis can you generate?
   - How do the different library sources support or complement each other?
   - Try to formulate one hypothesis based on the most central terms and how they possibly related to one another.
   
4. **Estimated Experimental Context**
   - Can you hypothesize on what was the experiment that generated this result?

**RESPONSE STRUCTURE:**
- **Executive Summary** (2-3 sentences)
- **Key Biological Processes Identified** (by library source)
- **Network Topology Insights**
- **Proposed Biological Hypothesis**
- **Estimated Experimental Context**

**IMPORTANT NOTES:**
- Focus on biological interpretation, not statistical significance or iteration number
- Consider the functional relationships between connected genes and terms
- Look for unexpected connections that might reveal novel or unexpected biology
- Consider the broader biological context and literature
- Note that protein interaction library could shed light on physical interaction between genes while terms indicate functional interaction
- Different library sources may provide complementary biological insights"""
    
    return prompt


def generate_structured_network_analysis(dot_content: str) -> str:
    """
    Generate a structured table format of network data for ChatGPT analysis.
    This provides the network information in a more parseable, tabular format.
    
    :param dot_content: The DOT network content
    :return: Structured table format for AI analysis
    """
    import re
    from collections import defaultdict
    
    # Parse the DOT content
    nodes = {}
    edges = []
    gene_connections = defaultdict(set)
    term_connections = defaultdict(set)
    
    # Parse DOT content
    lines = dot_content.split('\n')
    for line in lines:
        line = line.strip()
        
        # Parse node definitions
        if '[' in line and ']' in line and '--' not in line:
            match = re.match(r'"([^"]+)"\s*\[([^\]]+)\]', line)
            if match:
                node_id = match.group(1)
                attrs_str = match.group(2)
                
                # Parse attributes
                attrs = {}
                for attr in attrs_str.split(','):
                    attr = attr.strip()
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attrs[key.strip()] = value.strip().strip('"')
                
                nodes[node_id] = attrs
        
        # Parse edges
        elif '--' in line:
            match = re.match(r'"([^"]+)"\s*--\s*"([^"]+)"', line)
            if match:
                source = match.group(1)
                target = match.group(2)
                edges.append((source, target))
                
                # Track connections
                gene_connections[source].add(target)
                gene_connections[target].add(source)
                term_connections[source].add(target)
                term_connections[target].add(source)
    
    # Analyze network structure
    gene_nodes = {k: v for k, v in nodes.items() if v.get('type') == 'gene'}
    term_nodes = {k: v for k, v in nodes.items() if v.get('type') == 'term'}
    
    # Calculate centrality metrics
    gene_centrality = {}
    for gene_id, gene_data in gene_nodes.items():
        connections = len(gene_connections.get(gene_id, set()))
        gene_centrality[gene_data.get('label', gene_id)] = connections
    
    # Group terms by iteration
    term_iterations = defaultdict(list)
    for term_id, term_data in term_nodes.items():
        match = re.match(r'term_(\d+)_', term_id)
        if match:
            iteration = int(match.group(1))
            term_iterations[iteration].append({
                'id': term_id,
                'name': term_data.get('label', term_id),
                'color': term_data.get('fillcolor', 'unknown'),
                'connections': len(term_connections.get(term_id, set()))
            })
    
    # Create structured output
    output = f"""# ITERATIVE GENE SET ENRICHMENT NETWORK ANALYSIS
# Structured Data Format for AI Analysis

## NETWORK SUMMARY
Total Genes: {len(gene_nodes)}
Total Biological Terms: {len(term_nodes)}
Total Connections: {len(edges)}
Total Iterations: {len(term_iterations)}
Average Connections per Gene: {len(edges) / len(gene_nodes) if gene_nodes else 0:.2f}

## GENE CENTRALITY TABLE
# Rank | Gene Name | Connections | Biological Terms
"""
    
    # Add gene centrality table
    sorted_genes = sorted(gene_centrality.items(), key=lambda x: x[1], reverse=True)
    for i, (gene_name, centrality) in enumerate(sorted_genes, 1):
        # Get the terms this gene is connected to
        gene_id = f"gene_{gene_name}"
        connected_terms = []
        for term_id in gene_connections.get(gene_id, set()):
            if term_id in term_nodes:
                term_name = term_nodes[term_id].get('label', term_id)
                connected_terms.append(term_name)
        
        output += f"{i:4d} | {gene_name:15s} | {centrality:11d} | {', '.join(connected_terms[:3])}"
        if len(connected_terms) > 3:
            output += f" (+{len(connected_terms)-3} more)"
        output += "\n"
    
    output += "\n## ITERATION SEQUENCE TABLE\n"
    output += "# Iteration | Term Name | Genes Connected | Significance Level\n"
    
    # Add iteration sequence table
    for iteration in sorted(term_iterations.keys()):
        terms = term_iterations[iteration]
        for term in terms:
            significance = "HIGH" if iteration == 1 else "MEDIUM" if iteration <= 3 else "LOW"
            output += f"{iteration:9d} | {term['name']:30s} | {term['connections']:14d} | {significance}\n"
    
    output += "\n## CONNECTION MATRIX\n"
    output += "# Gene -> Biological Terms\n"
    
    # Add connection matrix
    for gene_name, centrality in sorted_genes[:10]:  # Top 10 genes
        gene_id = f"gene_{gene_name}"
        connected_terms = []
        for term_id in gene_connections.get(gene_id, set()):
            if term_id in term_nodes:
                term_name = term_nodes[term_id].get('label', term_id)
                connected_terms.append(term_name)
        
        output += f"{gene_name} -> {', '.join(connected_terms)}\n"
    
    output += "\n## BIOLOGICAL TERMS SUMMARY\n"
    output += "# Term Name | Iteration | Gene Count | Connected Genes\n"
    
    # Add biological terms summary
    for iteration in sorted(term_iterations.keys()):
        terms = term_iterations[iteration]
        for term in terms:
            term_id = term['id']
            connected_genes = []
            for gene_id in term_connections.get(term_id, set()):
                if gene_id in gene_nodes:
                    gene_name = gene_nodes[gene_id].get('label', gene_id)
                    connected_genes.append(gene_name)
            
            output += f"{term['name']:30s} | {iteration:9d} | {term['connections']:10d} | {', '.join(connected_genes[:5])}"
            if len(connected_genes) > 5:
                output += f" (+{len(connected_genes)-5} more)"
            output += "\n"
    
    output += "\n## ANALYSIS INSTRUCTIONS\n"
    output += """Please analyze this structured network data and provide:

1. **Executive Summary**: Key biological themes and experimental context
2. **Hub Gene Analysis**: Roles of highly connected genes in biological processes  
3. **Pathway Relationships**: How biological terms relate to each other
4. **Iteration Patterns**: Biological significance of the iteration sequence
5. **Biological Hypothesis**: Potential experimental conditions and follow-up studies
6. **Clinical Relevance**: Therapeutic targets and biomarkers

Focus on biological interpretation and actionable insights."""
    
    return output


def generate_json_network_analysis(dot_content: str) -> str:
    """
    Generate a JSON format of network data for ChatGPT analysis.
    This provides the network information in a structured JSON format that's easy to parse.
    
    :param dot_content: The DOT network content
    :return: JSON format for AI analysis
    """
    import re
    import json
    from collections import defaultdict
    
    # Parse the DOT content
    nodes = {}
    edges = []
    gene_connections = defaultdict(set)
    term_connections = defaultdict(set)
    
    # Parse DOT content
    lines = dot_content.split('\n')
    for line in lines:
        line = line.strip()
        
        # Parse node definitions
        if '[' in line and ']' in line and '--' not in line:
            match = re.match(r'"([^"]+)"\s*\[([^\]]+)\]', line)
            if match:
                node_id = match.group(1)
                attrs_str = match.group(2)
                
                # Parse attributes
                attrs = {}
                for attr in attrs_str.split(','):
                    attr = attr.strip()
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attrs[key.strip()] = value.strip().strip('"')
                
                nodes[node_id] = attrs
        
        # Parse edges
        elif '--' in line:
            match = re.match(r'"([^"]+)"\s*--\s*"([^"]+)"', line)
            if match:
                source = match.group(1)
                target = match.group(2)
                edges.append((source, target))
                
                # Track connections
                gene_connections[source].add(target)
                gene_connections[target].add(source)
                term_connections[source].add(target)
                term_connections[target].add(source)
    
    # Analyze network structure
    gene_nodes = {k: v for k, v in nodes.items() if v.get('type') == 'gene'}
    term_nodes = {k: v for k, v in nodes.items() if v.get('type') == 'term'}
    
    # Calculate centrality metrics
    gene_centrality = {}
    for gene_id, gene_data in gene_nodes.items():
        connections = len(gene_connections.get(gene_id, set()))
        gene_centrality[gene_data.get('label', gene_id)] = connections
    
    # Group terms by iteration
    term_iterations = defaultdict(list)
    for term_id, term_data in term_nodes.items():
        match = re.match(r'term_(\d+)_', term_id)
        if match:
            iteration = int(match.group(1))
            term_iterations[iteration].append({
                'id': term_id,
                'name': term_data.get('label', term_id),
                'color': term_data.get('fillcolor', 'unknown'),
                'connections': len(term_connections.get(term_id, set()))
            })
    
    # Create JSON structure
    network_data = {
        "analysis_type": "iterative_gene_set_enrichment",
        "network_summary": {
            "total_genes": len(gene_nodes),
            "total_biological_terms": len(term_nodes),
            "total_connections": len(edges),
            "total_iterations": len(term_iterations),
            "average_connections_per_gene": len(edges) / len(gene_nodes) if gene_nodes else 0
        },
        "gene_centrality": [
            {
                "gene_name": gene_name,
                "connections": centrality,
                "rank": i + 1,
                "connected_terms": [
                    {
                        "term_name": term_nodes[term_id].get('label', term_id),
                        "iteration": int(re.match(r'term_(\d+)_', term_id).group(1)) if re.match(r'term_(\d+)_', term_id) else None
                    }
                    for term_id in gene_connections.get(f"gene_{gene_name}", set())
                    if term_id in term_nodes
                ]
            }
            for i, (gene_name, centrality) in enumerate(sorted(gene_centrality.items(), key=lambda x: x[1], reverse=True))
        ],
        "iteration_sequence": [
            {
                "iteration": iteration,
                "significance_level": "HIGH" if iteration == 1 else "MEDIUM" if iteration <= 3 else "LOW",
                "terms": [
                    {
                        "term_name": term['name'],
                        "gene_count": term['connections'],
                        "connected_genes": [
                            gene_nodes[gene_id].get('label', gene_id)
                            for gene_id in term_connections.get(term['id'], set())
                            if gene_id in gene_nodes
                        ]
                    }
                    for term in terms
                ]
            }
            for iteration, terms in sorted(term_iterations.items())
        ],
        "biological_terms": [
            {
                "term_name": term_data.get('label', term_id),
                "iteration": int(re.match(r'term_(\d+)_', term_id).group(1)) if re.match(r'term_(\d+)_', term_id) else None,
                "gene_count": len(term_connections.get(term_id, set())),
                "connected_genes": [
                    gene_nodes[gene_id].get('label', gene_id)
                    for gene_id in term_connections.get(term_id, set())
                    if gene_id in gene_nodes
                ]
            }
            for term_id, term_data in term_nodes.items()
        ],
        "hub_genes": [
            {
                "gene_name": gene_name,
                "centrality_score": centrality,
                "biological_roles": [
                    term_nodes[term_id].get('label', term_id)
                    for term_id in gene_connections.get(f"gene_{gene_name}", set())
                    if term_id in term_nodes
                ]
            }
            for gene_name, centrality in sorted(gene_centrality.items(), key=lambda x: x[1], reverse=True)[:10]
        ]
    }
    
    # Create the output with instructions
    output = f"""# ITERATIVE GENE SET ENRICHMENT NETWORK ANALYSIS
# JSON Format for AI Analysis

## ANALYSIS INSTRUCTIONS
You are analyzing an iterative gene set enrichment network. The data below is provided in JSON format for easy parsing.

Please provide:
1. **Executive Summary**: Key biological themes and experimental context
2. **Hub Gene Analysis**: Roles of highly connected genes in biological processes  
3. **Pathway Relationships**: How biological terms relate to each other
4. **Iteration Patterns**: Biological significance of the iteration sequence
5. **Biological Hypothesis**: Potential experimental conditions and follow-up studies
6. **Clinical Relevance**: Therapeutic targets and biomarkers

Focus on biological interpretation and actionable insights.

## NETWORK DATA (JSON Format)
```json
{json.dumps(network_data, indent=2)}
```

## KEY INSIGHTS TO LOOK FOR:
- Genes with high centrality scores are likely key regulators or biomarkers
- Earlier iterations represent the most statistically significant enrichments
- Connected biological terms may represent related pathways or processes
- The iteration sequence reveals the hierarchical importance of biological processes
"""
    
    return output


def generate_regular_enrichment_json_analysis(enrichment_results: dict) -> str:
    """
    Generate a lightweight JSON format of regular enrichment results for AI analysis.
    This provides only essential data: Library, Term, Rank, Genes.
    
    :param enrichment_results: Dictionary of enrichment results from user-selected libraries
    :return: JSON format for AI analysis
    """
    import json
    
    # Extract data from enrichment results
    results = []
    
    for library_name, enrichment_obj in enrichment_results.items():
        # Get the DataFrame with all results
        df = enrichment_obj.to_dataframe()
        
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
    
    # Create the JSON structure
    analysis_data = {
        "analysis_type": "over_representation_analysis",
        "results": results
    }
    
    # Create the output with instructions
    output = f"""# OVER-REPRESENTATION ANALYSIS (ORA)
# Ranked List of Enriched Biological Terms

## ANALYSIS CONTEXT
You are a computational biologist analyzing over-representation analysis (ORA) results. This analysis identifies biological pathways, processes, and functions that are significantly enriched in a given gene set compared to a background gene set. The results are provided as a ranked list of enriched terms, where each term contains a list of genes that are over-represented in that biological process.

**DATA STRUCTURE:**
- **Ranked List**: Terms are ranked by statistical significance (lower rank = more significant)
- **Gene Lists**: Each term contains the specific genes that are enriched in that biological process
- **Library Sources**: Different databases provide complementary biological context
- **Table Format**: Results are organized by ranked terms per library source

## ENRICHMENT DATA (JSON Format)
```json
{json.dumps(analysis_data, indent=2)}
```

## LIBRARY SOURCES AND THEIR BIOLOGICAL CONTEXT:
- **H: Hallmark gene sets**: Curated gene sets representing well-defined biological states or processes
- **C2: BioCarta**: Canonical pathways from BioCarta database
- **C2: KEGG MEDICUS**: Metabolic and signaling pathways from KEGG
- **C2: Pathway Interaction Database**: Curated human signaling pathways
- **C2: Reactome Pathways**: Expert-curated biological pathways
- **C2: WikiPathways**: Community-curated biological pathways
- **C5: Gene Ontology: Biological Process**: Biological processes from Gene Ontology
- **C5: Gene Ontology: Cellular Component**: Cellular components from Gene Ontology
- **C5: Gene Ontology: Molecular Function**: Molecular functions from Gene Ontology
- **C5: Human Phenotype Ontology**: Human phenotypes and diseases
- **Protein Interaction**: Protein-protein interaction networks

## ANALYSIS REQUEST:

Please analyze this ranked list of enriched terms and provide:

1. **Key Biological Insights:**
   - What are the most significant biological processes/pathways identified?
   - Which genes appear may be central to the biological response?

2. **Ranked List Analysis:**
   - How do the terms relate to each other biologically?
   - Are there patterns in the gene composition across different terms and processes?

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
- Focus on biological interpretation, not just ranking
- Consider the functional relationships between genes within and across terms
- Consider that some genes may appear in multiple terms due to term redundancy
- Consider the broader biological context and literature
- Different library sources may provide complementary biological insights
- Build a coherent story that explains as many terms as possible rather than treating them in isolation
- The goal is to find a unifying biological hypothesis that connects the enriched terms into a meaningful narrative
"""
    
    return output
