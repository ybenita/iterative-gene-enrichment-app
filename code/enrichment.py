import json
import logging
import multiprocessing as mp
from datetime import datetime
from typing import Any, Dict, List, Tuple, Set

import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact, hypergeom
from statsmodels.stats.multitest import multipletests

from background_gene_set import BackgroundGeneSet
from gene_set import GeneSet
from gene_set_library import GeneSetLibrary

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def format_term_name(term_name: str) -> str:
    """
    Format term name by replacing underscores with spaces and adding colon after library shortcut.
    
    Args:
        term_name: Original term name (e.g., "GOBP_CELLULAR_RESPONSE_TO_STRESS")
        
    Returns:
        Formatted term name (e.g., "GOBP: CELLULAR RESPONSE TO STRESS")
    """
    # Replace underscores with spaces
    formatted = term_name.replace("_", " ")
    
    # Check if the term already contains a colon (indicating library name might be prepended)
    if ":" in formatted:
        # If it already has a colon, just return the formatted version
        return formatted
    
    # Split by space to get the first word (library shortcut)
    words = formatted.split()
    if len(words) > 1:
        # Add colon after the first word
        library_shortcut = words[0]
        rest_of_term = " ".join(words[1:])
        return f"{library_shortcut}: {rest_of_term}"
    else:
        # If only one word, just return it
        return formatted


def compute_pvalue(
    args: Tuple[GeneSet, BackgroundGeneSet, dict, str, int, Set[str]],
) -> Tuple[str, str, str, List[str], float]:
    """
    Computes the p-value for a given term using Fisher's exact test.
    This function is intended to be used with multiprocessing.Pool.map(),
    which requires functions to take a single argument. Therefore, the
    inputs are passed as a single tuple.

    Args:
        args: A tuple containing the following elements:
            - gene_set (GeneSet): The input gene set
            - background_gene_set (BackgroundGeneSet): The background gene set
            - term (dict): A dictionary representing a term in the gene set library
            - p-value method (str): The name of the method to calculate p-value
            - library_background_size (int): The size of the intersection between background and library genes
            - filtered_unique_genes (Set[str]): The unique genes from filtered terms

    Returns:
        A tuple containing the following elements:
            - term name (str): The name of the term
            - overlap size (int): The size of the overlap between the gene set and the term
            - term description (str): The description of the term
            - overlap genes (list): A list of genes in the overlap
            - p_value (float): The p-value computed by Fisher's exact test
    """
    gene_set, background_gene_set, term, p_value_method_name, library_background_size, filtered_unique_genes = args
    term_genes = set(term["genes"])
    n_term_genes = len(term_genes)
    
    # Intersect input gene set with library-specific background for statistical consistency
    # This ensures we only test genes that exist in the library-specific background
    filtered_gene_set = gene_set.genes & filtered_unique_genes
    filtered_gene_set_size = len(filtered_gene_set)
    
    # Calculate overlap using the filtered gene set
    overlap = filtered_gene_set & term_genes
    n_overlap = len(overlap)

    # Build contingency table for Fisher's exact test using library-specific background and filtered gene set
    contingency_table = [
        [n_overlap, n_term_genes - n_overlap],
        [
            filtered_gene_set_size - n_overlap,
            library_background_size - n_term_genes - filtered_gene_set_size + n_overlap,
        ],
    ]

    if p_value_method_name == "Fisher's Exact Test":
        _, p_value = fisher_exact(contingency_table)
    elif p_value_method_name == "Chi-squared Test":
        chi2, p_value, _, _ = chi2_contingency(contingency_table)
    elif p_value_method_name == "Hypergeometric Test":
        p_value = hypergeom.sf(
            n_overlap - 1, library_background_size, n_term_genes, filtered_gene_set_size
        )
    else:
        logger.error(f"Unsupported p_value_method: {p_value_method_name}")
        raise ValueError(f"Unsupported p_value_method: {p_value_method_name}")
    return (
        term["name"],
        f'{len(overlap)}/{len(term["genes"])}',
        term["description"],
        sorted(list(overlap)),
        p_value,
    )


class Enrichment:
    """
    Class for gene set enrichment analysis results.
    """

    def __init__(
        self,
        gene_set: GeneSet,
        gene_set_library: GeneSetLibrary,
        background_gene_set: BackgroundGeneSet,
        min_term_size: int = 10,
        max_term_size: int = 600,
        p_value_method_name="Fisher's Exact Test",
        name: str = None,
    ):
        """
        Initialize the class with gene set, gene set library, and background gene set.

        Args:
            gene_set: Input gene set
            gene_set_library: Gene set library
            background_gene_set: Background gene set
        """
        self.gene_set = gene_set
        self.gene_set_library = gene_set_library
        self.min_term_size = min_term_size
        self.max_term_size = max_term_size
        self.background_gene_set = background_gene_set
        self.p_value_method_name = p_value_method_name
        self.name = (
            name
            if name
            else f"{gene_set.name}_{gene_set_library.name}_{min_term_size}-{max_term_size}_{background_gene_set.name}_{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
        )
        self._results: List[Dict[str, Any]] = self._compute_enrichment()

    @property
    def results(self) -> List[Dict[str, Any]]:
        """
        Getter for _results.

        Returns:
            A list containing dictionaries of enrichment results
        """
        return self._results

    @results.setter
    def results(self, value: List[Dict[str, Any]]) -> None:
        """
        Setter for _results.

        Args:
            value: A list containing dictionaries of enrichment results
        """
        self._results = value

    def _compute_enrichment(self) -> List[Dict[str, Any]]:
        """
        Computes gene set enrichment analysis.

        Returns:
            A list containing dictionaries of enrichment results
        """
        results = []
        logger.info(f"Calculating p-values for {self.gene_set_library.name}")
        cpu_count = mp.cpu_count() - 2
        parallel_results = []  # Initialize outside try block
        
        with mp.Pool(cpu_count) as pool:
            logger.info(f"Initializing the MP pool with {cpu_count} CPUs")
            try:
                # Filter terms by size range first
                filtered_terms = [term for term in self.gene_set_library.library 
                                 if self.min_term_size <= term["size"] <= self.max_term_size]
                
                # Calculate unique genes from filtered terms only
                filtered_unique_genes = set()
                for term in filtered_terms:
                    filtered_unique_genes.update(term["genes"])
                
                # Calculate the size of the intersection between the background gene set and the filtered library's unique genes
                library_background_size = len(self.background_gene_set.genes & filtered_unique_genes)
                logger.info(f"Library-specific background size: {library_background_size} genes (intersection of {self.background_gene_set.size} background genes and {len(filtered_unique_genes)} filtered library genes from {len(filtered_terms)} terms within size range [{self.min_term_size}, {self.max_term_size}])")
                
                # Log the impact of filtering the input gene set
                original_gene_set_size = self.gene_set.size
                filtered_gene_set_size = len(self.gene_set.genes & filtered_unique_genes)
                if original_gene_set_size != filtered_gene_set_size:
                    logger.info(f"Input gene set filtered: {original_gene_set_size} → {filtered_gene_set_size} genes (intersected with library-specific background)")
                else:
                    logger.info(f"Input gene set size: {original_gene_set_size} genes (all genes present in library-specific background)")
                
                parallel_results = pool.map(
                    compute_pvalue,
                    [
                        (
                            self.gene_set,
                            self.background_gene_set,
                            term,
                            self.p_value_method_name,
                            library_background_size,
                            filtered_unique_genes,
                        )
                        for term in filtered_terms
                    ],
                )
            except Exception as e:
                logging.exception("An error occurred: %s", e)
                return results  # Return empty results if computation failed
            finally:
                pool.close()
                pool.join()
                logger.info(f"Releasing {cpu_count} CPUs from the MP pool")

        # Check if we have results to process
        if not parallel_results:
            logger.warning(f"No results obtained for {self.gene_set_library.name}")
            logger.info(f"Library has {len(self.gene_set_library.library)} total terms")
            logger.info(f"Terms within size range [{self.min_term_size}, {self.max_term_size}]: {len(filtered_terms)}")
            logger.info(f"Input gene set size: {self.gene_set.size}")
            logger.info(f"Background gene set size: {self.background_gene_set.size}")
            return results

        # Separate results and p_values for convenience
        p_values = [result[-1] for result in parallel_results]
        # Adjust p-values for multiple testing
        _, p_values_adjusted, _, _ = multipletests(p_values, method="fdr_bh")
        # Rank terms based on their p-values
        ranked_terms = sorted(
            list(enumerate(parallel_results)), key=lambda x: p_values[x[0]]
        )

        # Format results into a sorted list
        for i, result in ranked_terms:
            term_name, overlap_size, term_description, overlap_genes, _ = result
            results.append(
                {
                    "term": term_name,
                    "rank": i + 1,
                    "description": term_description,
                    "overlap": overlap_genes,
                    "overlap_size": overlap_size,
                    "p-value": p_values[i],
                    "fdr": p_values_adjusted[i],
                }
            )
        return results

    def to_dataframe(self):
        """Return the enrichment results as a pandas dataframe."""
        import math
        
        df = pd.DataFrame(
            {
                "Library": [self.gene_set_library.name for _ in self.results],
                "Rank": [result["rank"] for result in self.results],
                "Term": [format_term_name(result["term"]) for result in self.results],
                "Description": [result.get("description", "") for result in self.results],
                "Overlap size": [result["overlap_size"] for result in self.results],
                "Genes": [", ".join(result["overlap"]) for result in self.results],
                "p-value": [result["p-value"] for result in self.results],
                "-log(p-value)": [-math.log10(result["p-value"]) if result["p-value"] > 0 else 0 for result in self.results],
                "FDR": [result["fdr"] for result in self.results],
            }
        )
        # Reorder columns to put Library first, then add -log(p-value) after p-value
        column_order = ["Library", "Rank", "Term", "Description", "Overlap size", "p-value", "-log(p-value)", "FDR", "Genes"]
        return df[column_order]

    def to_json(self):
        """Return the enrichment results as a JSON string."""
        return json.dumps(self.results, indent=4, separators=(",", ": "))

    def to_html(self):
        """Return the enrichment results as an HTML page."""
        return self.to_dataframe().to_html()

    def to_tsv(self):
        """Return the enrichment results as a TSV spreadsheet."""
        return self.to_dataframe().to_csv(sep="\t")

    def to_snapshot(self) -> Dict:
        """Return the snapshot of input parameters and the enrichment results as a JSON string."""
        # Calculate library-specific background size for the snapshot
        # Filter terms by size range first
        filtered_terms = [term for term in self.gene_set_library.library 
                         if self.min_term_size <= term["size"] <= self.max_term_size]
        
        # Calculate unique genes from filtered terms only
        filtered_unique_genes = set()
        for term in filtered_terms:
            filtered_unique_genes.update(term["genes"])
        
        # Calculate the size of the intersection between the background gene set and the filtered library's unique genes
        library_background_size = len(self.background_gene_set.genes & filtered_unique_genes)
        return {
            "input_gene_set": list(self.gene_set.genes),
            "background": self.background_gene_set.name,
            "background_size": self.background_gene_set.size,
            "library_background_size": library_background_size,
            "library_size": self.gene_set_library.size,
            self.gene_set_library.name: self.results,
        }
