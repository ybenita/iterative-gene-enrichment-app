from pathlib import Path
from typing import Set, Optional

import streamlit as st

from gene_converter import GeneConverter


@st.cache_data(show_spinner="Loading background gene set...")
def _cached_load_background(background_file_path: str, input_format: str = "symbols", skip_validation: bool = False) -> Set[str]:
    """Load and cache background gene set from file. Returns set of gene symbols."""
    import logging
    logger = logging.getLogger(__name__)
    
    with open(background_file_path, "r") as f:
        raw_lines = [line.strip() for line in f.readlines() if line.strip()]
    
    if input_format == "entrez_ids":
        converter = GeneConverter()
        converted_symbols = []
        unrecognized_entrez = []
        for line in raw_lines:
            gene_id = line.strip()
            if not gene_id:
                continue
            symbol = converter.get_symbol(gene_id)
            if symbol:
                converted_symbols.append(symbol)
            else:
                unrecognized_entrez.append(gene_id)
        if unrecognized_entrez:
            logger.warning(f"{len(unrecognized_entrez)} Entrez IDs not found")
        return set(converted_symbols)
    else:
        if skip_validation:
            return set(line.strip().upper() for line in raw_lines if line.strip())
        
        converter = GeneConverter()
        stats = converter.get_stats()
        has_gene_data = stats.get('symbol_mappings', 0) > 0
        
        if not has_gene_data:
            logger.warning("GeneConverter has no loaded data. Accepting all symbols without validation.")
            return set(line.strip().upper() for line in raw_lines if line.strip())
        
        valid_symbols = []
        invalid_symbols = []
        for line in raw_lines:
            gene_id = line.strip()
            if not gene_id:
                continue
            if converter.is_symbol(gene_id):
                valid_symbols.append(gene_id)
            else:
                invalid_symbols.append(gene_id)
        
        if invalid_symbols:
            logger.warning(f"{len(invalid_symbols)} symbols not found in database")
        if len(valid_symbols) == 0 and len(raw_lines) > 0:
            logger.error(f"All {len(raw_lines)} symbols were rejected by GeneConverter")
        
        return set(valid_symbols)


class BackgroundGeneSet:
    """
    A class to store a list of genes and their size.
    """

    def __init__(
        self, 
        background_file_path: str, 
        name: str = "", 
        organism: str = "Homo Sapiens",
        input_format: str = "symbols",
        skip_validation: bool = False
    ) -> None:
        """
        Initialize BackgroundGeneList object with a list of genes.

        Args:
            background_file_path: Path to the background file
            name: Name for the background gene list
            organism: Organism name
            input_format: Either 'symbols' or 'entrez_ids'
            skip_validation: If True, skip gene symbol validation (faster, use for trusted sources)
        """
        self.genes: Set[str] = self._load_from_file(background_file_path, input_format, skip_validation)
        self.size: int = len(self.genes)
        self.name = name if name else Path(background_file_path).stem
        self.organism = organism
        self.input_format = input_format

    def _load_from_file(self, background_file_path: str, input_format: str = "symbols", skip_validation: bool = False) -> Set[str]:
        """
        Load background genes from a file (cached across Streamlit reruns).
        
        Args:
            background_file_path: Path to the background file
            input_format: Either 'symbols' or 'entrez_ids'
            skip_validation: If True, skip gene symbol validation

        Returns:
            Set of gene symbols representing the background
        """
        return _cached_load_background(background_file_path, input_format, skip_validation)

    def has_gene(self, gene: str) -> bool:
        """
        Check if the given gene is present in the BackgroundGenes.

        Args:
            gene: A gene name.

        Returns:
            True if the gene is present, False otherwise.
        """
        return gene in self.genes
