from pathlib import Path
from typing import Set, Optional

from gene_converter import GeneConverter


class BackgroundGeneSet:
    """
    A class to store a list of genes and their size.
    """

    def __init__(
        self, 
        background_file_path: str, 
        name: str = "", 
        organism: str = "Homo Sapiens",
        input_format: str = "symbols"
    ) -> None:
        """
        Initialize BackgroundGeneList object with a list of genes.

        Args:
            background_file_path: Path to the background file
            name: Name for the background gene list
            organism: Organism name
            input_format: Either 'symbols' or 'entrez_ids'
        """
        self.genes: Set[str] = self._load_from_file(background_file_path, input_format)
        self.size: int = len(self.genes)
        self.name = name if name else Path(background_file_path).stem
        self.organism = organism
        self.input_format = input_format

    def _load_from_file(self, background_file_path: str, input_format: str = "symbols") -> Set[str]:
        """
        Load background genes from a file with optional Entrez ID conversion
        
        Args:
            background_file_path: Path to the background file
            input_format: Either 'symbols' or 'entrez_ids'

        Returns:
            Set of gene symbols representing the background
        """
        # Read raw lines from file
        with open(background_file_path, "r") as f:
            raw_lines = [line.strip() for line in f.readlines() if line.strip()]
        
        if input_format == "entrez_ids":
            # Convert Entrez IDs to symbols
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
            
            # Log conversion results
            if unrecognized_entrez:
                print(f"Warning: {len(unrecognized_entrez)} Entrez IDs not found in database: {unrecognized_entrez[:10]}{'...' if len(unrecognized_entrez) > 10 else ''}")
            
            return set(converted_symbols)
        else:
            # Validate symbols
            converter = GeneConverter()
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
            
            # Log validation results
            if invalid_symbols:
                print(f"Warning: {len(invalid_symbols)} symbols not found in database: {invalid_symbols[:10]}{'...' if len(invalid_symbols) > 10 else ''}")
            
            return set(valid_symbols)

    def has_gene(self, gene: str) -> bool:
        """
        Check if the given gene is present in the BackgroundGenes.

        Args:
            gene: A gene name.

        Returns:
            True if the gene is present, False otherwise.
        """
        return gene in self.genes
