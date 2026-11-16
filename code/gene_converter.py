import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import streamlit as st

logger = logging.getLogger(__name__)


@st.cache_resource(show_spinner="Loading gene info database...")
def _cached_load_gene_info(gene_info_path: str) -> Tuple[Dict[str, str], Dict[str, str], Dict[str, str]]:
    """Load and cache gene info parsing. Returns (entrez_to_symbol, symbol_to_entrez, synonyms_to_symbol)."""
    entrez_to_symbol: Dict[str, str] = {}
    symbol_to_entrez: Dict[str, str] = {}
    synonyms_to_symbol: Dict[str, str] = {}
    
    with GeneConverter._open(gene_info_path) as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#'):
                continue
            try:
                fields = line.strip().split('\t')
                if len(fields) < 3:
                    continue
                tax_id = fields[0]
                if tax_id != '9606':
                    continue
                entrez_id = fields[1]
                symbol = fields[2]
                if not entrez_id or not symbol or symbol == '-':
                    continue
                entrez_to_symbol[entrez_id] = symbol
                symbol_to_entrez[symbol.upper()] = entrez_id
                if len(fields) > 4 and fields[4] != '-':
                    synonyms = fields[4].split('|')
                    for synonym in synonyms:
                        if synonym and synonym != '-':
                            synonyms_to_symbol[synonym.upper()] = symbol
            except Exception as e:
                logger.warning(f"Error parsing line {line_num}: {e}")
                continue
    
    logger.info(f"Loaded {len(entrez_to_symbol)} human Entrez ID mappings")
    return entrez_to_symbol, symbol_to_entrez, synonyms_to_symbol


@st.cache_resource(show_spinner="Loading gene history...")
def _cached_load_gene_history(gene_history_path: str) -> Dict[str, str]:
    """Load and cache gene history parsing. Returns old_to_current_symbol mapping."""
    old_to_current_symbol: Dict[str, str] = {}
    
    with GeneConverter._open(gene_history_path) as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#'):
                continue
            try:
                fields = line.strip().split('\t')
                if len(fields) < 4:
                    continue
                tax_id = fields[0]
                if tax_id != '9606':
                    continue
                old_symbol = fields[2]
                current_symbol = fields[3]
                if not old_symbol or not current_symbol or old_symbol == '-' or current_symbol == '-':
                    continue
                if old_symbol.upper() not in old_to_current_symbol:
                    old_to_current_symbol[old_symbol.upper()] = current_symbol.upper()
                elif old_to_current_symbol[old_symbol.upper()] != current_symbol.upper():
                    del old_to_current_symbol[old_symbol.upper()]
            except Exception as e:
                logger.warning(f"Error parsing gene_history line {line_num}: {e}")
                continue
    
    logger.info(f"Loaded {len(old_to_current_symbol)} unique old-to-current symbol mappings")
    return old_to_current_symbol


class GeneConverter:
    """
    Converter between Entrez IDs and gene symbols using NCBI gene info data.
    Enhanced with gene history support for mapping old symbols to current symbols.
    """
    
    def __init__(self, gene_info_path: Optional[str] = None, gene_history_path: Optional[str] = None):
        """
        Initialize the gene converter.
        
        Args:
            gene_info_path: Path to the Homo_sapiens.gene_info file. 
                          If None, will look for it in the data/recent_release directory.
            gene_history_path: Path to the gene_history file.
                             If None, will look for it in the data/recent_release directory.
        """
        if gene_info_path is None:
            gene_info_path = self._find_data_file(["Homo_sapiens.gene_info", "Homo_sapiens.gene_info.gz"])
            if gene_info_path is None:
                raise FileNotFoundError(
                    "NCBI gene info not found. Expected data/recent_release/Homo_sapiens.gene_info(.gz) in the "
                    "repository. Gene symbols cannot be validated without it, and results would depend on the "
                    "working directory.")
        self.gene_info_path = Path(gene_info_path)

        if gene_history_path is None:
            gene_history_path = self._find_data_file(["gene_history", "gene_history_9606.gz", "gene_history.gz"])
            if gene_history_path is None:
                logger.warning("NCBI gene history not found in data/recent_release; outdated symbols will not be mapped")
                gene_history_path = "gene_history_not_found"
        self.gene_history_path = Path(gene_history_path)
        
        # Track conversions for this session
        self.conversions: List[str] = []
        
        # Load from cached functions
        if self.gene_info_path.exists():
            self.entrez_to_symbol, self.symbol_to_entrez, self.synonyms_to_symbol = (
                _cached_load_gene_info(str(self.gene_info_path))
            )
        else:
            logger.warning(f"Gene info file not found at {self.gene_info_path}")
            self.entrez_to_symbol: Dict[str, str] = {}
            self.symbol_to_entrez: Dict[str, str] = {}
            self.synonyms_to_symbol: Dict[str, str] = {}
        
        if self.gene_history_path.exists():
            self.old_to_current_symbol = _cached_load_gene_history(str(self.gene_history_path))
        else:
            logger.warning(f"Gene history file not found at {self.gene_history_path}")
            self.old_to_current_symbol: Dict[str, str] = {}
    
    @staticmethod
    def _find_data_file(names):
        """First existing file among ``names`` in data/recent_release (repository first, then the
        current directory and its parent, then /data/recent_release)."""
        here = Path(__file__).resolve().parent.parent
        for d in (here / "data" / "recent_release", Path.cwd() / "data" / "recent_release",
                  Path.cwd().parent / "data" / "recent_release", Path("/data") / "recent_release"):
            for n in names:
                if (d / n).exists():
                    return d / n
        return None

    @staticmethod
    def _open(path):
        import gzip
        path = Path(path)
        if path.suffix == ".gz":
            return gzip.open(path, "rt", encoding="utf-8")
        return open(path, "r", encoding="utf-8")

    def _load_gene_info(self) -> None:
        """Load gene information from the NCBI gene info file."""
        logger.info(f"Loading gene info from {self.gene_info_path}")
        
        with self._open(self.gene_info_path) as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#'):
                    continue
                
                try:
                    fields = line.strip().split('\t')
                    if len(fields) < 3:
                        continue
                    
                    tax_id = fields[0]
                    if tax_id != '9606':  # Only human genes
                        continue
                    
                    entrez_id = fields[1]
                    symbol = fields[2]
                    
                    # Skip entries with missing or invalid data
                    if not entrez_id or not symbol or symbol == '-':
                        continue
                    
                    # Store primary mappings
                    self.entrez_to_symbol[entrez_id] = symbol
                    self.symbol_to_entrez[symbol.upper()] = entrez_id
                    
                    # Store synonyms if available
                    if len(fields) > 4 and fields[4] != '-':
                        synonyms = fields[4].split('|')
                        for synonym in synonyms:
                            if synonym and synonym != '-':
                                self.synonyms_to_symbol[synonym.upper()] = symbol
                
                except Exception as e:
                    logger.warning(f"Error parsing line {line_num}: {e}")
                    continue
        
        logger.info(f"Loaded {len(self.entrez_to_symbol)} human Entrez ID mappings (taxonomy ID 9606)")
        logger.info(f"Loaded {len(self.synonyms_to_symbol)} synonym mappings")
        
        # Debug: Show a few examples
        if self.entrez_to_symbol:
            sample_entrez = list(self.entrez_to_symbol.keys())[:3]
            logger.info(f"Sample Entrez IDs: {sample_entrez}")
            logger.info(f"Sample symbols: {[self.entrez_to_symbol[e] for e in sample_entrez]}")
    
    def _load_gene_history(self) -> None:
        """Load gene history and create mapping from old symbols to current symbols."""
        logger.info(f"Loading gene history from {self.gene_history_path}")
        
        try:
            with self._open(self.gene_history_path) as f:
                for line_num, line in enumerate(f, 1):
                    if line.startswith('#'):
                        continue
                    
                    try:
                        fields = line.strip().split('\t')
                        if len(fields) < 4:
                            continue
                        
                        tax_id = fields[0]
                        if tax_id != '9606':  # Only human genes
                            continue
                        
                        old_symbol = fields[2]
                        current_symbol = fields[3]
                        
                        # Skip entries with missing or invalid data
                        if not old_symbol or not current_symbol or old_symbol == '-' or current_symbol == '-':
                            continue
                        
                        # Only add if there's a unique mapping (no conflicts)
                        if old_symbol.upper() not in self.old_to_current_symbol:
                            self.old_to_current_symbol[old_symbol.upper()] = current_symbol.upper()
                        elif self.old_to_current_symbol[old_symbol.upper()] != current_symbol.upper():
                            # Conflict found - remove this mapping
                            del self.old_to_current_symbol[old_symbol.upper()]
                            logger.debug(f"Removed conflicting mapping for {old_symbol}")
                    
                    except Exception as e:
                        logger.warning(f"Error parsing gene_history line {line_num}: {e}")
                        continue
            
            logger.info(f"Loaded {len(self.old_to_current_symbol)} unique old-to-current symbol mappings")
            
        except Exception as e:
            logger.error(f"Error loading gene_history: {e}")
    
    def validate_and_map_symbol(self, symbol: str) -> Optional[str]:
        """
        Validate a gene symbol and map it to current symbol if needed.
        
        Args:
            symbol: Gene symbol to validate
            
        Returns:
            Current gene symbol if valid, None otherwise
        """
        symbol_upper = symbol.upper()
        
        # First try direct match
        if symbol_upper in self.symbol_to_entrez:
            return symbol_upper
        
        # Try synonyms
        if symbol_upper in self.synonyms_to_symbol:
            mapped_symbol = self.synonyms_to_symbol[symbol_upper].upper()
            if mapped_symbol in self.symbol_to_entrez:
                # Track conversion
                if symbol_upper != mapped_symbol:
                    self.conversions.append(f"{symbol_upper}→{mapped_symbol}")
                return mapped_symbol
        
        # Try gene history mapping
        if symbol_upper in self.old_to_current_symbol:
            mapped_symbol = self.old_to_current_symbol[symbol_upper]
            if mapped_symbol in self.symbol_to_entrez:
                logger.debug(f"Mapped {symbol_upper} -> {mapped_symbol} using gene history")
                # Track conversion
                if symbol_upper != mapped_symbol:
                    self.conversions.append(f"{symbol_upper}→{mapped_symbol}")
                return mapped_symbol
        
        return None
    
    def convert_input(self, input_text: str) -> Tuple[List[str], List[str], List[str]]:
        """
        Convert input text to gene symbols, handling both Entrez IDs and gene symbols.
        Enhanced with gene history support for mapping old symbols.
        
        Args:
            input_text: Text containing gene identifiers (one per line)
            
        Returns:
            Tuple of (converted_symbols, unrecognized_entrez, unrecognized_symbols)
        """
        lines = [line.strip() for line in input_text.split('\n') if line.strip()]
        
        converted_symbols = []
        unrecognized_entrez = []
        unrecognized_symbols = []
        
        for line in lines:
            gene_id = line.strip()
            if not gene_id:
                continue
            
            # Check if it's an Entrez ID
            if self.is_entrez_id(gene_id):
                symbol = self.get_symbol(gene_id)
                if symbol:
                    converted_symbols.append(symbol.upper())
                else:
                    unrecognized_entrez.append(gene_id)
            else:
                # Try to validate and map the symbol
                mapped_symbol = self.validate_and_map_symbol(gene_id)
                if mapped_symbol:
                    converted_symbols.append(mapped_symbol)
                else:
                    unrecognized_symbols.append(gene_id)
        
        return converted_symbols, unrecognized_entrez, unrecognized_symbols
    
    def get_symbol(self, entrez_id: str) -> Optional[str]:
        """Get gene symbol for an Entrez ID."""
        return self.entrez_to_symbol.get(entrez_id)
    
    def get_entrez(self, symbol: str) -> Optional[str]:
        """Get Entrez ID for a gene symbol."""
        return self.symbol_to_entrez.get(symbol.upper())
    
    def is_entrez_id(self, gene_id: str) -> bool:
        """Check if the gene ID is an Entrez ID (numeric)."""
        return gene_id.isdigit()
    
    def is_symbol(self, gene_id: str) -> bool:
        """Check if the gene ID is a valid gene symbol."""
        return self.validate_and_map_symbol(gene_id) is not None
    
    def get_stats(self) -> Dict[str, int]:
        """Get statistics about loaded data."""
        return {
            'entrez_mappings': len(self.entrez_to_symbol),
            'symbol_mappings': len(self.symbol_to_entrez),
            'synonym_mappings': len(self.synonyms_to_symbol),
            'history_mappings': len(self.old_to_current_symbol)
        }
    
    def get_conversions(self) -> List[str]:
        """Get list of conversions made in this session."""
        return self.conversions.copy()
    
    def clear_conversions(self) -> None:
        """Clear the conversion tracking list."""
        self.conversions.clear() 