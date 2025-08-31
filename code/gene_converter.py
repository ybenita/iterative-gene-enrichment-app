import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

logger = logging.getLogger(__name__)


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
                          If None, will look for it in the project root.
            gene_history_path: Path to the gene_history file.
                             If None, will look for it in the recent_release directory.
        """
        if gene_info_path is None:
            # Look for the file in the project root
            project_root = Path(__file__).resolve().parent.parent
            gene_info_path = project_root / "Homo_sapiens.gene_info"
        
        self.gene_info_path = Path(gene_info_path)
        
        if gene_history_path is None:
            # Look for the file in the recent_release directory
            project_root = Path(__file__).resolve().parent.parent
            gene_history_path = project_root / "data" / "recent_release" / "gene_history"
        
        self.gene_history_path = Path(gene_history_path)
        
        self.entrez_to_symbol: Dict[str, str] = {}
        self.symbol_to_entrez: Dict[str, str] = {}
        self.synonyms_to_symbol: Dict[str, str] = {}
        self.old_to_current_symbol: Dict[str, str] = {}
        
        # Track conversions for this session
        self.conversions: List[str] = []
        
        if self.gene_info_path.exists():
            self._load_gene_info()
        else:
            logger.warning(f"Gene info file not found at {self.gene_info_path}")
        
        if self.gene_history_path.exists():
            self._load_gene_history()
        else:
            logger.warning(f"Gene history file not found at {self.gene_history_path}")
    
    def _load_gene_info(self) -> None:
        """Load gene information from the NCBI gene info file."""
        logger.info(f"Loading gene info from {self.gene_info_path}")
        
        with open(self.gene_info_path, 'r', encoding='utf-8') as f:
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
            with open(self.gene_history_path, 'r', encoding='utf-8') as f:
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