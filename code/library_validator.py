#!/usr/bin/env python3
"""
Library validator module for ensuring consistency between gene_info database and MSigDB libraries.
This module validates gene symbols in GMT files and updates them to use consistent symbols.
"""

import logging
import os
import shutil
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Set, Tuple

from code.gene_converter import GeneConverter

logger = logging.getLogger(__name__)


class LibraryValidator:
    """
    Validates and updates gene symbols in GMT files to ensure consistency with gene_info database.
    """
    
    def __init__(self, gene_info_path: str = None, gene_history_path: str = None):
        """
        Initialize the library validator.
        
        Args:
            gene_info_path: Path to the Homo_sapiens.gene_info file
            gene_history_path: Path to the gene_history file
        """
        self.converter = GeneConverter(gene_info_path, gene_history_path)
        self.stats = {
            'total_genes': 0,
            'validated_genes': 0,
            'unchanged_genes': 0,
            'invalid_genes': 0,
            'terms_processed': 0,
            'files_processed': 0
        }
    
    def validate_gmt_file(self, gmt_file_path: str, output_path: str = None, 
                         create_backup: bool = True) -> Dict[str, int]:
        """
        Validate and update gene symbols in a single GMT file.
        
        Args:
            gmt_file_path: Path to the input GMT file
            output_path: Path for the output file (if None, overwrites input)
            create_backup: Whether to create a backup of the original file
            
        Returns:
            Dictionary with validation statistics
        """
        gmt_path = Path(gmt_file_path)
        if not gmt_path.exists():
            raise FileNotFoundError(f"GMT file not found: {gmt_file_path}")
        
        # Create backup if requested
        if create_backup and output_path is None:
            backup_path = gmt_path.with_suffix(f'.backup_{datetime.now().strftime("%Y%m%d_%H%M%S")}')
            shutil.copy2(gmt_path, backup_path)
            logger.info(f"Created backup: {backup_path}")
        
        # Set output path
        if output_path is None:
            output_path = gmt_file_path
        
        # Process the file
        file_stats = self._process_gmt_file(gmt_path, output_path)
        
        # Update global stats
        for key in file_stats:
            self.stats[key] += file_stats[key]
        
        self.stats['files_processed'] += 1
        
        logger.info(f"Processed {gmt_path.name}: {file_stats['terms_processed']} terms, "
                   f"{file_stats['validated_genes']} genes validated, "
                   f"{file_stats['invalid_genes']} genes unchanged")
        
        return file_stats
    
    def validate_library_directory(self, library_dir: str, output_dir: str = None,
                                  file_pattern: str = "*.gmt", create_backup: bool = True) -> Dict[str, int]:
        """
        Validate all GMT files in a directory.
        
        Args:
            library_dir: Directory containing GMT files
            output_dir: Output directory (if None, overwrites original files)
            file_pattern: Pattern to match GMT files
            create_backup: Whether to create backups of original files
            
        Returns:
            Dictionary with overall validation statistics
        """
        library_path = Path(library_dir)
        if not library_path.exists():
            raise FileNotFoundError(f"Library directory not found: {library_dir}")
        
        if output_dir:
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)
        
        # Find all GMT files
        gmt_files = list(library_path.glob(file_pattern))
        if not gmt_files:
            logger.warning(f"No GMT files found in {library_dir} matching pattern '{file_pattern}'")
            return self.stats
        
        logger.info(f"Found {len(gmt_files)} GMT files to process")
        
        # Process each file
        for gmt_file in gmt_files:
            if output_dir:
                output_file = output_path / gmt_file.name
            else:
                output_file = None
            
            try:
                self.validate_gmt_file(str(gmt_file), str(output_file) if output_file else None, create_backup)
            except Exception as e:
                logger.error(f"Error processing {gmt_file.name}: {e}")
                continue
        
        logger.info(f"Completed processing {self.stats['files_processed']} files")
        return self.stats
    
    def _process_gmt_file(self, input_path: Path, output_path: str) -> Dict[str, int]:
        """
        Process a single GMT file and write the validated version.
        
        Args:
            input_path: Path to input GMT file
            output_path: Path to output GMT file
            
        Returns:
            Dictionary with file-specific statistics
        """
        file_stats = {
            'total_genes': 0,
            'validated_genes': 0,
            'unchanged_genes': 0,
            'invalid_genes': 0,
            'terms_processed': 0
        }
        
        validated_lines = []
        
        with open(input_path, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                
                # Parse GMT line: term_name\tdescription\tgene1\tgene2\t...
                parts = line.split('\t')
                if len(parts) < 3:
                    logger.warning(f"Line {line_num}: Invalid GMT format, skipping")
                    validated_lines.append(line)
                    continue
                
                term_name = parts[0]
                description = parts[1]
                genes = parts[2:]
                
                # Validate genes
                validated_genes = []
                for gene in genes:
                    file_stats['total_genes'] += 1
                    
                    # Try to validate the gene symbol
                    validated_symbol = self.converter.validate_and_map_symbol(gene)
                    
                    if validated_symbol:
                        if validated_symbol != gene.upper():
                            # Gene was updated
                            validated_genes.append(validated_symbol)
                            file_stats['validated_genes'] += 1
                            logger.debug(f"Updated gene: {gene} -> {validated_symbol}")
                        else:
                            # Gene was already correct (just case change)
                            validated_genes.append(validated_symbol)
                            file_stats['unchanged_genes'] += 1
                    else:
                        # Gene could not be validated, keep original
                        validated_genes.append(gene)
                        file_stats['invalid_genes'] += 1
                        logger.debug(f"Could not validate gene: {gene}")
                
                # Reconstruct the line
                validated_line = f"{term_name}\t{description}\t{'\t'.join(validated_genes)}"
                validated_lines.append(validated_line)
                file_stats['terms_processed'] += 1
        
        # Write the validated file
        with open(output_path, 'w', encoding='utf-8') as f:
            for line in validated_lines:
                f.write(line + '\n')
        
        return file_stats
    
    def get_validation_report(self) -> str:
        """
        Generate a detailed validation report.
        
        Returns:
            Formatted report string
        """
        report = []
        report.append("=" * 60)
        report.append("LIBRARY VALIDATION REPORT")
        report.append("=" * 60)
        report.append(f"Files processed: {self.stats['files_processed']}")
        report.append(f"Terms processed: {self.stats['terms_processed']}")
        report.append(f"Total genes: {self.stats['total_genes']}")
        report.append(f"Validated genes: {self.stats['validated_genes']} ({self.stats['validated_genes']/max(self.stats['total_genes'], 1)*100:.1f}%)")
        report.append(f"Unchanged genes: {self.stats['unchanged_genes']} ({self.stats['unchanged_genes']/max(self.stats['total_genes'], 1)*100:.1f}%)")
        report.append(f"Invalid genes: {self.stats['invalid_genes']} ({self.stats['invalid_genes']/max(self.stats['total_genes'], 1)*100:.1f}%)")
        report.append("=" * 60)
        
        return '\n'.join(report)
    
    def reset_stats(self):
        """Reset validation statistics."""
        self.stats = {
            'total_genes': 0,
            'validated_genes': 0,
            'unchanged_genes': 0,
            'invalid_genes': 0,
            'terms_processed': 0,
            'files_processed': 0
        }


def main():
    """
    Command-line interface for library validation.
    """
    import argparse
    
    parser = argparse.ArgumentParser(description="Validate and update gene symbols in GMT files")
    parser.add_argument("input", help="Input GMT file or directory")
    parser.add_argument("-o", "--output", help="Output file or directory")
    parser.add_argument("--no-backup", action="store_true", help="Don't create backup files")
    parser.add_argument("--pattern", default="*.gmt", help="File pattern for directory processing")
    parser.add_argument("--gene-info", help="Path to Homo_sapiens.gene_info file")
    parser.add_argument("--gene-history", help="Path to gene_history file")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    
    args = parser.parse_args()
    
    # Set up logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=log_level, format='%(asctime)s - %(levelname)s - %(message)s')
    
    # Initialize validator
    validator = LibraryValidator(args.gene_info, args.gene_history)
    
    try:
        input_path = Path(args.input)
        
        if input_path.is_file():
            # Process single file
            stats = validator.validate_gmt_file(
                str(input_path), 
                args.output, 
                create_backup=not args.no_backup
            )
        elif input_path.is_dir():
            # Process directory
            stats = validator.validate_library_directory(
                str(input_path), 
                args.output, 
                args.pattern
            )
        else:
            logger.error(f"Input path does not exist: {args.input}")
            return 1
        
        # Print report
        print(validator.get_validation_report())
        
        return 0
        
    except Exception as e:
        logger.error(f"Validation failed: {e}")
        return 1


if __name__ == "__main__":
    exit(main())
