from pathlib import Path
from typing import List, Optional
import json
import logging
from datetime import datetime

import typer
import pandas as pd

from code.background_gene_set import BackgroundGeneSet
from code.enrichment import Enrichment
from code.gene_set import GeneSet
from code.gene_set_library import GeneSetLibrary
from code.iter_enrichment import IterativeEnrichment

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parent.parent
app = typer.Typer(
    help="Gene Set Enrichment Analysis CLI",
    add_completion=False,
    no_args_is_help=True,
    context_settings={"help_option_names": ["-h", "--help"]}
)


def run_enrichment(
    gene_sets: List[Path],
    background: Path,
    libraries: List[Path],
    gene_format: str,
    p_value_method: str,
    mode: str,
    p_threshold: float,
    min_overlap: int,
    min_term_size: int,
    max_term_size: int,
    max_iterations: int,
    output_dir: Path,
):
    """Run enrichment analysis for the given parameters."""
    
    # Ensure output directory exists
    output_dir.mkdir(exist_ok=True)
    
    # Load background gene set
    logger.info(f"Loading background gene set: {background}")
    background_gene_set = BackgroundGeneSet(str(background))
    
    # Load gene set libraries
    logger.info(f"Loading {len(libraries)} gene set libraries")
    gene_set_libraries = []
    for lib_path in libraries:
        lib_name = lib_path.stem
        library = GeneSetLibrary(str(lib_path), name=lib_name)
        # Filter terms by size (same as Streamlit app)
        filtered_terms = [
            t for t in library.library if min_term_size <= t["size"] <= max_term_size
        ]
        library.library = filtered_terms
        library.num_terms = len(filtered_terms)
        library.unique_genes = library.compute_unique_genes()
        library.size = len(library.unique_genes)
        gene_set_libraries.append(library)
        logger.info(f"Loaded library {lib_name}: {library.num_terms} terms")
    
    # Process each gene set
    for gene_set_path in gene_sets:
        logger.info(f"Processing gene set: {gene_set_path}")
        
        # Load gene set
        gene_set_name = gene_set_path.stem
        with open(gene_set_path, 'r') as f:
            gene_input = [line.strip() for line in f if line.strip()]
        
        # Handle gene format conversion
        if gene_format == "entrez_ids":
            # Convert Entrez IDs to gene symbols
            from code.gene_converter import GeneConverter
            converter = GeneConverter()
            gene_symbols = []
            unrecognized_entrez = []
            
            for gene_id in gene_input:
                symbol = converter.get_symbol(gene_id)
                if symbol:
                    gene_symbols.append(symbol)
                else:
                    unrecognized_entrez.append(gene_id)
            
            logger.info(f"Converted {len(gene_symbols)} Entrez IDs to gene symbols")
            if unrecognized_entrez:
                logger.warning(f"Warning: {len(unrecognized_entrez)} Entrez IDs not found in database")
        else:
            # Use gene symbols directly
            from code.gene_converter import GeneConverter
            converter = GeneConverter()
            gene_symbols = []
            unrecognized_symbols = []
            
            for gene_id in gene_input:
                mapped_symbol = converter.validate_and_map_symbol(gene_id)
                if mapped_symbol:
                    gene_symbols.append(mapped_symbol)
                else:
                    unrecognized_symbols.append(gene_id)
            
            if unrecognized_symbols:
                logger.warning(f"Warning: {len(unrecognized_symbols)} gene symbols not found in database")
        
        # Create GeneSet object
        gene_set = GeneSet(
            gene_symbols,
            background_gene_set.genes,
            gene_set_name,
            hgcn=False,
            format=False,
        )
        
        logger.info(f"Gene set {gene_set_name}: {gene_set.size} genes")
        
        # Check gene list size limit (800 genes maximum)
        if gene_set.size > 800:
            logger.error(f"❌ **Gene list too large!** Your input contains {gene_set.size} genes, but the maximum allowed is 800 genes. Please reduce your gene list size.")
            raise typer.Exit(code=1)
        
        # Create unique results directory for this run with timestamp
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        unique_folder_name = f"{gene_set_name}_{timestamp}"
        gene_set_output_dir = output_dir / unique_folder_name
        gene_set_output_dir.mkdir(exist_ok=True)
        logger.info(f"Created unique output directory: {unique_folder_name}")
        
        if mode == "regular":
            run_regular_enrichment(
                gene_set, gene_set_libraries, background_gene_set,
                p_value_method, p_threshold, min_overlap,
                gene_set_output_dir
            )
        elif mode == "iterative":
            run_iterative_enrichment(
                gene_set, gene_set_libraries, background_gene_set,
                p_value_method, p_threshold, min_overlap, min_term_size, max_term_size,
                max_iterations, gene_set_output_dir
            )


def run_regular_enrichment(
    gene_set: GeneSet,
    gene_set_libraries: List[GeneSetLibrary],
    background_gene_set: BackgroundGeneSet,
    p_value_method: str,
    p_threshold: float,
    min_overlap: int,
    output_dir: Path,
):
    """Run regular enrichment analysis."""
    logger.info("Running regular enrichment analysis")
    
    all_results = []
    
    for library in gene_set_libraries:
        logger.info(f"Processing library: {library.name}")
        
        try:
            enrich = Enrichment(
                gene_set,
                library,
                background_gene_set,
                min_term_size=library.library[0]["size"] if library.library else 10,
                max_term_size=library.library[-1]["size"] if library.library else 600,
                p_value_method_name=p_value_method,
            )
            
            # Filter results by p-value threshold and minimum overlap
            filtered_results = [
                result for result in enrich.results
                if (result.get("p-value", 1.0) <= p_threshold and
                    result.get("overlap_size", "").split("/")[0].isdigit() and 
                    int(result.get("overlap_size", "").split("/")[0]) >= min_overlap)
            ]
            
            # Add library name to each result for combined results
            for result in filtered_results:
                result["library"] = library.name
            
            # Update the enrichment object with filtered results to maintain compatibility
            enrich.results = filtered_results
            
            # Use the same formatting as Streamlit app
            if filtered_results:
                # Get the properly formatted DataFrame using the same method as Streamlit
                results_df = enrich.to_dataframe()
                output_file = output_dir / f"{library.name}_regular_results.tsv"
                results_df.to_csv(output_file, sep='\t', index=False)
                logger.info(f"Saved {len(filtered_results)} results to {output_file}")
                
                # Store for combined results (keep original format for compatibility)
                all_results.extend(filtered_results)
            
        except Exception as e:
            logger.error(f"Error processing library {library.name}: {e}")
            continue
    
    # Save combined results
    if all_results:
        # Create combined results in the same format as Streamlit
        combined_rows = []
        for result in all_results:
            import math
            # Get the library name from the result
            library_name = result.get("library", "")
            # If library name is not in the result, try to get it from the enrichment object
            if not library_name and hasattr(result, 'gene_set_library'):
                library_name = result.gene_set_library.name
            
            combined_rows.append({
                "Library": library_name,
                "Rank": result.get("rank", ""),
                "Term": result.get("term", "").replace("_", " ") if result.get("term") else "",
                "Description": result.get("description", ""),
                "Overlap size": result.get("overlap_size", ""),
                "p-value": result.get("p-value", ""),
                "-log(p-value)": -math.log10(result.get("p-value", 1)) if result.get("p-value", 1) > 0 else 0,
                "FDR": result.get("fdr", ""),
                "Genes": ", ".join(result.get("overlap", [])) if result.get("overlap") else "",
            })
        
        # Reorder columns to match Streamlit format
        column_order = ["Library", "Rank", "Term", "Description", "Overlap size", "p-value", "-log(p-value)", "FDR", "Genes"]
        combined_df = pd.DataFrame(combined_rows)
        if not combined_df.empty:
            combined_df = combined_df[column_order]
        
        combined_file = output_dir / "combined_regular_results.tsv"
        combined_df.to_csv(combined_file, sep='\t', index=False)
        logger.info(f"Saved combined results to {combined_file}")
        
        # Save JSON snapshot
        json_file = output_dir / "regular_enrichment_snapshot.json"
        with open(json_file, 'w') as f:
            json.dump({
                "gene_set_name": gene_set.name,
                "gene_set_size": gene_set.size,
                "background_name": background_gene_set.name,
                "background_size": background_gene_set.size,
                "libraries": [lib.name for lib in gene_set_libraries],
                "parameters": {
                    "p_value_method": p_value_method,
                    "p_threshold": p_threshold,
                    "min_overlap": min_overlap,
                },
                "total_results": len(all_results),
                "timestamp": datetime.now().isoformat()
            }, f, indent=2)
        logger.info(f"Saved metadata to {json_file}")


def _combine_dot_files(all_iter_results: dict) -> str:
    """
    Combine DOT files from multiple libraries into a single DOT file.
    
    Args:
        all_iter_results: Dictionary with library names as keys and dict with 'results' and 'dot_content' as values
        
    Returns:
        Combined DOT content as string
    """
    # Start with the DOT header
    combined_dot = "graph iterative_enrichment {\n"
    combined_dot += "  graph [layout=neato];\n"
    combined_dot += "  node [shape=ellipse];\n"
    
    # Collect all nodes and edges
    all_nodes = set()
    all_edges = set()
    
    for lib_name, data in all_iter_results.items():
        dot_content = data['dot_content']
        
        # Parse the DOT content to extract nodes and edges
        lines = dot_content.split('\n')
        for line in lines:
            line = line.strip()
            if line.startswith('"') and '[' in line and ']' in line:
                # This is a node definition
                all_nodes.add(line)
            elif line.startswith('"') and ' -- ' in line:
                # This is an edge definition
                all_edges.add(line)
    
    # Add all nodes
    for node in sorted(all_nodes):
        combined_dot += f"  {node}\n"
    
    # Add all edges
    for edge in sorted(all_edges):
        combined_dot += f"  {edge}\n"
    
    # Close the graph
    combined_dot += "}\n"
    
    return combined_dot


def run_iterative_enrichment(
    gene_set: GeneSet,
    gene_set_libraries: List[GeneSetLibrary],
    background_gene_set: BackgroundGeneSet,
    p_value_method: str,
    p_threshold: float,
    min_overlap: int,
    min_term_size: int,
    max_term_size: int,
    max_iterations: int,
    output_dir: Path,
):
    """Run iterative enrichment analysis."""
    logger.info("Running iterative enrichment analysis")
    
    # Generate a shared run ID for all libraries
    shared_run_id = datetime.now().strftime("%Y%m%d_%H%M%S_%f")[:-3]
    logger.info(f"Generated shared run ID: {shared_run_id}")
    
    all_iter_results = {}
    
    for library in gene_set_libraries:
        logger.info(f"Processing library: {library.name}")
        
        try:
            # Create progress callback for logging
            def progress_callback(message: str):
                logger.info(f"Library {library.name}: {message}")
            
            it = IterativeEnrichment(
                gene_set=gene_set,
                gene_set_library=library,
                background_gene_set=background_gene_set,
                min_term_size=min_term_size,
                max_term_size=max_term_size,
                p_value_method_name=p_value_method,
                p_threshold=p_threshold,
                max_iterations=None if max_iterations == 0 else max_iterations,
                min_overlap=min_overlap,
                progress_callback=progress_callback,
                run_id=shared_run_id,
            )
            
            # Save individual library results
            if it.results:
                # Use the same formatting as Streamlit app
                results_df = it.to_dataframe()
                
                output_file = output_dir / f"{library.name}_iterative_results.tsv"
                results_df.to_csv(output_file, sep='\t', index=False)
                logger.info(f"Saved {len(it.results)} iterations to {output_file}")
                
                # Store DOT content for combined network
                dot_content = it.to_dot()
                all_iter_results[library.name] = {
                    'results': it.results,
                    'dot_content': dot_content
                }
            
        except Exception as e:
            logger.error(f"Error processing library {library.name}: {e}")
            continue
    
    # Save combined results
    if all_iter_results:
        # Create combined TSV in the same format as Streamlit
        combined_rows = []
        for lib_name, data in all_iter_results.items():
            results = data['results']
            for result in results:
                import math
                combined_rows.append({
                    "Library": lib_name,
                    "Iteration": result.get("Iteration", ""),
                    "Term": result.get("Term", "").replace("_", " ") if result.get("Term") else "",
                    "Description": result.get("Description", ""),
                    "Overlap size": result.get("Overlap size", ""),
                    "p-value": result.get("p-value", ""),
                    "-log(p-value)": -math.log10(result.get("p-value", 1)) if result.get("p-value", 1) > 0 else 0,
                    "Genes": ", ".join(result.get('Genes', [])) if result.get('Genes') else "",
                })
        
        if combined_rows:
            # Reorder columns to match Streamlit format
            column_order = ["Library", "Iteration", "Term", "Description", "Overlap size", "p-value", "-log(p-value)", "Genes"]
            combined_df = pd.DataFrame(combined_rows)
            if not combined_df.empty:
                combined_df = combined_df[column_order]
            
            combined_file = output_dir / "combined_iterative_results.tsv"
            combined_df.to_csv(combined_file, sep='\t', index=False)
            logger.info(f"Saved combined iterative results to {combined_file}")
        
        # Generate combined DOT file
        try:
            combined_dot_content = _combine_dot_files(all_iter_results)
            combined_dot_file = output_dir / "combined_network.dot"
            with open(combined_dot_file, 'w') as f:
                f.write(combined_dot_content)
            logger.info(f"Saved combined network DOT to {combined_dot_file}")
            
            # Generate and save AI analysis prompt for combined network
            try:
                from code.ui.rendering import generate_ai_analysis_prompt
                
                # Enhanced prompt format
                enhanced_prompt = generate_ai_analysis_prompt(combined_dot_content)
                enhanced_file = output_dir / "combined_ai_analysis_prompt_enhanced.txt"
                with open(enhanced_file, 'w') as f:
                    f.write(enhanced_prompt)
                logger.info(f"Saved combined enhanced AI prompt to {enhanced_file}")
                
            except ImportError as e:
                logger.warning(f"Could not import AI analysis functions: {e}")
            except Exception as e:
                logger.warning(f"Error generating AI prompts: {e}")
                
        except Exception as e:
            logger.warning(f"Error generating combined DOT file: {e}")
        
        # Save metadata
        json_file = output_dir / "iterative_enrichment_snapshot.json"
        with open(json_file, 'w') as f:
            json.dump({
                "gene_set_name": gene_set.name,
                "gene_set_size": gene_set.size,
                "background_name": background_gene_set.name,
                "background_size": background_gene_set.size,
                "libraries": list(all_iter_results.keys()),
                "parameters": {
                    "p_value_method": p_value_method,
                    "p_threshold": p_threshold,
                    "min_overlap": min_overlap,
                    "min_term_size": min_term_size,
                    "max_term_size": max_term_size,
                    "max_iterations": max_iterations,
                },
                "total_iterations": sum(len(data['results']) for data in all_iter_results.values()),
                "run_id": shared_run_id,
                "timestamp": datetime.now().isoformat()
            }, f, indent=2)
        logger.info(f"Saved metadata to {json_file}")


@app.command(help="Run gene set enrichment analysis from the command line")
def main(
    gene_sets: List[Path] = typer.Option(
        None,
        "--genelist",
        "-g",
        exists=True,
        file_okay=True,
        dir_okay=False,
        help="Paths to gene set files.",
    ),
    gene_format: str = typer.Option(
        "symbols",
        "--gene-format",
        help="Gene input format: 'symbols' (Gene Symbols) or 'entrez_ids' (Entrez IDs)",
    ),
    background: Optional[Path] = typer.Option(
        None,
        "--background",
        "-b",
        exists=True,
        file_okay=True,
        dir_okay=False,
        help="Path to the background gene set file. Default: all_genes.txt",
        show_default=False,
    ),
    libraries: List[Path] = typer.Option(
        None,
        "--libraries",
        "-l",
        exists=True,
        file_okay=True,
        dir_okay=False,
        help="Paths to gene set library files. Default: active libraries from alias.json.",
        show_default=False,
    ),
    mode: str = typer.Option(
        "regular",
        "--mode",
        "-m",
        help="Analysis mode: 'regular' or 'iterative'",
    ),
    p_value_method: str = typer.Option(
        "Fisher's Exact Test",
        "--method",
        help="P-value calculation method: 'fisher' (Fisher's Exact Test), 'hga' (Hypergeometric Test), or 'chi' (Chi-squared Test)",
    ),
    p_threshold: float = typer.Option(
        0.01,
        "--p-threshold",
        "-p",
        help="Raw p-value threshold for including terms (not FDR-corrected)",
    ),
    min_overlap: int = typer.Option(
        3,
        "--min-overlap",
        help="Minimum overlap size required for terms",
    ),
    min_term_size: int = typer.Option(
        10,
        "--min-term-size",
        help="Minimum term size",
    ),
    max_term_size: int = typer.Option(
        600,
        "--max-term-size",
        help="Maximum term size",
    ),
    max_iterations: int = typer.Option(
        10,
        "--max-iterations",
        help="Maximum iterations for iterative mode (0 = no limit)",
    ),
    output_dir: Path = typer.Option(
        Path("cli_results"),
        "--output-dir",
        "-o",
        help="Output directory for results",
    ),
):
    """
    Run gene set enrichment analysis from the command line.
    
    This CLI provides the same core functionality as the Streamlit app,
    supporting both regular and iterative enrichment analysis modes.
    """
    
    # Convert method shortcuts to full names
    method_mapping = {
        "fisher": "Fisher's Exact Test",
        "hga": "Hypergeometric Test", 
        "chi": "Chi-squared Test"
    }
    
    if p_value_method.lower() in method_mapping:
        p_value_method = method_mapping[p_value_method.lower()]
    
    # Validate gene format
    if gene_format not in ["symbols", "entrez_ids"]:
        typer.echo("Error: Gene format must be 'symbols' or 'entrez_ids'", err=True)
        raise typer.Exit(code=1)
    
    # Validate mode
    if mode not in ["regular", "iterative"]:
        typer.echo("Error: Mode must be 'regular' or 'iterative'", err=True)
        raise typer.Exit(code=1)
    
    # Default values handling
    if not gene_sets:
        gene_sets = list((ROOT / "data/gene_lists").glob("*.txt"))
        if not gene_sets:
            typer.echo("Error: No gene set files found in data/gene_lists/", err=True)
            raise typer.Exit(code=1)

    if not background:
        # Try to use all_genes.txt as default
        all_genes_file = ROOT / "data/backgrounds/all_genes.txt"
        if all_genes_file.exists():
            background = all_genes_file
            typer.echo(f"🎯 Using default background: all_genes.txt")
        else:
            # Fallback to first available background file
            background_files = list((ROOT / "data/backgrounds").glob("*.txt"))
            background = background_files[0] if background_files else None
            if not background:
                typer.echo("Error: No background files found in data/backgrounds/", err=True)
                raise typer.Exit(code=1)
            typer.echo(f"⚠️  all_genes.txt not found, using: {background.name}")

    if not libraries:
        # Load active libraries from alias file
        alias_file = ROOT / "data/libraries/alias.json"
        if not alias_file.exists():
            typer.echo("Error: alias.json file not found in data/libraries/", err=True)
            raise typer.Exit(code=1)
        
        try:
            import json
            with open(alias_file, 'r') as f:
                aliases = json.load(f)
            
            # Get active libraries
            active_libraries = []
            for alias in aliases:
                if alias.get("active", False):
                    lib_file = ROOT / "data/libraries" / alias["file"]
                    if lib_file.exists():
                        active_libraries.append(lib_file)
                    else:
                        typer.echo(f"Warning: Active library file not found: {alias['file']}", err=True)
            
            if not active_libraries:
                typer.echo("Error: No active library files found in alias.json", err=True)
                raise typer.Exit(code=1)
            
            libraries = active_libraries
            typer.echo(f"📚 Using active libraries from alias.json: {len(libraries)} files")
            for lib in libraries:
                typer.echo(f"   - {lib.name}")
            typer.echo("")
            
        except Exception as e:
            typer.echo(f"Error reading alias.json: {e}", err=True)
            raise typer.Exit(code=1)

    # Ensure that the resulting variables are not empty
    if not gene_sets or not libraries:
        typer.echo("Error: Gene sets and libraries cannot be empty.", err=True)
        raise typer.Exit(code=1)

    # Display analysis parameters
    typer.echo(f"Analysis Mode: {mode}")
    typer.echo(f"Gene Sets: {len(gene_sets)} files")
    typer.echo(f"Background: {background}")
    typer.echo(f"Libraries: {len(libraries)} files")
    typer.echo(f"P-value Method: {p_value_method}")
    typer.echo(f"P-value Threshold: {p_threshold}")
    typer.echo(f"Min Overlap: {min_overlap}")
    typer.echo(f"Term Size Range: {min_term_size}-{max_term_size}")
    if mode == "iterative":
        typer.echo(f"Max Iterations: {max_iterations}")
    typer.echo(f"Output Directory: {output_dir}")
    typer.echo("")

    # Run the enrichment analysis
    try:
        run_enrichment(
            gene_sets=gene_sets,
            background=background,
            libraries=libraries,
            gene_format=gene_format,
            p_value_method=p_value_method,
            mode=mode,
            p_threshold=p_threshold,
            min_overlap=min_overlap,
            min_term_size=min_term_size,
            max_term_size=max_term_size,
            max_iterations=max_iterations,
            output_dir=output_dir,
        )
        typer.echo(f"✅ Analysis completed successfully! Results saved to {output_dir}")
    except Exception as e:
        typer.echo(f"❌ Error during analysis: {e}", err=True)
        raise typer.Exit(code=1)


if __name__ == "__main__":
    app()
