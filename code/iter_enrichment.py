import logging
import re
from datetime import datetime
from typing import Any, Dict, List, Optional, Set, Callable

import pandas as pd

from code.background_gene_set import BackgroundGeneSet
from code.enrichment import Enrichment, format_term_name
from code.gene_set import GeneSet
from code.gene_set_library import GeneSetLibrary

logger = logging.getLogger(__name__)


class IterativeEnrichment:
    """
    Perform iterative enrichment analysis by repeatedly removing genes from the top hit.
    """

    def __init__(
        self,
        gene_set: GeneSet,
        gene_set_library: GeneSetLibrary,
        background_gene_set: BackgroundGeneSet,
        min_term_size: int = 10,
        max_term_size: int = 600,
        p_value_method_name: str = "Fisher's Exact Test",
        name: str = None,
        p_threshold: float = 0.01,
        max_iterations: Optional[int] = None,
        min_overlap: int = 1,
        progress_callback: Optional[Callable[[str], None]] = None,
        run_id: Optional[str] = None,
    ) -> None:
        """
        Initialize iterative enrichment.

        :param gene_set: Input gene set
        :param gene_set_library: Gene set library
        :param background_gene_set: Background gene set
        :param min_term_size: Minimum term size
        :param max_term_size: Maximum term size
        :param p_value_method_name: P-value calculation method
        :param name: Name for the enrichment
        :param p_threshold: P-value cutoff for including terms
        :param max_iterations: Maximum number of iterations (None for no limit)
        :param min_overlap: Minimum overlap size required for terms
        :param progress_callback: Optional callback function to report progress
        """
        self.gene_set = gene_set
        self.gene_set_library = gene_set_library
        self.min_term_size = min_term_size
        self.max_term_size = max_term_size
        self.background_gene_set = background_gene_set
        self.p_value_method_name: str = p_value_method_name
        self.p_threshold: float = p_threshold
        self.max_iterations: Optional[int] = max_iterations
        self.min_overlap: int = min_overlap
        self.progress_callback = progress_callback
        from datetime import datetime
        self.name = (
            name
            if name
            else f"{gene_set.name}_{gene_set_library.name}_{background_gene_set.name}_{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
        )
        self._iteration_enrichments: List[Enrichment] = []  # Store full enrichment results for each iteration
        
        # Use provided run ID or generate unique run ID for this session
        from datetime import datetime
        self._run_id = run_id if run_id is not None else datetime.now().strftime("%Y%m%d_%H%M%S_%f")[:-3]  # Include milliseconds
        
        self._results: List[Dict[str, Any]] = self._compute_enrichment()

    def _get_run_results_dir(self) -> "Path":
        """
        Get the unique results directory for this run.
        
        Returns:
            Path to the run-specific results directory
        """
        from pathlib import Path
        # Use absolute path to project root, same as regular enrichment
        ROOT = Path(__file__).resolve().parent.parent
        results_dir = ROOT / "results" / f"run_{self._run_id}"
        return results_dir

    @property
    def results(self) -> List[Dict[str, Any]]:
        """
        Get the iterative enrichment results.

        :returns: List of iteration records
        :rtype: List[Dict[str, Any]]
        """
        return self._results

    @results.setter
    def results(self, value: List[Dict[str, Any]]) -> None:
        """
        Set the iterative enrichment results.

        :param value: List of iteration records
        :type value: List[Dict[str, Any]]
        """
        self._results = value

    def _compute_enrichment(self) -> List[Dict[str, Any]]:
        """
        Perform iterative enrichment, peeling off top terms until no further terms pass p-value threshold.

        :returns: List of iteration records
        :rtype: List[Dict[str, Any]]
        """
        remaining: Set[str] = set(self.gene_set.genes)
        iteration: int = 1
        records: List[Dict[str, Any]] = []
        
        # Report initial status
        if self.progress_callback:
            self.progress_callback(f"Starting iterative enrichment with {len(remaining)} genes")

        while True:
            if not remaining:
                logger.info("No genes left; stopping iterative enrichment.")
                if self.progress_callback:
                    self.progress_callback("No genes left; stopping iterative enrichment.")
                break
            if self.max_iterations is not None and iteration > self.max_iterations:
                logger.warning("Reached max_iterations; stopping iterative enrichment.")
                if self.progress_callback:
                    self.progress_callback(f"Reached max iterations ({self.max_iterations}); stopping.")
                break

            # Report iteration progress
            if self.progress_callback:
                self.progress_callback(f"Iteration {iteration}: Processing {len(remaining)} remaining genes")

            # Create a new GeneSet with the remaining genes, but avoid re-validation
            # since these genes were already validated when the original gene set was created
            current_set = GeneSet(
                list(remaining), 
                set(self.background_gene_set.genes),
                hgcn=False,  # Don't re-validate against background
                format=False  # Don't re-format genes
            )
            try:
                enr = Enrichment(
                    gene_set=current_set,
                    gene_set_library=self.gene_set_library,
                    min_term_size=self.min_term_size,
                    max_term_size=self.max_term_size,
                    background_gene_set=self.background_gene_set,
                    p_value_method_name=self.p_value_method_name,
                )
            except Exception as e:
                logger.error(f"Enrichment failed at iteration {iteration}: {e}")
                if self.progress_callback:
                    self.progress_callback(f"Enrichment failed at iteration {iteration}: {e}")
                break

            results = enr.results
            if not results:
                logger.info("No enrichment results; terminating.")
                if self.progress_callback:
                    self.progress_callback("No enrichment results; terminating.")
                break

            # Filter results by minimum overlap (same as regular mode)
            filtered_results = [
                result for result in results
                if (result.get("overlap_size", "").split("/")[0].isdigit() and 
                    int(result.get("overlap_size", "").split("/")[0]) >= self.min_overlap)
            ]
            
            if not filtered_results:
                logger.info(f"No results meet minimum overlap requirement ({self.min_overlap}); terminating.")
                if self.progress_callback:
                    self.progress_callback(f"No results meet minimum overlap requirement ({self.min_overlap}); terminating.")
                break

            top = filtered_results[0]
            pval = top.get("p-value")
            if pval is None or pval >= self.p_threshold:
                logger.info("Top term p-value >= threshold; terminating.")
                if self.progress_callback:
                    self.progress_callback(f"Top term p-value ({pval}) >= threshold ({self.p_threshold}); terminating.")
                break

            # Get overlap size from overlap_size field (format: "3/50")
            overlap_size_str = top.get("overlap_size", "0/0")
            try:
                overlap_count = int(overlap_size_str.split("/")[0])
            except (ValueError, IndexError):
                logger.warning(f"Could not parse overlap_size: {overlap_size_str}")
                overlap_count = 0
            
            # Get overlap genes for the record
            overlap_data = top.get("overlap", [])
            if isinstance(overlap_data, str):
                # If overlap is a string, try to parse it
                genes_in_term = set(overlap_data.split(',') if overlap_data else [])
            elif isinstance(overlap_data, list):
                genes_in_term = set(overlap_data)
            else:
                logger.warning(f"Unexpected overlap data type: {type(overlap_data)}")
                genes_in_term = set()
            
            # Debug logging
            logger.info(f"Iteration {iteration}: Top term '{top.get('term', '')}' has p-value {pval} and overlap size {overlap_count}")
            logger.info(f"Overlap genes: {genes_in_term}")
            logger.info(f"Minimum overlap requirement: {self.min_overlap}")
            
            # Report iteration result
            if self.progress_callback:
                self.progress_callback(f"Iteration {iteration}: Found term '{top.get('term', '')}' (p={pval:.4f}, overlap={overlap_count})")
            
            # Note: We already filtered by minimum overlap above, so this check is redundant but kept for safety
            if overlap_count < self.min_overlap:
                logger.info(f"Top term overlap size ({overlap_count}) < minimum overlap requirement ({self.min_overlap}); terminating.")
                if self.progress_callback:
                    self.progress_callback(f"Top term overlap size ({overlap_count}) < minimum overlap requirement ({self.min_overlap}); terminating.")
                break
            
            record: Dict[str, Any] = {
                "Iteration": iteration,
                "Term": format_term_name(top.get("term", "")),
                "Description": top.get("description", ""),
                "Library": self.gene_set_library.name,
                "p-value": pval,
                "Overlap size": top.get("overlap_size", "0/0"),
                "Genes": sorted(genes_in_term),
            }
            records.append(record)
            
            # Save this iteration's enrichment results
            self._save_iteration_results(enr, iteration)
            
            # Store the enrichment object for later export
            self._iteration_enrichments.append(enr)
            
            remaining -= genes_in_term
            iteration += 1

        # Report final status
        if self.progress_callback:
            self.progress_callback(f"Completed iterative enrichment: {len(records)} iterations, {len(remaining)} genes remaining")

        return records

    def _save_iteration_results(self, enrichment: Enrichment, iteration: int) -> None:
        """
        Save individual iteration enrichment results to a file.
        
        Args:
            enrichment: The Enrichment object for this iteration
            iteration: The iteration number
        """
        import json
        from pathlib import Path
        
        # Ensure results directory exists
        results_dir = self._get_run_results_dir()
        results_dir.mkdir(exist_ok=True)
        
        # Create filename with library name and iteration number
        # Replace problematic file name characters with dots
        library_name = self.gene_set_library.name
        for char in ['/', '\\', ':', '*', '?', '"', '<', '>', '|', ' ']:
            library_name = library_name.replace(char, ".")
        
        # Collapse multiple consecutive dots into a single dot
        import re
        library_name = re.sub(r'\.+', '.', library_name)
        filename = f"{library_name}_iteration_{iteration:03d}.json"
        filepath = self._get_run_results_dir() / filename
        
        # Create a snapshot with iteration information
        snapshot = enrichment.to_snapshot()
        snapshot["iteration"] = iteration
        # Get the term from the top result of this enrichment
        top_result = enrichment.results[0] if enrichment.results else {"term": "Unknown"}
        snapshot["iteration_term"] = top_result.get("term", "Unknown")
        
        # Add iteration number to each result
        for result in snapshot[self.gene_set_library.name]:
            result["iteration"] = iteration
        
        # Save to file
        with open(filepath, "w") as f:
            json.dump(snapshot, f, indent=2)
        
        # Also save as TSV with iteration number as first column
        tsv_filename = f"{library_name}_iteration_{iteration:03d}.tsv"
        tsv_filepath = self._get_run_results_dir() / tsv_filename
        
        # Create TSV data with library name as first column
        import math
        tsv_data = []
        
        # Get the genes from the top result of this iteration (same as used in main results)
        top_result = enrichment.results[0] if enrichment.results else {}
        genes_in_term = []
        if top_result:
            overlap_data = top_result.get("overlap", [])
            if isinstance(overlap_data, str):
                # If overlap is a string, try to parse it
                genes_in_term = sorted(overlap_data.split(',') if overlap_data else [])
            elif isinstance(overlap_data, list):
                genes_in_term = sorted(overlap_data)
            else:
                genes_in_term = []
        
        for result in enrichment.results:
            p_value = result.get("p-value", 0)
            if isinstance(p_value, str):
                try:
                    p_value = float(p_value)
                except ValueError:
                    p_value = 0
            
            tsv_data.append({
                "Library": self.gene_set_library.name,
                "Iteration": iteration,
                "Rank": result.get("rank", ""),
                "Term": format_term_name(result.get("term", "")),
                "Description": result.get("description", ""),
                "Overlap size": result.get("overlap_size", ""),
                "Genes": ", ".join(genes_in_term),
                "p-value": result.get("p-value", ""),
                "-log(p-value)": -math.log10(p_value) if p_value > 0 else 0,
            })
        
        # Save TSV file
        import pandas as pd
        df = pd.DataFrame(tsv_data)
        # Reorder columns to put Library first, then Iteration, then Term, then Description, then add -log(p-value) after p-value
        column_order = ["Library", "Iteration", "Rank", "Term", "Description", "Overlap size", "p-value", "-log(p-value)", "Genes"]
        df = df[column_order]
        df.to_csv(tsv_filepath, sep="\t", index=False)
        
        logger.info(f"Saved iteration {iteration} results to {filepath} and {tsv_filepath}")

    def create_iteration_results_archive(self) -> str:
        """
        Create a tar.gz archive of all iteration results for this library.
        
        Returns:
            Path to the created archive file
        """
        import tarfile
        from pathlib import Path
        from datetime import datetime
        
        # Ensure results directory exists
        results_dir = self._get_run_results_dir()
        results_dir.mkdir(exist_ok=True)
        
        # Create archive filename
        library_name = self.gene_set_library.name.replace(" ", "_").replace("/", "_")
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        archive_name = f"{library_name}_iterations_{timestamp}.tar.gz"
        archive_path = results_dir / archive_name
        
        # Create tar.gz archive
        with tarfile.open(archive_path, "w:gz") as tar:
            # Find all iteration files for this library (both JSON and TSV)
            json_pattern = f"{library_name}_iteration_*.json"
            tsv_pattern = f"{library_name}_iteration_*.tsv"
            
            for file_path in results_dir.glob(json_pattern):
                tar.add(file_path, arcname=file_path.name)
            for file_path in results_dir.glob(tsv_pattern):
                tar.add(file_path, arcname=file_path.name)
        
        logger.info(f"Created iteration results archive: {archive_path}")
        return str(archive_path)

    def export_iteration_results_tsv(self) -> str:
        """
        Export all iteration results as TSV with iteration number as first column.
        
        Returns:
            TSV string with iteration results
        """
        import pandas as pd
        
        if not self.results or not self._iteration_enrichments:
            return ""
        
        # Create a list to hold all iteration data
        all_iteration_data = []
        
        # For each iteration, get the full enrichment results and add iteration number
        import math
        for i, (iteration_record, enrichment) in enumerate(zip(self.results, self._iteration_enrichments)):
            iteration_num = iteration_record["Iteration"]
            
            # Get the genes from the iteration record (same as used in main results)
            genes_in_term = iteration_record.get("Genes", [])
            
            # Get all results from this iteration's enrichment
            for result in enrichment.results:
                p_value = result.get("p-value", 0)
                if isinstance(p_value, str):
                    try:
                        p_value = float(p_value)
                    except ValueError:
                        p_value = 0
                
                all_iteration_data.append({
                    "Library": self.gene_set_library.name,
                    "Iteration": iteration_num,
                    "Rank": result.get("rank", ""),
                    "Term": format_term_name(result.get("term", "")),
                    "Description": result.get("description", ""),
                    "Overlap size": result.get("overlap_size", ""),
                    "Genes": ", ".join(genes_in_term),
                    "p-value": result.get("p-value", ""),
                    "-log(p-value)": -math.log10(p_value) if p_value > 0 else 0,
                })
        
        # Create DataFrame and export to TSV
        df = pd.DataFrame(all_iteration_data)
        # Reorder columns to put Library first, then Iteration, then Rank, then Term, then Description, then add -log(p-value) after p-value
        column_order = ["Library", "Iteration", "Rank", "Term", "Description", "Overlap size", "p-value", "-log(p-value)", "Genes"]
        df = df[column_order]
        return df.to_csv(sep="\t", index=False)

    def save_to_results_folder(self) -> None:
        """
        Save the main iterative enrichment summary files to the results folder.
        This creates files similar to regular mode enrichment for consistency.
        """
        import json
        from pathlib import Path
        
        # Ensure results directory exists
        results_dir = self._get_run_results_dir()
        results_dir.mkdir(exist_ok=True)
        
        # Create filename with library name and timestamp
        # Replace problematic file name characters with dots
        library_name = self.gene_set_library.name
        for char in ['/', '\\', ':', '*', '?', '"', '<', '>', '|', ' ']:
            library_name = library_name.replace(char, ".")
        
        # Collapse multiple consecutive dots into a single dot
        import re
        library_name = re.sub(r'\.+', '.', library_name)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # 1. Save main iterative enrichment summary as JSON
        summary_filename = f"{library_name}_iterative_enrichment_{timestamp}.json"
        summary_filepath = results_dir / summary_filename
        
        summary_data = {
            "input_gene_set": list(self.gene_set.genes),
            "background": self.background_gene_set.name,
            "background_size": self.background_gene_set.size,
            "library": self.gene_set_library.name,
            "library_size": self.gene_set_library.size,
            "p_value_method": self.p_value_method_name,
            "p_threshold": self.p_threshold,
            "max_iterations": self.max_iterations,
            "min_overlap": self.min_overlap,
            "min_term_size": self.min_term_size,
            "max_term_size": self.max_term_size,
            "total_iterations": len(self.results),
            "iterations": self.results,
            "final_remaining_genes": list(set(self.gene_set.genes) - set().union(*[set(record.get("Genes", [])) for record in self.results]))
        }
        
        with open(summary_filepath, "w") as f:
            json.dump(summary_data, f, indent=2)
        
        # 2. Save main iterative enrichment summary as TSV
        tsv_filename = f"{library_name}_iterative_enrichment_{timestamp}.tsv"
        tsv_filepath = results_dir / tsv_filename
        
        # Create TSV with iteration summary data
        import math
        tsv_data = []
        for record in self.results:
            p_value = record.get("p-value", 0)
            if isinstance(p_value, str):
                try:
                    p_value = float(p_value)
                except ValueError:
                    p_value = 0
            
            tsv_data.append({
                "Library": record.get("Library", ""),
                "Iteration": record.get("Iteration", ""),
                "Term": record.get("Term", ""),
                "Description": record.get("Description", ""),
                "Overlap size": record.get("Overlap size", ""),
                "p-value": record.get("p-value", ""),
                "-log(p-value)": -math.log10(p_value) if p_value > 0 else 0,
                "Genes": ", ".join(record.get("Genes", [])),
            })
        
        import pandas as pd
        df = pd.DataFrame(tsv_data)
        # Reorder columns to put Library first, then Iteration, then Term, then Description, then add -log(p-value) after p-value
        column_order = ["Library", "Iteration", "Term", "Description", "Overlap size", "p-value", "-log(p-value)", "Genes"]
        df = df[column_order]
        df.to_csv(tsv_filepath, sep="\t", index=False)
        
        logger.info(f"Saved iterative enrichment summary to {summary_filepath} and {tsv_filepath}")

    def to_dataframe(self) -> pd.DataFrame:
        """
        Convert iteration records to a pandas DataFrame.

        :returns: DataFrame of iteration records
        :rtype: pandas.DataFrame
        """
        import math
        
        df = pd.DataFrame(self.results)
        
        # Convert genes from list format to comma-separated format
        if not df.empty and 'Genes' in df.columns:
            df['Genes'] = df['Genes'].apply(
                lambda x: ', '.join(x) if isinstance(x, list) else str(x)
            )
        
        # Ensure correct column order if columns exist
        if not df.empty:
            expected_columns = ["Library", "Iteration", "Term", "Description", "Overlap size", "p-value", "-log(p-value)", "Genes"]
            existing_columns = [col for col in expected_columns if col in df.columns]
            if existing_columns:
                df = df[existing_columns]
        return df

    def to_tsv(self) -> str:
        """
        Serialize iteration records as a TSV string.

        :returns: TSV-formatted string
        :rtype: str
        """
        return self.to_dataframe().to_csv(sep="\t", index=False)

    def to_json(self) -> str:
        """
        Serialize iteration records as a JSON-formatted string.

        :returns: JSON-formatted string
        :rtype: str
        """
        import json

        return json.dumps(self.results, indent=2)

    def to_dot(self) -> str:
        """
        Generate a valid Graphviz DOT for the iterative enrichment network,
        with sanitized, quoted IDs, semicolons, and no duplicates.
        """

        def _sanitize_id(raw: str) -> str:
            """
            Convert raw label into a valid DOT node ID: replace non-alphanumeric with underscores.
            Collapse multiple underscores and strip leading/trailing underscores.
            """
            # replace non-word characters with underscore
            s = re.sub(r"\W+", "_", raw)
            # collapse underscores
            s = re.sub(r"_+", "_", s)
            return s.strip("_")

        def _format_term_name(term_name: str) -> str:
            """
            Format term name by replacing underscores with spaces and adding colon after library shortcut.
            """
            # Import the format_term_name function from enrichment module
            from code.enrichment import format_term_name
            return format_term_name(term_name)

        nodes: Set[str] = set()
        edges: Set[str] = set()

        # Build nodes and edges
        for rec in self.results:
            term_label = rec.get("Term", "")
            # Format term name for display (convert underscores to spaces)
            formatted_term_label = _format_term_name(term_label)
            # sanitize and quote term ID
            raw_id = f"term_{rec['Iteration']}_{term_label}"
            term_id = _sanitize_id(raw_id)
            term_node = (
                f'"{term_id}" '
                f'[label="{formatted_term_label}", '
                f'style=filled, fontcolor="white", type="term"];'
            )
            nodes.add(term_node)

            for gene in rec.get("Genes", []):
                gene_id = _sanitize_id(f"gene_{gene}")
                gene_node = f'"{gene_id}" [label="{gene}", type="gene"];'
                nodes.add(gene_node)
                # create edge with quoted IDs and semicolon
                edge = f'"{gene_id}" -- "{term_id}";'
                edges.add(edge)

        # Assemble DOT
        lines: List[str] = []
        lines.append("graph iterative_enrichment {")
        lines.append("  graph [layout=neato];")
        lines.append("  node [shape=ellipse];")

        # add nodes
        for node in sorted(nodes):
            lines.append(f"  {node}")
        # add edges
        for edge in sorted(edges):
            lines.append(f"  {edge}")

        lines.append("}")
        return "\n".join(lines)
