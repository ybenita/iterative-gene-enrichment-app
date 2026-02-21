"""
Root-level entry point for Streamlit Cloud.
Adds the code/ directory to sys.path so all internal imports resolve,
then runs the real application from code/streamlit_app.py.
"""
import importlib.util
import sys
from pathlib import Path

_code_dir = str(Path(__file__).resolve().parent / "code")

# Add code/ directory to Python path so bare imports work
# (background_gene_set, enrichment, gene_set, ui.*, etc.)
sys.path.insert(0, _code_dir)

# Load code/streamlit_app.py by file path to avoid self-import
_spec = importlib.util.spec_from_file_location(
    "_inner_app", str(Path(_code_dir) / "streamlit_app.py")
)
_app = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_app)

_app.main()
