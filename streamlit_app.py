"""
Root-level entry point for Streamlit Cloud.
Adds the code/ directory to sys.path so all internal imports resolve,
then runs the real application.
"""
import sys
from pathlib import Path

# Add code/ directory to Python path so bare imports work
sys.path.insert(0, str(Path(__file__).resolve().parent / "code"))

# Import and run the actual app
from streamlit_app import *  # noqa: E402, F401, F403
