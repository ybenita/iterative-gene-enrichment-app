# Iterative Gene Enrichment Analysis

A web-based tool for gene set enrichment analysis with support for both standard over-representation analysis (ORA) and a novel **iterative enrichment** mode. Deployed on [Streamlit Community Cloud](https://streamlit.io/cloud).

---

## 🚀 Live App

Access the application directly in your browser — no installation required:

**[Launch Enrichment Analysis →](https://gene-enrichment-analysis.streamlit.app)**

---

## Quick Start

1. Prepare a gene list in **official human gene symbols** (e.g., TP53, BRCA1) or **Entrez IDs** (e.g., 7157, 672).
2. Open the app and paste your genes or upload a file.
3. Select gene set libraries and analysis mode (Regular or Iterative).
4. Click **Run Enrichment** and explore interactive results.

---

## Features

### Gene Input & Validation
- **Flexible input**: Paste a newline-separated gene list or upload a `.txt` file.
- **Dual format support**: Gene symbols and Entrez IDs with automatic validation.
- **Smart conversion**: Entrez IDs → official symbols; outdated symbols → current symbols via NCBI gene information.
- **Duplicate detection**: Automatically flags and handles duplicate entries.

### Analysis Modes

| | Regular (ORA) | Iterative Enrichment |
|---|---|---|
| **Method** | Standard over-representation analysis | Multi-iteration refinement of gene-term associations |
| **P-value** | Fisher's exact, hypergeometric, or chi-squared | Fisher's exact, hypergeometric, or chi-squared |
| **Correction** | FDR-adjusted p-values | Raw p-values per iteration (FDR not applicable) |
| **Best for** | Quick, standard enrichment | Deep, iterative discovery of enriched pathways |

### Gene Set Libraries (MSigDB v2025.1)

The app ships with 12 pre-loaded gene set libraries from [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/):

| Collection | Library |
|---|---|
| **Hallmark** | H: Hallmark Gene Sets |
| **Canonical Pathways (C2)** | BioCarta · KEGG Legacy · KEGG MEDICUS · PID · Reactome · WikiPathways |
| **Ontology (C5)** | GO Biological Process · GO Cellular Component · GO Molecular Function · Human Phenotype Ontology |
| **Protein Interaction** | StringDB Protein Interactions |

Custom libraries in `.gmt` format can be uploaded through the UI.

### Background Gene Sets
- **All genes** (default): Uses the full set of annotated human genes.
- **Custom background**: Upload a `.txt` file with all measured genes for more accurate statistics.

### Results & Export
- **Interactive bar charts** for each library with adjustable display settings.
- **Export options**: Download results as TSV, JSON, or combined archives (`.tar.gz`).
- **Network visualization**: Explore gene-term relationships as interactive network graphs.
- **Unique filenames**: All outputs are saved with timestamped, unique identifiers.

### Statistical Methods
- **Fisher's Exact Test**
- **Hypergeometric Test**
- **Chi-Squared Test**

P-value filtering uses both raw and FDR-adjusted thresholds. In Iterative mode, only raw p-values are used.

---

## Project Structure

```
iterative_enrichment/
├── code/
│   ├── streamlit_app.py          # Main Streamlit application
│   ├── enrichment.py             # ORA enrichment engine
│   ├── iter_enrichment.py        # Iterative enrichment engine
│   ├── gene_converter.py         # Gene symbol / Entrez ID converter
│   ├── gene_set.py               # Gene set data model
│   ├── gene_set_library.py       # GMT library parser
│   ├── background_gene_set.py    # Background gene set loader
│   └── ui/                       # UI components and rendering
├── data/
│   ├── libraries/                # GMT gene set library files
│   │   └── alias.json            # Library metadata & activation config
│   ├── backgrounds/              # Background gene set files
│   ├── gene_lists/               # Example gene lists
│   └── recent_release/           # NCBI gene info files
├── .streamlit/
│   └── config.toml               # Streamlit theme & server config
├── requirements.txt              # Python dependencies
├── packages.txt                  # System packages (apt) for Cloud
├── pyproject.toml                # Project metadata
└── README.md
```

---

## Deployment on Streamlit Community Cloud

This repository is configured for **zero-setup deployment** on Streamlit Community Cloud.

### Required files

| File | Purpose |
|---|---|
| `requirements.txt` | Python dependencies (installed via `pip`) |
| `packages.txt` | System-level packages (installed via `apt`) — includes `graphviz` |
| `.streamlit/config.toml` | Theme, upload limits, and server settings |

### Deploy your own instance

1. Fork this repository to your GitHub account.
2. Go to [share.streamlit.io](https://share.streamlit.io) and sign in.
3. Click **New app** → select your fork → set main file to `code/streamlit_app.py`.
4. Click **Deploy**. Streamlit will install dependencies and launch the app.

### Running locally

```bash
# Clone the repository
git clone https://github.com/<your-org>/iterative-enrichment.git
cd iterative-enrichment

# Create and activate a virtual environment (Python 3.12.4+)
python -m venv .venv
source .venv/bin/activate

# Install dependencies
pip install -r requirements.txt

# Launch the app
streamlit run code/streamlit_app.py
```

> **Note**: Graphviz must be installed on your system for network visualization (`brew install graphviz` on macOS, `sudo apt install graphviz graphviz-dev` on Linux).

---

## Adding Custom Gene Set Libraries

1. Place your `.gmt` file in `data/libraries/`.
2. Add an entry to `data/libraries/alias.json`:
   ```json
   {
       "name": "My Custom Library",
       "file": "my_library.gmt",
       "active": true,
       "color": "#4CAF50"
   }
   ```
3. Restart the app — the new library will appear in the selection panel.

Alternatively, upload a `.gmt` file directly through the app's **Upload Library** feature.

---

## License

This project is licensed under the **GNU General Public License v3.0** — see [LICENSE](LICENSE) for details.
