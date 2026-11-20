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
3. Select gene set collections and analysis mode (Regular or Iterative).
4. Click **Submit** and explore interactive results.

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
| **Correction** | FDR-adjusted p-values | Empirical q-value per term and network significance, both from random gene lists matched to your run |
| **Best for** | Quick, standard enrichment | Deep, iterative discovery of enriched pathways |

### Gene Set Collections (MSigDB v2025.1)

The app ships with the 10 gene set collections from [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/) used in the iGEA paper:

| MSigDB collection | Gene set collections |
|---|---|
| **Hallmark** | H: Hallmark Gene Sets |
| **Canonical Pathways (C2)** | BioCarta · KEGG Legacy · KEGG MEDICUS · PID · Reactome · WikiPathways |
| **Ontology (C5)** | GO Biological Process · GO Cellular Component · GO Molecular Function |

Other collections in `.gmt` format can be added (see *Adding Custom Gene Set Collections*).

### Background Gene Sets
- **All genes** (default): Uses the full set of annotated human genes.
- **Custom background**: Upload a `.txt` file with all measured genes for more accurate statistics.

### Results & Export
- **Interactive bar charts** for each collection.
- **Export options**: Download results as TSV, JSON, or combined archives (`.tar.gz`).
- **Network visualization**: the network of terms that passed the q-value filter, coloured by community or collection, with hub and bridge genes; the same network is exported for AI analysis.
- **Unique filenames**: All outputs are saved with timestamped, unique identifiers.

### Statistical Methods
- **Fisher's exact test** (one-sided: over-representation; default; same p-values as the hypergeometric test)
- **Hypergeometric Test**
- **Chi-Squared Test**

P-value filtering uses both raw and FDR-adjusted thresholds. Input genes that are not in the background are
listed and removed, and gene sets are restricted to the background, so the input, the gene sets and the random
lists share one universe.

### Network statistics, communities and AI interpretation (iterative mode)
- **Two-level empirical statistics.** Random gene lists of the same size, drawn from your background and analysed
  with the same collections and parameters (1,000 by default in the hosted app, up to 10,000): a q-value for
  each term (terms with q < 0.10 form the network; use 1 to keep all) and an empirical p-value for the largest
  connected cluster of the network.
- **Communities and hub genes.** Leiden communities of terms linked by shared genes; gene degree (number of
  collections placing a gene in the network) and bridge genes between communities.
- **AI-ready export.** A network file for drawing, a network file for AI analysis (communities and statistics as
  attributes) and a complete prompt for a language model: statistics, communities, network and instructions
  to interpret each community, combine them, estimate the experiment and cite network genes.

---

## Project Structure

```
iterative-gene-enrichment-app/
├── streamlit_app.py              # Entry point for Community Cloud (runs code/streamlit_app.py)
├── code/
│   ├── streamlit_app.py          # Main Streamlit application
│   ├── enrichment.py             # ORA enrichment engine
│   ├── iter_enrichment.py        # Iterative enrichment engine
│   ├── fast_igea.py              # Fast engine for random gene lists (identical results)
│   ├── per_run_null.py           # Term q-values and network statistics from random lists
│   ├── network_modules.py        # Communities, hub and bridge genes
│   ├── llm_export.py             # Network file and prompt for AI analysis
│   ├── network_connectivity_benchmark.py  # Largest connected cluster
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
│   └── recent_release/           # NCBI gene info and gene history (.gz)
├── requirements.txt              # Python dependencies
├── pyproject.toml                # Project metadata
└── README.md
```

---

## Deployment on Streamlit Community Cloud

This repository is configured for **zero-setup deployment** on Streamlit Community Cloud.

### Required files

| File | Purpose |
|---|---|
| `requirements.txt` | Python dependencies (installed by Community Cloud; no system packages are needed) |

### Deploy your own instance

1. Fork this repository to your GitHub account.
2. Go to [share.streamlit.io](https://share.streamlit.io) and sign in.
3. Click **New app** → select your fork → set main file to `streamlit_app.py`.
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


---

## Adding Custom Gene Set Collections

1. Place your `.gmt` file in `data/libraries/`.
2. Add an entry to `data/libraries/alias.json`:
   ```json
   {
       "name": "My Custom Collection",
       "file": "my_collection.gmt",
       "active": true,
       "color": "#4CAF50"
   }
   ```
3. Restart the app — the new collection will appear in the selection panel.

Use gene symbols that match the bundled NCBI gene information (the full repository has `code/libraries.py add` to harmonise symbols).

---

## License

This project is licensed under the **GNU General Public License v3.0** — see [LICENSE](LICENSE) for details.
