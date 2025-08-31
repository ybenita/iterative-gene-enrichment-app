# Gene Enrichment Analysis - Streamlit Cloud Deployment

This repository contains the Streamlit Cloud deployment of the Iterative Gene Enrichment Analysis tool.

## Live Application

The application is deployed on Streamlit Cloud and can be accessed at the URL provided by Streamlit Cloud after deployment.

## Features

- **Gene Set Enrichment Analysis**: Perform enrichment analysis using multiple statistical methods
- **Iterative Enrichment**: Advanced iterative analysis with network visualization
- **Multiple Libraries**: Support for MSigDB, Reactome, and other gene set libraries
- **Interactive Visualizations**: Plotly charts and network diagrams
- **Export Options**: Download results in CSV, TSV, and JSON formats
- **Gene ID Conversion**: Support for various gene identifier formats

## Deployment Structure

```
/
├── streamlit_app.py          # Main Streamlit application
├── requirements.txt          # Python dependencies
├── README.md                 # Project documentation
├── LICENSE                   # GPL-3.0 license
├── .gitignore               # Git ignore rules
├── code/                    # Core application code
│   ├── ui/                 # UI components
│   ├── enrichment.py       # Enrichment analysis
│   ├── gene_converter.py   # Gene ID conversion
│   └── ...                # Other modules
└── data/                   # Data files
    ├── libraries/         # Gene set libraries
    ├── gene_lists/        # Example gene lists
    └── backgrounds/       # Background gene sets
```

## Streamlit Cloud Configuration

- **Main file path**: `streamlit_app.py`
- **Python version**: 3.12+
- **Dependencies**: `requirements.txt`

## Usage

1. **Upload Gene List**: Paste or upload a list of gene identifiers
2. **Select Libraries**: Choose from available gene set libraries
3. **Configure Analysis**: Set parameters for enrichment analysis
4. **Run Analysis**: Execute the enrichment analysis
5. **View Results**: Explore interactive visualizations and download results

## Included Libraries

- **H: Hallmark Gene Sets** - Curated gene sets representing well-defined biological states
- **C1: Positional Gene Sets** - Gene sets by chromosome location  
- **C2: Curated Gene Sets** - Gene sets from various sources including Reactome
- **C5: GO Gene Sets** - Gene Ontology biological processes

## Example Gene Lists

- `example_gene_list.txt` - General example gene list
- `HIV.InputGeneList.txt` - HIV-related genes
- `hypoxia-genes.symbols.txt` - Hypoxia-related genes

## Technical Details

- **Framework**: Streamlit for web interface
- **Analysis**: Custom Python modules for statistical analysis
- **Visualization**: Plotly and NetworkX for interactive plots
- **Data Format**: GMT format for gene set libraries

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Support

For technical support or questions about the analysis methods, please refer to the main project documentation or contact the development team.
