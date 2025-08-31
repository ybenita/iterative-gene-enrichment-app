# Iterative Gene Enrichment Analysis

## Prerequisites
- Have a list of genes for enrichment analysis in official human gene symbols or entrez IDs.
- Unless you use all genes as a background, have a background file with all measured genes ready in the human gene symbols or entrez IDs.
- You may add new gene set libraries in GMT format. Make sure they are properly listed in the data folder and added to the alias.json file.  
- The NCBI gene information database (`Homo_sapiens.gene_info`) is automatically downloaded for Entrez ID support.

## Run the Application

### Web Interface (Streamlit)
The application provides a user-friendly web interface that runs in your browser:
- **To run locally**: `streamlit run code/streamlit_app.py`
- **To run in Docker**: `docker build -t cc_enrichment:latest .` and `docker run -p 8080:8080 cc_enrichment:latest`

The Streamlit interface offers:
- Interactive web-based UI accessible through your browser
- Real-time parameter adjustment
- Visual results display with charts and graphs
- Easy file upload and download capabilities

### Command Line Interface (CLI)
For automated analysis, batch processing, or integration into workflows, use the command-line interface:

```bash
# Basic usage
python code/cli.py --genelist your_genes.txt --libraries "H: Hallmark Gene Sets"

# Run with all active libraries
python code/cli.py --genelist your_genes.txt

# Iterative enrichment mode
python code/cli.py --genelist your_genes.txt --mode iterative

# See all options
python code/cli.py --help
```

**CLI Features:**
- Full feature parity with the web interface
- Support for both regular and iterative enrichment modes
- Batch processing capabilities
- Automated output generation (TSV, JSON, DOT files)
- Gene format support (symbols or Entrez IDs)
- Configurable statistical parameters

For detailed CLI documentation, see [CLI_README.md](CLI_README.md).

### Choosing Between Web Interface and CLI

| Feature | Web Interface (Streamlit) | Command Line Interface (CLI) |
|---------|---------------------------|------------------------------|
| **Ease of Use** | ✅ User-friendly, no command line knowledge required | ⚠️ Requires basic command line skills |
| **Visualization** | ✅ Interactive charts and graphs | ❌ Text-based output only |
| **Real-time Adjustments** | ✅ Instant parameter changes | ❌ Requires re-running commands |
| **Batch Processing** | ❌ Manual operation only | ✅ Automated processing |
| **Integration** | ❌ Limited to web browser | ✅ Easy integration into workflows |
| **Resource Usage** | ⚠️ Higher memory usage | ✅ Lower resource footprint |
| **Remote Access** | ✅ Accessible from any browser | ⚠️ Requires SSH/remote access |

**Recommendation:**
- **Use Streamlit** for interactive exploration, visualization, and one-off analyses
- **Use CLI** for batch processing, automation, and integration into computational pipelines

## Overview of Application Components
### Key Features

#### Gene Set Submission
Users have the flexibility to input their gene set in multiple ways:
- **Direct Input**: Paste a newline-separated list of gene identifiers and assign a name to the gene set.
- **File Selection**: Choose an existing gene list from the local `data` folder.
- **Multiple Formats**: Support for gene symbols (e.g., TP53, BRCA1) and Entrez IDs (e.g., 7157, 672).
- **Format Selection**: Users must specify their input format before entering gene lists.

#### Gene Set Validation
The app performs crucial validation steps to ensure the accuracy of the analysis:
- **Format Validation**: Validates gene identifiers against the selected input format (symbols or Entrez IDs).
- **Conversion**: Converts Entrez IDs to official gene symbols when using Entrez ID format. In addition, converts old symbols to updated symbols when possible. 
- **Duplicate Check**: Checks for any duplicate gene names within the input gene set.
- **Background Check**: Validates the input gene set against a predefined background gene list, which helps in identifying any genes that might not be present in the background list.

#### Results Display and Interaction
- Interactive Bar Graph: An interactive bar graph provides a visual representation of the results for each gene set library.

#### Data Export Options
- Individual Library Export: Users can save the results for each gene set library in TSV and JSON formats.
- Consolidated Export: There is an option to save results for all libraries in a single TSV file.

#### Results Management
All results, along with the complementary metadata, are stored in the `results` folder and are assigned unique filenames.

### Advanced Features

#### Customization of Results Display
Users can customize the number of results to display in both the chart and the bar graph in the Streamlit app.

#### P-value Calculation Methods
The app offers various statistical methods to calculate p-values:
- Fisher's Exact Test.
- Hypergeometric Test.
- Chi-Squared Test.

**Important Note on P-values**: The p-value threshold used for filtering results are both **raw p-value** and FDR adjusted p-values. In Iterative mode, adjusted p-values are not relevant and thus not calculated. 

#### Custom Reference Data Upload
Background Gene List Upload: Users can upload their background gene lists as .txt files containing newline-separated gene names or Entrez IDs.
Gene Set Libraries Upload: The app allows for the upload of gene set libraries in the .gmt file format, enabling the use of custom gene sets in the analysis.
