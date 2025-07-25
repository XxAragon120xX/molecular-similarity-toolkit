# Molecular Analysis Tools

A comprehensive collection of Python tools for molecular similarity analysis and Morgan fingerprint calculations using RDKit.

## Table of Contents

- [Description](#description)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
  - [Molecular Similarity Search](#molecular-similarity-search)
  - [Morgan Fingerprint Analysis](#morgan-fingerprint-analysis)
- [Project Structure](#project-structure)
- [Examples](#examples)
- [Dataset Configuration](#dataset-configuration)
- [Result Interpretation](#result-interpretation)
- [Contributing](#contributing)
- [License](#license)

## Description

This repository contains two main tools for computational molecular analysis:

1. **molecular_similarity_search.py**: A comprehensive tool for molecular similarity searching using the Tanimoto coefficient
2. **Huellas_morgan.py**: A utility for calculating and analyzing Morgan fingerprints of molecules

These tools are designed for computational chemists, drug discovery researchers, and chemistry students working with molecular similarity analysis.

## Features

### Molecular Similarity Search
- Molecular similarity search in databases using Tanimoto coefficient
- Interactive visualization of results with molecular structures
- Detailed similarity statistics
- CSV export functionality
- Automatic generation of molecular grid images
- Customizable parameter configuration
- Batch processing for multiple queries

### Morgan Fingerprint Analysis
- Morgan fingerprint calculation
- Active bit analysis
- Similarity comparison between bit lists
- Identification of responsible molecular fragments

## Installation

### Prerequisites

Ensure you have Python 3.7+ installed on your system.

### Dependencies

Install the required dependencies:

```bash
pip install pandas matplotlib numpy pillow rdkit-pypi
```

Or using conda (recommended for RDKit):

```bash
conda install -c conda-forge rdkit pandas matplotlib numpy pillow
```

### Clone Repository

```bash
git clone https://github.com/your-username/molecular-analysis-tools.git
cd molecular-analysis-tools
```

## Usage

### Molecular Similarity Search

#### Basic Usage (Interactive Mode)

```bash
python molecular_similarity_search.py
```

The script will prompt you for:
1. Dataset path confirmation (CSV format with columns: 'Smiles', 'Name', 'ChEMBL ID')
2. Query molecule SMILES string

#### Programmatic Usage

```python
from molecular_similarity_search import MolecularSimilaritySearch

# Initialize the tool
search = MolecularSimilaritySearch()

# Load dataset
search.load_dataset('path/to/your/dataset.csv')

# Set query molecule
search.set_query_molecule('CCO')  # Ethanol as example

# Perform search
similarities = search.calculate_similarities()
results = search.filter_and_rank_results(similarities)

# Visualize results
img = search.visualize_top_molecules()
search.display_results(img)
```

#### Custom Configuration

```python
from molecular_similarity_search import Config, MolecularSimilaritySearch

# Create custom configuration
config = Config()
config.SIMILARITY_THRESHOLD = 0.5
config.NUM_MOLECULES_TO_VISUALIZE = 10
config.FINGERPRINT_RADIUS = 3

# Use custom configuration
search = MolecularSimilaritySearch(config)
```

### Morgan Fingerprint Analysis

#### Individual Molecule Analysis

```bash
python Huellas_morgan.py
```

Enter a SMILES string when prompted (example: `CCO` for ethanol).

#### Similarity Comparison

The script also allows comparing similarity between two bit lists:

```python
# Programmatic usage example
from Huellas_morgan import calcular_similitud

lista1 = [1, 2, 3, 4, 5]
lista2 = [3, 4, 5, 6, 7]
similitud = calcular_similitud(lista1, lista2)
print(f"Similarity: {similitud:.2f}%")
```

## Project Structure

```
molecular-analysis-tools/
│
├── molecular_similarity_search.py    # Main search tool
├── Huellas_morgan.py                 # Morgan fingerprint analysis
├── README.md                         # This file
├── requirements.txt                  # Project dependencies
│
├── RT/                              # Output folder (auto-generated)
│   ├── filtered_compounds.csv       # Filtered results
│   └── displayed_compounds.csv      # Displayed compounds
│
└── structures_png/                  # Generated images (auto-created)
    └── molecules.png                # Similar molecules grid
```

## Examples

### Example 1: Paracetamol Similarity Search

```python
# Paracetamol SMILES: CC(=O)Nc1ccc(O)cc1
search = MolecularSimilaritySearch()
search.load_dataset('compounds_database.csv')
search.set_query_molecule('CC(=O)Nc1ccc(O)cc1')

similarities = search.calculate_similarities()
results = search.filter_and_rank_results(similarities, threshold=0.3)

print(f"Found {len(results)} similar molecules")
```

### Example 2: Batch Analysis

```python
from molecular_similarity_search import batch_similarity_search

# List of SMILES to analyze
query_molecules = [
    'CCO',                    # Ethanol
    'CC(=O)Nc1ccc(O)cc1',    # Paracetamol
    'CC(C)Cc1ccc(C(C)C(=O)O)cc1'  # Ibuprofen
]

results = batch_similarity_search(
    query_molecules, 
    dataset_path='database.csv',
    threshold=0.2
)
```

### Example 3: Advanced Configuration

```python
# Configuration for high-precision analysis
config = Config()
config.SIMILARITY_THRESHOLD = 0.7
config.FINGERPRINT_RADIUS = 3
config.FINGERPRINT_BITS = 4096
config.NUM_MOLECULES_TO_VISUALIZE = 30

search = MolecularSimilaritySearch(config)
```

## Dataset Configuration

Your CSV file should contain the following columns:
- `Smiles`: SMILES strings of molecules
- `Name`: Compound names
- `ChEMBL ID`: ChEMBL identifiers (or unique IDs)

Example:
```csv
Smiles,Name,ChEMBL ID
CCO,Ethanol,CHEMBL545
CC(=O)Nc1ccc(O)cc1,Paracetamol,CHEMBL112
```

## Result Interpretation

### Tanimoto Coefficient
- **1.0**: Identical molecules
- **0.7-0.9**: Very similar
- **0.5-0.7**: Moderately similar
- **0.3-0.5**: Slightly similar
- **< 0.3**: Poorly similar

### Output Files
- `filtered_compounds.csv`: All compounds exceeding similarity threshold
- `displayed_compounds.csv`: Detailed information of displayed molecules
- `molecules.png`: Grid visualization of most similar molecules

## Configuration Parameters

### MolecularSimilaritySearch Configuration

| Parameter | Default | Description |
|-----------|---------|-------------|
| `SIMILARITY_THRESHOLD` | 0.1 | Minimum similarity threshold |
| `NUM_MOLECULES_TO_VISUALIZE` | 20 | Number of molecules to display |
| `FINGERPRINT_RADIUS` | 2 | Morgan fingerprint radius |
| `FINGERPRINT_BITS` | 2048 | Number of bits in fingerprint |
| `SUBIMG_SIZE` | (250, 250) | Size of individual molecule images |
| `MOLS_PER_ROW` | 3 | Molecules per row in grid |

## Performance Considerations

- For large datasets (>100k compounds), consider increasing `FINGERPRINT_BITS` to 4096
- Batch processing is more efficient for multiple queries on the same dataset
- Memory usage scales with dataset size and fingerprint bit count

## Error Handling

The tools include comprehensive error handling for:
- Invalid SMILES strings
- Missing dataset columns
- File I/O errors
- RDKit molecule parsing failures

## Dependencies

### Core Dependencies
- **RDKit**: Cheminformatics toolkit
- **pandas**: Data manipulation and analysis
- **matplotlib**: Plotting library
- **numpy**: Numerical computing
- **PIL (Pillow)**: Image processing

### Optional Dependencies
- **jupyter**: For notebook-based analysis
- **seaborn**: Enhanced statistical plotting

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.


## Contact

For questions or suggestions, please create an issue on GitHub.
