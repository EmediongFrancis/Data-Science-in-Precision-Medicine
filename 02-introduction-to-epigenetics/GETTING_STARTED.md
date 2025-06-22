# Getting Started with Epigenetics Analysis

**Quick start guide for the Introduction to Epigenetics module**

## Installation

1. **Install Python dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

2. **Verify installation:**
   ```bash
   python3 -c "import pandas, numpy, scipy, sklearn; print('All dependencies installed!')"
   ```

## Quick Start

### Option 1: Run the Demo Script
```bash
python demo.py
```

This will generate sample data and create several visualization files demonstrating epigenetic analysis capabilities.

### Option 2: Use the Jupyter Notebook
```bash
jupyter notebook notebooks/epigenetic_data_analysis.ipynb
```

The notebook provides an interactive introduction to epigenetic data analysis with detailed explanations.

### Option 3: Use the Python API
```python
from scripts.epigenetic_utils import EpigenomeDataProcessor
from scripts.config import generate_sample_methylation_data

# Generate sample data
data = generate_sample_methylation_data(n_regions=1000, n_samples=20)

# Analyze methylation patterns
processor = EpigenomeDataProcessor()
stats = processor.calculate_methylation_statistics(data)
print(f"Analyzed {stats['total_regions']} regions across {stats['total_samples']} samples")
```

## What You'll Learn

- **DNA Methylation Analysis**: Understand how cytosine methylation regulates gene expression
- **Histone Modifications**: Explore the histone code and chromatin states
- **Chromatin State Classification**: Learn to identify active, repressed, and poised genomic regions
- **Epigenetic Regulation**: Discover how epigenetics controls development and disease
- **Data Integration**: Combine multiple epigenomic datasets for comprehensive analysis

## Key Concepts

### DNA Methylation
- **CpG dinucleotides**: Primary targets for methylation
- **Gene silencing**: Methylation typically represses gene expression
- **Tissue specificity**: Different cell types have distinct methylation patterns
- **Disease relevance**: Aberrant methylation contributes to cancer and other diseases

### Histone Modifications
- **H3K4me3**: Active promoter mark
- **H3K27ac**: Active enhancer mark
- **H3K27me3**: Repressive mark (Polycomb)
- **H3K9me3**: Heterochromatin mark
- **Combinatorial code**: Multiple marks define chromatin states

### Chromatin States
- **Active Promoter**: H3K4me3 + H3K27ac, low methylation
- **Strong Enhancer**: H3K4me1 + H3K27ac, low methylation
- **Polycomb Repressed**: H3K27me3, variable methylation
- **Heterochromatin**: H3K9me3, high methylation
- **Quiescent**: Low signal for all marks

## Sample Analyses

### 1. Methylation Distribution Analysis
```python
from scripts.epigenetic_utils import EpigenomeVisualizer
from scripts.config import generate_sample_methylation_data

# Generate and visualize methylation data
data = generate_sample_methylation_data(n_regions=500)
visualizer = EpigenomeVisualizer()
fig = visualizer.plot_methylation_distribution(data)
```

### 2. Chromatin State Classification
```python
from scripts.epigenetic_utils import HistoneModificationAnalyzer
from scripts.config import generate_sample_histone_data

# Generate histone data and classify states
histone_data = generate_sample_histone_data(n_regions=100)
analyzer = HistoneModificationAnalyzer()

# Summarize by gene and classify
gene_summary = histone_data.groupby('gene')[['H3K4me3', 'H3K4me1', 'H3K27ac', 'H3K27me3', 'H3K9me3']].mean()
classified = analyzer.classify_chromatin_states(gene_summary)
print(classified['chromatin_state'].value_counts())
```

### 3. Gene-Specific Analysis
```python
from scripts.config import DEMO_GENES, get_gene_info

# Analyze pluripotency genes
for gene_name in DEMO_GENES['pluripotency']:
    gene_info = get_gene_info(gene_name, 'pluripotency')
    print(f"{gene_info['name']}: {gene_info['function']}")
```

## Generated Files

Running the demo creates several visualization files:

- `methylation_distribution.png` - Distribution of methylation levels
- `methylation_pca.png` - PCA analysis of methylation patterns
- `histone_profile.png` - Histone modification profiles around genes
- `chromatin_states.png` - Chromatin state distribution
- `epigenetic_correlations.png` - Correlations between epigenetic marks
- `methylation_vs_h3k4me3.png` - Relationship between methylation and H3K4me3

## Troubleshooting

### Common Issues

**Import Errors:**
```bash
# Install missing dependencies
pip install pandas numpy matplotlib seaborn scipy scikit-learn
```

**Visualization Issues:**
```bash
# For macOS users, you might need:
pip install --upgrade matplotlib
```

**Memory Issues:**
```python
# Reduce dataset size for testing
data = generate_sample_methylation_data(n_regions=100, n_samples=10)
```

### Getting Help

1. **Check the notebook** - Contains detailed explanations and examples
2. **Review the demo script** - Shows complete analysis workflows
3. **Examine the config file** - Contains sample data and parameters
4. **Read the documentation** - Each function has detailed docstrings

## Next Steps

1. **Explore real data** - Download datasets from the NIH Roadmap Epigenomics Program
2. **Analyze your own data** - Adapt the tools for your specific research questions
3. **Extend the analysis** - Add new histone marks or chromatin states
4. **Apply machine learning** - Use the features for predictive modeling

## Learning Path

1. **Start with the notebook** - Interactive introduction to concepts
2. **Run the demo script** - See complete analysis workflows
3. **Experiment with parameters** - Modify analysis settings
4. **Try real datasets** - Apply to published epigenomic data
5. **Develop new methods** - Extend the toolkit for your research

---

*Ready to explore epigenetics? Start with the Jupyter notebook and then run the demo script to see the full analysis pipeline in action!* 