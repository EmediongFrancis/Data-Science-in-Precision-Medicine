# Getting Started with Transcriptomic Data Analysis

Welcome to the Introduction to Transcriptomics module! This guide will help you quickly get up and running with RNA-seq data analysis using Python.

## Quick Setup

### 1. Install Dependencies

```bash
# Install basic requirements
pip install -r requirements.txt

# Optional: Install bioinformatics tools
pip install scanpy anndata gseapy bioservices

# Optional: Install advanced visualization tools
pip install plotly bokeh umap-learn
```

### 2. Test Your Installation

```bash
# Run the demonstration script
python demo.py

# Start Jupyter notebook
jupyter notebook notebooks/transcriptomic_data_analysis.ipynb
```

## Your First Analysis

### Load and Explore Data

```python
from scripts.transcriptomic_utils import TranscriptomeDataProcessor
from scripts.config import generate_sample_count_matrix, get_sample_metadata_df

# Initialize processor
processor = TranscriptomeDataProcessor()

# Generate sample data (or load your own)
count_data = generate_sample_count_matrix(n_genes=5000, n_samples=12)
metadata = get_sample_metadata_df()

print(f"Data shape: {count_data.shape}")
print(f"Conditions: {metadata['condition'].unique()}")
```

### Quality Control and Normalization

```python
# Calculate library sizes
library_sizes = processor.calculate_library_sizes(count_data)
print(f"Library size range: {library_sizes.min():,.0f} - {library_sizes.max():,.0f}")

# Filter low-expression genes
filtered_data = processor.filter_low_expression_genes(
    count_data, min_count=10, min_samples=3
)

# Normalize for analysis
normalized_data = processor.normalize_counts(filtered_data, method='log_cpm')
```

### Differential Expression Analysis

```python
from scripts.transcriptomic_utils import DifferentialExpressionAnalyzer

# Initialize analyzer
de_analyzer = DifferentialExpressionAnalyzer()

# Define sample groups
treatment_samples = ['OHT_1', 'OHT_2', 'OHT_3', 'OHT_4']
control_samples = ['Control_1', 'Control_2', 'Control_3', 'Control_4']

# Perform analysis
results = de_analyzer.perform_deseq2_like_analysis(
    filtered_data, treatment_samples, control_samples
)

# Identify significant genes
significant = de_analyzer.identify_top_genes(
    results, p_threshold=0.05, fc_threshold=1.5
)

print(f"Found {len(significant)} significantly differentially expressed genes")
```

### Pathway Enrichment

```python
from scripts.transcriptomic_utils import PathwayEnrichmentAnalyzer

# Initialize pathway analyzer
pathway_analyzer = PathwayEnrichmentAnalyzer()
pathway_analyzer.load_gene_sets()

# Get list of differentially expressed genes
de_genes = significant['gene'].tolist()

# Perform enrichment analysis
enrichment = pathway_analyzer.perform_overrepresentation_analysis(de_genes)
print(f"Found {len(enrichment)} enriched pathways")
```

### Create Visualizations

```python
from scripts.transcriptomic_utils import TranscriptomeVisualizer

# Initialize visualizer
visualizer = TranscriptomeVisualizer()

# PCA plot
fig1 = visualizer.plot_pca(normalized_data, metadata=metadata, color_by='condition')

# Volcano plot
fig2 = visualizer.plot_volcano(results)

# Expression heatmap
fig3 = visualizer.plot_heatmap(normalized_data, top_n=50)
```

## Common Use Cases

### 1. Compare Treatment vs Control

```python
# Standard two-group comparison
treatment_samples = metadata[metadata['condition'] == 'OHT'].index.tolist()
control_samples = metadata[metadata['condition'] == 'Control'].index.tolist()

de_results = de_analyzer.perform_deseq2_like_analysis(
    filtered_data, treatment_samples, control_samples
)
```

### 2. Time Course Analysis

```python
# For time series data, compare each timepoint to baseline
timepoints = ['0h', '6h', '12h', '24h']
baseline = '0h'

for timepoint in timepoints[1:]:
    treatment = metadata[metadata['timepoint'] == timepoint].index.tolist()
    control = metadata[metadata['timepoint'] == baseline].index.tolist()
    
    results = de_analyzer.perform_deseq2_like_analysis(
        filtered_data, treatment, control
    )
    print(f"{timepoint} vs {baseline}: {len(results)} genes analyzed")
```

### 3. Multi-Factor Analysis

```python
# Account for batch effects or paired samples
# Use the patient factor in your experimental design
for condition in ['OHT', 'DPN']:
    treatment_samples = [s for s, info in metadata.iterrows() 
                        if info['condition'] == condition]
    control_samples = [s for s, info in metadata.iterrows() 
                      if info['condition'] == 'Control' and 
                         info['patient'] in [metadata.loc[t, 'patient'] for t in treatment_samples]]
```

## Working with Your Own Data

### Data Format Requirements

**Count Matrix (TSV format):**
```
gene_id    Sample_1    Sample_2    Sample_3    ...
GENE1      1523        1876        1234        ...
GENE2      234         456         789         ...
...
```

**Sample Metadata (TSV format):**
```
sample_id    condition    patient    batch
Sample_1     Control      Patient1   Batch1
Sample_2     Treatment    Patient1   Batch1
...
```

### Loading Your Data

```python
# Load count matrix
processor = TranscriptomeDataProcessor()
count_data = processor.load_count_matrix('path/to/counts.tsv')

# Load metadata
metadata = processor.load_sample_metadata('path/to/metadata.tsv')

# Validate data
from scripts.config import validate_count_matrix
validation = validate_count_matrix(count_data)
if not validation['valid']:
    print("Data validation issues:", validation['errors'])
```

### Customizing Analysis Parameters

```python
from scripts.config import ANALYSIS_PARAMS

# Modify parameters
ANALYSIS_PARAMS['differential_expression']['p_value_threshold'] = 0.01
ANALYSIS_PARAMS['differential_expression']['fold_change_threshold'] = 2.0

# Use custom gene sets
custom_pathways = {
    'My_Pathway': ['GENE1', 'GENE2', 'GENE3'],
    'Another_Pathway': ['GENE4', 'GENE5', 'GENE6']
}
pathway_analyzer.load_gene_sets(custom_pathways)
```

## Troubleshooting

### Common Issues

**1. Import Errors**
```bash
# Make sure you're in the correct directory
cd 03-introduction-to-transcriptomics

# Check Python path
python -c "import sys; print('\n'.join(sys.path))"

# Install missing dependencies
pip install pandas numpy matplotlib seaborn scipy scikit-learn
```

**2. Memory Issues with Large Datasets**
```python
# Process data in chunks for large datasets
chunk_size = 1000
for i in range(0, len(count_data), chunk_size):
    chunk = count_data.iloc[i:i+chunk_size]
    # Process chunk
```

**3. Visualization Problems**
```python
# Set matplotlib backend if plots don't display
import matplotlib
matplotlib.use('Agg')  # For saving plots without display

# Or use interactive backend
matplotlib.use('TkAgg')  # For interactive plots
```

**4. Statistical Issues**
```python
# Check for sufficient samples
min_samples_per_group = 3
if len(treatment_samples) < min_samples_per_group:
    print("Warning: Insufficient samples for robust statistical analysis")

# Handle missing values
normalized_data = normalized_data.fillna(0)

# Check for extreme outliers
outliers = (normalized_data > normalized_data.quantile(0.99)).sum()
print(f"Potential outliers: {outliers.sum()} values")
```

### Performance Optimization

```python
# For large datasets, consider:
# 1. Filtering more aggressively
filtered_data = processor.filter_low_expression_genes(
    count_data, min_count=20, min_samples=5
)

# 2. Using fewer genes for visualization
top_variable_genes = normalized_data.var(axis=1).nlargest(1000).index
plot_data = normalized_data.loc[top_variable_genes]

# 3. Sampling for pathway analysis
if len(de_genes) > 1000:
    de_genes_sample = np.random.choice(de_genes, 1000, replace=False)
```

## Next Steps

### 1. Explore Advanced Features
- Single-cell RNA-seq analysis with scanpy
- Time-course and trajectory analysis
- Integration with other omics data
- Custom pathway databases

### 2. Apply to Your Research
- Adapt workflows to your experimental design
- Customize visualizations for publications
- Develop analysis pipelines for routine use
- Integrate with laboratory information systems

### 3. Learn More
- Read the comprehensive documentation in README.md
- Work through the Jupyter notebook examples
- Explore the demo script source code
- Check out the advanced analysis options

## Getting Help

### Resources
- **Documentation**: See README.md for comprehensive information
- **Examples**: Run `python demo.py` for working examples
- **Code**: All source code is documented and available in `scripts/`

### Common Questions

**Q: How do I handle batch effects?**
A: Include batch as a factor in your experimental design and use appropriate statistical models.

**Q: What if I don't have replicates?**
A: Biological replicates are essential for statistical analysis. Technical replicates alone are insufficient.

**Q: How do I choose normalization methods?**
A: log_cpm is good for most analyses, TPM for between-sample comparisons, raw counts for DESeq2-style analysis.

**Q: Can I analyze single-cell data?**
A: Yes, but you'll need additional tools like scanpy. This module focuses on bulk RNA-seq.

---

Ready to start your transcriptomic analysis journey? Begin with the Jupyter notebook for an interactive experience! 