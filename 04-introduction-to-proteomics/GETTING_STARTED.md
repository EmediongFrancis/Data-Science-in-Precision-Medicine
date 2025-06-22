# Getting Started with Proteomic Data Analysis

Welcome to the Introduction to Proteomics module! This guide will help you quickly get up and running with mass spectrometry-based protein analysis using Python.

## Quick Setup

### 1. Install Dependencies

```bash
# Install basic requirements
pip install -r requirements.txt

# Optional: Install proteomics-specific tools
pip install pyteomics pymzml spectrum_utils

# Optional: Install advanced visualization tools
pip install plotly bokeh networkx
```

### 2. Test Your Installation

```bash
# Run the demonstration script
python demo.py

# Start Jupyter notebook
jupyter notebook notebooks/proteomic_data_analysis.ipynb
```

## Your First Analysis

### Load and Explore Protein Data

```python
from scripts.proteomic_utils import ProteomeDataProcessor
from scripts.config import generate_sample_abundance_matrix, get_sample_metadata_df

# Initialize processor
processor = ProteomeDataProcessor()

# Generate sample data (or load your own)
abundance_data = generate_sample_abundance_matrix(n_proteins=2000, n_samples=12)
metadata = get_sample_metadata_df()

print(f"Data shape: {abundance_data.shape}")
print(f"Conditions: {metadata['condition'].unique()}")
print(f"Missing values: {abundance_data.isnull().sum().sum() / abundance_data.size * 100:.1f}%")
```

### Handle Missing Values and Normalize

```python
# Handle missing values (common in proteomics)
processed_data = processor.handle_missing_values(
    abundance_data, method='knn', missing_threshold=0.5
)

# Calculate protein statistics
protein_stats = processor.calculate_protein_statistics(processed_data)
print(f"Median CV: {protein_stats['cv'].median():.2f}")
print(f"Average detection frequency: {protein_stats['detection_freq'].mean():.2f}")

# Normalize data
normalized_data = processor.normalize_abundance_data(
    processed_data, method='median_centering'
)
```

### Differential Protein Expression Analysis

```python
from scripts.proteomic_utils import DifferentialProteinAnalyzer

# Initialize analyzer
de_analyzer = DifferentialProteinAnalyzer()

# Define sample groups
infection_samples = ['Patient_Z_Infection_Early', 'Patient_Z_Infection_Peak']
healthy_samples = ['Patient_Z_Baseline', 'Patient_Z_Recovered']

# Perform analysis
results = de_analyzer.perform_ttest_analysis(
    normalized_data, infection_samples, healthy_samples
)

# Identify significant proteins
significant = de_analyzer.identify_significant_proteins(
    results, p_threshold=0.05, fc_threshold=1.5
)

print(f"Found {len(significant)} significantly differentially expressed proteins")
```

### Pathway Enrichment Analysis

```python
from scripts.proteomic_utils import ProteinPathwayAnalyzer

# Initialize pathway analyzer
pathway_analyzer = ProteinPathwayAnalyzer()
pathway_analyzer.load_pathway_databases()

# Get list of differentially expressed proteins
de_proteins = significant['protein'].tolist()

# Perform enrichment analysis
enrichment = pathway_analyzer.perform_pathway_enrichment(de_proteins)
print(f"Found {len(enrichment)} enriched pathways")
```

### Create Visualizations

```python
from scripts.proteomic_utils import ProteomicsVisualizer

# Initialize visualizer
visualizer = ProteomicsVisualizer()

# PCA plot
fig1 = visualizer.plot_pca(normalized_data, metadata=metadata, color_by='condition')

# Volcano plot
fig2 = visualizer.plot_volcano(results)

# Protein abundance heatmap
fig3 = visualizer.plot_abundance_heatmap(normalized_data, top_n=50)
```

## Common Use Cases

### 1. Acute Phase Response Analysis

```python
# Focus on acute phase proteins
acute_phase_proteins = ['CRP', 'SAA1', 'HP', 'ALB', 'APOA1']

# Compare infection vs healthy
infection_samples = metadata[metadata['condition'] == 'Infection'].index.tolist()
healthy_samples = metadata[metadata['condition'] == 'Healthy'].index.tolist()

de_results = de_analyzer.perform_ttest_analysis(
    normalized_data, infection_samples, healthy_samples
)

# Analyze specific proteins
for protein in acute_phase_proteins:
    if protein in de_results['protein'].values:
        result = de_results[de_results['protein'] == protein].iloc[0]
        print(f"{protein}: {result['log2_fold_change']:.2f} fold change, p={result['p_adj']:.3f}")
```

### 2. Time Course Analysis

```python
# Analyze Patient Z's timeline
timepoints = ['Baseline', 'Early', 'Peak', 'Recovery_1', 'Recovery_2', 'Recovered']
baseline = 'Baseline'

for timepoint in timepoints[1:]:
    treatment_samples = [f'Patient_Z_Infection_{timepoint}' if 'Recovery' not in timepoint 
                        else f'Patient_Z_{timepoint}']
    control_samples = ['Patient_Z_Baseline']
    
    if all(s in normalized_data.columns for s in treatment_samples + control_samples):
        results = de_analyzer.perform_ttest_analysis(
            normalized_data, treatment_samples, control_samples
        )
        print(f"{timepoint} vs {baseline}: {len(results)} proteins analyzed")
```

### 3. Missing Value Pattern Analysis

```python
# Analyze missing value patterns
missing_per_protein = abundance_data.isnull().sum(axis=1)
missing_per_sample = abundance_data.isnull().sum(axis=0)

print(f"Proteins with >50% missing: {(missing_per_protein > 6).sum()}")
print(f"Samples with >20% missing: {(missing_per_sample > abundance_data.shape[0]*0.2).sum()}")

# Visualize missing patterns
fig = visualizer.plot_missing_values(abundance_data)
```

## Working with Your Own Data

### Data Format Requirements

**Protein Abundance Matrix (TSV format):**
```
Protein.IDs    Sample_1    Sample_2    Sample_3    ...
P12345         23.4        24.1        22.8        ...
P67890         18.2        17.9        18.5        ...
...
```

**Sample Metadata (TSV format):**
```
sample_id      condition    timepoint    patient
Sample_1       Control      Single       Patient1
Sample_2       Infection    Early        Patient2
...
```

### Loading Your Data

```python
# Load abundance matrix
processor = ProteomeDataProcessor()
abundance_data = processor.load_abundance_matrix('path/to/abundance.tsv')

# Load metadata
metadata = processor.load_sample_metadata('path/to/metadata.tsv')

# Validate data
from scripts.config import validate_abundance_matrix
validation = validate_abundance_matrix(abundance_data)
if not validation['valid']:
    print("Data validation issues:", validation['errors'])
```

### Customizing Analysis Parameters

```python
from scripts.config import ANALYSIS_PARAMS

# Modify parameters
ANALYSIS_PARAMS['differential_expression']['p_value_threshold'] = 0.01
ANALYSIS_PARAMS['differential_expression']['fold_change_threshold'] = 2.0

# Use custom protein pathways
custom_pathways = {
    'My_Pathway': ['PROTEIN1', 'PROTEIN2', 'PROTEIN3'],
    'Another_Pathway': ['PROTEIN4', 'PROTEIN5', 'PROTEIN6']
}
pathway_analyzer.load_pathway_databases(custom_pathways)
```

## Understanding Proteomics Data

### Mass Spectrometry Basics

```python
# Typical abundance values (log2 transformed)
print("Typical protein abundance ranges:")
print("- High abundance (e.g., Albumin): 25-30")
print("- Medium abundance (e.g., Enzymes): 18-25") 
print("- Low abundance (e.g., Cytokines): 10-18")
print("- Missing values: Often due to detection limits")

# Check your data range
print(f"Your data range: {abundance_data.min().min():.1f} to {abundance_data.max().max():.1f}")
```

### Missing Value Strategies

```python
# Different strategies for missing values
strategies = ['knn', 'min', 'median', 'drop']

for strategy in strategies:
    processed = processor.handle_missing_values(abundance_data, method=strategy)
    print(f"{strategy}: {processed.shape[0]} proteins retained")
```

### Normalization Methods

```python
# Compare normalization methods
methods = ['median_centering', 'quantile', 'zscore']

for method in methods:
    normalized = processor.normalize_abundance_data(processed_data, method=method)
    print(f"{method}: range {normalized.min().min():.2f} to {normalized.max().max():.2f}")
```

## Troubleshooting

### Common Issues

**1. Import Errors**
```bash
# Make sure you're in the correct directory
cd 04-introduction-to-proteomics

# Check Python path
python -c "import sys; print('\n'.join(sys.path))"

# Install missing dependencies
pip install pandas numpy matplotlib seaborn scipy scikit-learn
```

**2. High Missing Value Percentage**
```python
# Check missing value distribution
missing_pct = abundance_data.isnull().sum(axis=1) / abundance_data.shape[1]
print(f"Proteins with <50% missing: {(missing_pct < 0.5).sum()}")

# Filter more aggressively
filtered_data = processor.handle_missing_values(
    abundance_data, missing_threshold=0.3  # More stringent
)
```

**3. No Significant Results**
```python
# Check sample sizes
print(f"Group 1 samples: {len(group1_samples)}")
print(f"Group 2 samples: {len(group2_samples)}")

# Relax thresholds if needed
significant = de_analyzer.identify_significant_proteins(
    results, p_threshold=0.1, fc_threshold=1.2  # Less stringent
)
```

**4. Visualization Issues**
```python
# Set matplotlib backend if plots don't display
import matplotlib
matplotlib.use('Agg')  # For saving plots without display

# Handle missing values before plotting
plot_data = normalized_data.fillna(0)
```

### Performance Optimization

```python
# For large datasets:
# 1. Filter proteins more aggressively
filtered_data = processor.filter_proteins(
    abundance_data, min_detection_freq=0.7, max_cv=0.8
)

# 2. Use subset for visualization
top_variable = abundance_data.var(axis=1).nlargest(500).index
plot_data = abundance_data.loc[top_variable]

# 3. Sample proteins for pathway analysis
if len(de_proteins) > 500:
    de_proteins_sample = np.random.choice(de_proteins, 500, replace=False)
```

## Understanding Results

### Interpreting Differential Expression

```python
# Understanding fold changes in proteomics
print("Fold change interpretation:")
print("- 2-fold (log2FC = 1): Moderate change")
print("- 4-fold (log2FC = 2): Strong change") 
print("- 10-fold (log2FC = 3.3): Dramatic change (e.g., CRP in infection)")

# Effect sizes
print("Effect sizes (Cohen's d):")
print("- 0.2: Small effect")
print("- 0.5: Medium effect")
print("- 0.8: Large effect")
```

### Pathway Enrichment Interpretation

```python
# Understanding enrichment results
print("Enrichment interpretation:")
print("- Enrichment ratio > 2: Strong enrichment")
print("- P-value < 0.05: Statistically significant")
print("- Overlap size: Number of proteins in both lists")

# Clinical relevance
print("Clinically relevant pathways:")
print("- Acute Phase Response: Infection/inflammation")
print("- Complement System: Immune activation")
print("- Blood Coagulation: Hemostasis changes")
```

## Next Steps

### 1. Explore Advanced Features
- Multi-timepoint analysis for longitudinal studies
- Integration with clinical parameters
- Protein-protein interaction networks
- Custom pathway databases

### 2. Apply to Your Research
- Adapt workflows to your experimental design
- Develop analysis pipelines for routine use
- Create custom visualizations for publications
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

**Q: How do I handle batch effects in proteomics data?**
A: Consider including batch as a factor in your design matrix and use appropriate statistical models.

**Q: What's the difference between missing values and low abundance?**
A: Missing values are below detection limits; low abundance values are quantified but near the limit.

**Q: How do I choose the right normalization method?**
A: Median centering is good for most analyses, quantile normalization for between-sample comparisons.

**Q: Can I analyze single-cell proteomics data?**
A: This module focuses on bulk proteomics. Single-cell proteomics requires specialized methods.

---

Ready to start your proteomics analysis journey? Begin with the Jupyter notebook for an interactive experience! 