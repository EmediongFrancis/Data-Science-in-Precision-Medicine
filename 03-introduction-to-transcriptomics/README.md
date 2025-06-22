# Introduction to Transcriptomics

**Stanford Data Ocean - Data Science in Precision Medicine**

## Overview

Transcriptomics represents the comprehensive study of RNA transcripts produced by the genome under specific conditions. Unlike the static genome we explored in Module 1 or the dynamic epigenome from Module 2, the transcriptome provides a real-time snapshot of gene expression - showing which genes are active, when, and to what extent.

This module provides hands-on experience with transcriptomic data analysis, focusing on bulk RNA-seq and single-cell RNA-seq analysis using Python-based tools for differential expression analysis, pathway enrichment, and biological interpretation.

## Learning Objectives

By completing this module, you will:

- **Understand transcriptomic technologies** - Learn the principles of RNA-seq, bulk vs single-cell approaches
- **Analyze real RNA-seq datasets** - Work with count matrices, perform normalization and quality control
- **Identify differentially expressed genes** - Statistical methods for comparing expression between conditions
- **Perform pathway analysis** - Connect gene expression changes to biological processes
- **Master visualization techniques** - Create publication-quality plots for transcriptomic data
- **Apply computational methods** - Use Python tools for comprehensive transcriptomic analysis

## Module Structure

```
03-introduction-to-transcriptomics/
├── scripts/
│   ├── transcriptomic_utils.py      # Core analysis functions
│   ├── config.py                    # Configuration and constants
│   └── visualization_utils.py       # Plotting and visualization tools
├── notebooks/
│   └── transcriptomic_data_analysis.ipynb  # Main analysis notebook
├── data/                            # Sample datasets
├── demo.py                          # Standalone demonstration script
├── requirements.txt                 # Python dependencies
├── setup.py                        # Package installation
└── README.md                        # This file
```

## Key Concepts Covered

### Transcriptomic Technologies
- **Bulk RNA-seq** - Population-level gene expression profiling
- **Single-cell RNA-seq (scRNA-seq)** - Individual cell transcriptome analysis
- **Single-nucleus RNA-seq (snRNA-seq)** - Nuclear transcriptome profiling
- **Spatial transcriptomics** - Location-aware gene expression analysis

### Analysis Workflows
- **Quality Control** - Assessing sequencing quality and sample integrity
- **Normalization** - Correcting for technical variation and sequencing depth
- **Differential Expression** - Statistical identification of expression changes
- **Functional Analysis** - Pathway enrichment and gene ontology analysis

### Biological Applications
- **Disease mechanisms** - Understanding pathological gene expression changes
- **Drug discovery** - Identifying therapeutic targets and biomarkers
- **Development** - Tracking gene expression during differentiation
- **Cell type identification** - Discovering and characterizing cell populations

## Datasets

This module uses data representing common transcriptomic study designs:

- **Hormone treatment study** - OHT vs DPN vs Control comparisons
- **Time-course experiments** - Gene expression dynamics over time
- **Single-cell datasets** - Cell type identification and trajectory analysis
- **Clinical samples** - Disease vs healthy tissue comparisons

## Technical Requirements

### Python Dependencies
```
pandas >= 1.3.0
numpy >= 1.21.0
matplotlib >= 3.4.0
seaborn >= 0.11.0
scipy >= 1.7.0
scikit-learn >= 1.0.0
statsmodels >= 0.12.0
```

### Bioinformatics Libraries
```
scanpy >= 1.8.0
anndata >= 0.8.0
gseapy >= 0.10.0
bioservices >= 1.8.0
```

## Quick Start

1. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

2. **Run the demonstration:**
   ```bash
   python demo.py
   ```

3. **Explore the notebook:**
   ```bash
   jupyter notebook notebooks/transcriptomic_data_analysis.ipynb
   ```

## Data Sources

- **Gene Expression Omnibus (GEO)**: https://www.ncbi.nlm.nih.gov/geo/
- **The Cancer Genome Atlas (TCGA)**: https://portal.gdc.cancer.gov/
- **Human Cell Atlas**: https://www.humancellatlas.org/
- **Reference**: Various published transcriptomic studies

## Key Analysis Workflows

### 1. Bulk RNA-seq Analysis
- Load and process count matrices
- Perform quality control and normalization
- Identify differentially expressed genes
- Visualize results with volcano plots and heatmaps

### 2. Pathway Enrichment Analysis
- Gene Ontology (GO) term enrichment
- KEGG pathway analysis
- Gene Set Enrichment Analysis (GSEA)
- Custom gene set testing

### 3. Single-cell Analysis
- Cell filtering and quality control
- Dimensionality reduction (PCA, UMAP)
- Cell clustering and annotation
- Trajectory analysis and pseudotime

### 4. Comparative Analysis
- Multi-condition comparisons
- Time-course analysis
- Integration with other omics data
- Clinical correlation analysis

## Example Applications

### Disease Research
Examine gene expression changes in disease:
- **Cancer transcriptomics**: Tumor vs normal tissue analysis
- **Inflammatory diseases**: Immune response characterization
- **Neurological disorders**: Brain region-specific expression changes
- **Metabolic diseases**: Tissue-specific metabolic pathway alterations

### Drug Discovery
Investigate therapeutic mechanisms:
- **Target identification**: Finding druggable genes and pathways
- **Biomarker discovery**: Predictive and prognostic signatures
- **Mechanism of action**: Understanding drug effects on gene expression
- **Resistance mechanisms**: Identifying pathways involved in treatment failure

## Visualization Capabilities

The module includes comprehensive visualization tools for:
- **Expression heatmaps** across samples and genes
- **Volcano plots** for differential expression results
- **PCA plots** for sample relationship visualization
- **Pathway enrichment plots** showing biological processes
- **Single-cell visualizations** including UMAP and t-SNE plots

## Educational Features

- **Progressive complexity** - Builds from basic concepts to advanced analysis
- **Real-world examples** - Uses actual research datasets and findings
- **Interactive exploration** - Hands-on analysis with immediate feedback
- **Biological context** - Connects computational results to biological meaning
- **Best practices** - Demonstrates proper statistical methods and interpretation

## Further Reading

### Core Papers
- Conesa et al. (2016). A survey of best practices for RNA-seq data analysis. *Genome Biology*
- Love et al. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*
- Wolf et al. (2018). SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*

### Methodological Resources
- RNA-seq analysis best practices from ENCODE
- Single-cell RNA-seq analysis workflows
- Pathway analysis methods and databases

## Support and Troubleshooting

### Common Issues
- **Memory requirements**: Large transcriptomic datasets may require significant RAM
- **Statistical considerations**: Multiple testing correction and effect size interpretation
- **Batch effects**: Technical variation between samples or experiments

### Getting Help
- Check the troubleshooting section in the notebook
- Review the example outputs in the demo script
- Consult the documentation in the utility modules

## Contributions

This module was developed as part of the Stanford Data Ocean certificate program in Data Science for Precision Medicine.

**Module Development Team:**
- **Content Development**: Based on Ramesh V Nair's original R module
- **Python Implementation**: Comprehensive Python translation and enhancement
- **Technical Implementation**: Advanced statistical methods and visualizations
- **Educational Design**: Progressive learning with biological context
- **Data Curation**: Curated datasets representing diverse study designs

---

*Ready to explore the dynamic world of gene expression? Start with the Jupyter notebook to begin your journey into transcriptomic data analysis!* 