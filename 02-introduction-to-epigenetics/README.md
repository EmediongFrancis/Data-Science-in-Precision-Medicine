# Introduction to Epigenetics

**Stanford Data Ocean - Data Science in Precision Medicine**

## Overview

Epigenetics represents one of the most dynamic and fascinating areas of modern genomics. Unlike the static DNA sequence we explored in Module 1, epigenetic modifications are dynamic changes that control gene expression without altering the underlying DNA sequence. These modifications play crucial roles in development, aging, disease, and responses to environmental factors.

This module provides hands-on experience with epigenomic data analysis, focusing on DNA methylation patterns and histone modifications across human reference epigenomes from the NIH Roadmap Epigenomics Program.

## Learning Objectives

By completing this module, you will:

- **Understand epigenetic mechanisms** - Learn how DNA methylation and histone modifications regulate gene expression
- **Analyze real epigenomic datasets** - Work with data from 111 reference human epigenomes
- **Explore development and disease** - Understand how epigenetic changes drive cell differentiation and contribute to disease
- **Master visualization techniques** - Create meaningful plots of epigenomic landscapes
- **Apply computational methods** - Use Python tools for epigenomic data analysis

## Module Structure

```
02-introduction-to-epigenetics/
├── scripts/
│   ├── epigenetic_utils.py          # Core analysis functions
│   ├── config.py                    # Configuration and constants
│   └── visualization_utils.py       # Plotting and visualization tools
├── notebooks/
│   └── epigenetic_data_analysis.ipynb  # Main analysis notebook
├── data/                            # Sample datasets
├── demo.py                          # Standalone demonstration script
├── requirements.txt                 # Python dependencies
├── setup.py                        # Package installation
└── README.md                        # This file
```

## Key Concepts Covered

### Epigenetic Mechanisms
- **DNA Methylation** - How cytosine methylation silences genes
- **Histone Modifications** - Chemical changes that affect chromatin structure
- **Chromatin States** - Active vs. repressive chromatin configurations

### Experimental Techniques
- **Whole Genome Bisulfite Sequencing (WGBS)** - Genome-wide methylation mapping
- **ChIP-Sequencing** - Chromatin immunoprecipitation for histone modifications
- **Differentially Methylated Regions (DMRs)** - Identifying regulatory regions

### Biological Applications
- **Development** - How epigenetics drives cell differentiation
- **Disease** - Epigenetic changes in cancer and other diseases
- **Environmental Response** - How lifestyle factors affect the epigenome

## Datasets

This module uses data from the **NIH Roadmap Epigenomics Program**, including:

- **111 reference human epigenomes** representing diverse cell types and tissues
- **Whole Genome Bisulfite Sequencing data** for DNA methylation analysis
- **ChIP-seq data** for five core histone modification marks
- **Differentially Methylated Regions** identified across all epigenomes

## Technical Requirements

### Python Dependencies
```
pandas >= 1.3.0
numpy >= 1.21.0
matplotlib >= 3.4.0
seaborn >= 0.11.0
scipy >= 1.7.0
scikit-learn >= 1.0.0
```

### Optional (for advanced analysis)
```
pybedtools >= 0.8.0
pysam >= 0.16.0
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
   jupyter notebook notebooks/epigenetic_data_analysis.ipynb
   ```

## Data Sources

- **NIH Roadmap Epigenomics Program**: https://egg2.wustl.edu/roadmap/web_portal/index.html
- **Reference**: Roadmap Epigenomics Consortium. Integrative analysis of 111 reference human epigenomes. *Nature* 518, 317–330 (2015)

## Key Analysis Workflows

### 1. DNA Methylation Analysis
- Load and process WGBS data from 111 epigenomes
- Identify differentially methylated regions (DMRs)
- Analyze methylation patterns across cell types
- Visualize methylation landscapes

### 2. Histone Modification Analysis
- Process ChIP-seq peak data for core histone marks
- Map modifications around transcription start sites
- Identify active vs. repressive chromatin states
- Compare modifications between cell types

### 3. Integrative Analysis
- Combine methylation and histone modification data
- Identify regulatory elements and their states
- Analyze gene expression relationships
- Explore development and disease patterns

## Example Applications

### Stem Cell Biology
Examine how epigenetic marks change during differentiation:
- **ESC markers**: High H3K4me3, low DNA methylation at pluripotency genes
- **Differentiation**: Progressive methylation and chromatin changes
- **Cell-type specific patterns**: Unique epigenomic signatures

### Disease Analysis
Investigate epigenetic changes in disease:
- **Cancer epigenomes**: Hypermethylation of tumor suppressors
- **Aging signatures**: Progressive methylation changes
- **Environmental responses**: Lifestyle-induced epigenetic modifications

## Visualization Capabilities

The module includes comprehensive visualization tools for:
- **Methylation heatmaps** across samples and regions
- **Chromatin state maps** showing histone modification patterns
- **Gene-centric views** of epigenetic landscapes
- **Comparative analysis plots** between cell types and conditions

## Educational Features

- **Progressive complexity** - Builds from basic concepts to advanced analysis
- **Real-world examples** - Uses actual research datasets and findings
- **Interactive exploration** - Hands-on analysis with immediate feedback
- **Biological context** - Connects computational results to biological meaning

## Further Reading

### Core Papers
- Roadmap Epigenomics Consortium (2015). Integrative analysis of 111 reference human epigenomes. *Nature*
- Ernst & Kellis (2012). ChromHMM: automating chromatin-state discovery and characterization. *Nature Methods*
- Lister et al. (2013). Global epigenomic reconfiguration during mammalian brain development. *Science*

### Methodological Resources
- ENCODE Project Guidelines for ChIP-seq analysis
- Bisulfite sequencing data analysis best practices
- Chromatin state annotation methods

## Support and Troubleshooting

### Common Issues
- **Memory requirements**: Large epigenomic datasets may require significant RAM
- **Data access**: Some analyses require downloading large reference datasets
- **Visualization**: Complex plots may take time to render

### Getting Help
- Check the troubleshooting section in the notebook
- Review the example outputs in the demo script
- Consult the documentation in the utility modules

## Contributions

This module was developed as part of the Stanford Data Ocean certificate program in Data Science for Precision Medicine.

**Module Development Team:**
- **Content Development**: Abtin Tondar, Mohan Babu
- **Technical Implementation**: Amit Dixit
- **Educational Design**: Multiple contributors
- **Data Curation**: NIH Roadmap Epigenomics Consortium

---

*Ready to explore the dynamic world of epigenetics? Start with the Jupyter notebook to begin your journey into epigenomic data analysis!* 