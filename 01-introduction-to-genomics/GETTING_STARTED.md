# Getting Started with Genomic Data Analysis

Welcome to the **Introduction to Genomic Data Analysis** module! This guide will help you get up and running with the genomic analysis tools.

## 🚀 Quick Start

### 1. Install Dependencies

```bash
cd 01-introduction-to-genomics
pip install -r requirements.txt
```

### 2. Run the Demonstration

```bash
python demo.py
```

This will showcase all the key features without requiring AWS access.

### 3. Explore the Jupyter Notebook

```bash
jupyter lab notebooks/genomic_data_analysis.ipynb
```

**Note**: The notebook requires AWS Athena access for live data analysis.

---

## 📁 Module Structure

```
01-introduction-to-genomics/
├── scripts/
│   ├── genomic_utils.py      # Core analysis classes
│   └── config.py             # Configuration and constants
├── notebooks/
│   └── genomic_data_analysis.ipynb  # Interactive analysis
├── demo.py                   # Standalone demonstration
├── requirements.txt          # Python dependencies
├── setup.py                  # Module installation
├── README.md                 # Detailed documentation
└── GETTING_STARTED.md        # This file
```

---

## 🛠️ Core Components

### GenomicDataProcessor
- Extract population frequencies from VCF data
- Analyze variant quality scores
- Filter high-quality variants
- Generate chromosome summaries

### VariantAnnotator
- Classify variant types (SNV, Insertion, Deletion)
- Identify transitions vs transversions
- Calculate Ti/Tv ratios

### GenomicVisualizer
- Create population frequency plots
- Visualize chromosome distributions
- Plot quality score distributions

---

## 🧬 Example Usage

### Basic Population Analysis

```python
from scripts.genomic_utils import GenomicDataProcessor
from scripts.config import get_demo_snp_info

# Initialize processor
processor = GenomicDataProcessor()

# Analyze eye color SNP
info_text = "EAS_AF=0.002;EUR_AF=0.6362;AFR_AF=0.0023;AMR_AF=0.2017;SAS_AF=0.028"
frequencies = processor.extract_population_frequencies(info_text)

# Visualize results
processor.plot_population_frequencies(frequencies, "Eye Color SNP")
```

### Variant Quality Analysis

```python
import pandas as pd

# Sample genomic data
data = pd.DataFrame({
    'chrm': ['1', '2', '3'],
    'qual': [45.2, 67.8, 89.1],
    'reference_bases': ['A', 'C', 'G'],
    'alternate_bases': ['G', 'T', 'A']
})

# Analyze quality
stats = processor.analyze_variant_quality(data)
high_qual = processor.filter_high_quality_variants(data, min_quality=50.0)
```

---

## 🔗 AWS Athena Setup (Optional)

For live data analysis, you'll need AWS access:

### 1. Configure AWS Credentials

```bash
pip install awscli
aws configure
```

### 2. Set Environment Variables

```bash
export AWS_S3_STAGING_DIR=s3://your-bucket/
export AWS_DEFAULT_REGION=us-east-1
```

### 3. Available Datasets

- **1000 Genomes Project**: `default.g1000vcf_csv_int`
- **COSMIC Cancer Database**: `1000_genomes.hg19_cosmic68_int`
- **UCSC RefGene**: `1000_genomes.hg19_ucsc_refgene_int`

---

## 📊 Demo SNPs

The module includes several interesting SNPs for demonstration:

| SNP ID | Name | Description |
|--------|------|-------------|
| rs12913832 | Eye Color | Blue vs brown eye color genetics |
| rs53576 | Empathy | Social behavior and empathy traits |
| rs72921001 | Cilantro Taste | Love it or hate it taste preference |
| rs28936679 | Sleep Patterns | Circadian rhythm and sleep duration |
| rs1805009 | Skin Cancer | Melanoma and skin cancer risk |

---

## 🎯 Learning Path

1. **Start with the demo**: `python demo.py`
2. **Explore the utilities**: Review `scripts/genomic_utils.py`
3. **Try the notebook**: Work through the Jupyter notebook
4. **Read the documentation**: Check the main README for detailed analysis examples
5. **Experiment**: Try analyzing different SNPs and populations

---

## 🔧 Troubleshooting

### Common Issues

**Import Error**
```
ImportError: No module named 'genomic_utils'
```
**Solution**: Make sure you're running from the module directory:
```bash
cd 01-introduction-to-genomics
python demo.py
```

**Missing Dependencies**
```
ModuleNotFoundError: No module named 'pandas'
```
**Solution**: Install requirements:
```bash
pip install -r requirements.txt
```

**AWS Connection Issues**
```
Error: Failed to connect to Athena
```
**Solution**: The demo works without AWS. For live analysis, configure AWS credentials.

---

## 🤝 Need Help?

- Check the main module README for detailed documentation
- Review the demo script for working examples
- Examine the utility classes for available methods
- Explore the configuration file for available datasets

---

*This module demonstrates practical genomic data analysis techniques as part of the Stanford Data Ocean Precision Medicine certificate program.* 