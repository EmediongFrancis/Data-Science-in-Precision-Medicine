#!/usr/bin/env python3
"""
Epigenetic Data Analysis Demonstration

This script demonstrates the key capabilities of the epigenetic analysis module,
including DNA methylation analysis, histone modification analysis, and 
chromatin state classification.

Usage:
    python demo.py

Author: Stanford Data Ocean - Data Science in Precision Medicine
"""

import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Add scripts directory to path
sys.path.append(str(Path(__file__).parent / 'scripts'))

from epigenetic_utils import (
    EpigenomeDataProcessor, 
    MethylationAnalyzer, 
    HistoneModificationAnalyzer,
    EpigenomeVisualizer
)
from config import (
    SAMPLE_METADATA, 
    HISTONE_MARKS, 
    DEMO_GENES,
    generate_sample_methylation_data,
    generate_sample_histone_data,
    get_gene_info,
    get_histone_mark_info
)


def print_header(title):
    """Print a formatted header."""
    print("\n" + "="*60)
    print(f" {title}")
    print("="*60)


def print_subheader(title):
    """Print a formatted subheader."""
    print(f"\n--- {title} ---")


def demonstrate_methylation_analysis():
    """Demonstrate DNA methylation analysis capabilities."""
    print_header("DNA METHYLATION ANALYSIS")
    
    # Initialize processors
    processor = EpigenomeDataProcessor()
    analyzer = MethylationAnalyzer()
    visualizer = EpigenomeVisualizer()
    
    print("Generating sample methylation data...")
    methylation_data = generate_sample_methylation_data(n_regions=500, n_samples=20)
    
    print(f"Dataset shape: {methylation_data.shape}")
    print(f"Chromosomes: {sorted(methylation_data['chr'].unique())}")
    
    # Basic statistics
    print_subheader("Methylation Statistics")
    stats = processor.calculate_methylation_statistics(methylation_data)
    print(f"Total regions: {stats['total_regions']}")
    print(f"Total samples: {stats['total_samples']}")
    print(f"Methylation range: {stats['methylation_range']['min']:.3f} - {stats['methylation_range']['max']:.3f}")
    
    # Sample means
    sample_means = pd.Series(stats['mean_methylation_per_sample'])
    print(f"Mean methylation across samples: {sample_means.mean():.3f} ± {sample_means.std():.3f}")
    
    # Methylation classification
    print_subheader("Methylation Classification")
    sample_col = [col for col in methylation_data.columns 
                  if col not in ['chr', 'start', 'end', 'region_id']][0]
    
    classifications = analyzer.classify_methylation_levels(methylation_data[sample_col])
    class_counts = pd.Series(classifications).value_counts()
    
    print("Methylation level distribution:")
    for category, count in class_counts.items():
        percentage = (count / len(classifications)) * 100
        print(f"  {category}: {count} regions ({percentage:.1f}%)")
    
    # Differential methylation analysis
    print_subheader("Differential Methylation Analysis")
    
    # Compare ESC vs ESC-derived samples
    esc_samples = [sample for sample, info in SAMPLE_METADATA.items() 
                   if info['category'] == 'ESC' and sample in methylation_data.columns]
    esc_derived_samples = [sample for sample, info in SAMPLE_METADATA.items() 
                          if info['category'] == 'ESC_Derived' and sample in methylation_data.columns]
    
    if len(esc_samples) >= 2 and len(esc_derived_samples) >= 2:
        dmrs = processor.identify_differential_regions(
            methylation_data, esc_samples, esc_derived_samples,
            min_diff=0.15, p_threshold=0.05
        )
        
        if len(dmrs) > 0:
            print(f"Found {len(dmrs)} differentially methylated regions")
            print("Top 5 DMRs:")
            print(dmrs[['region_id', 'methylation_diff', 'p_value']].head())
        else:
            print("No significant DMRs found with current thresholds")
    
    # Visualization
    print_subheader("Creating Visualizations")
    
    # Filter for most variable regions for visualization
    top_variable = processor.filter_high_variance_regions(methylation_data, top_n=100)
    
    try:
        # Methylation distribution plot
        fig1 = visualizer.plot_methylation_distribution(methylation_data)
        plt.savefig('methylation_distribution.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("✓ Saved methylation distribution plot")
        
        # PCA analysis
        fig2 = visualizer.plot_pca_analysis(top_variable)
        plt.savefig('methylation_pca.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("✓ Saved PCA analysis plot")
        
    except Exception as e:
        print(f"Visualization error: {e}")
    
    return methylation_data


def demonstrate_histone_analysis():
    """Demonstrate histone modification analysis."""
    print_header("HISTONE MODIFICATION ANALYSIS")
    
    # Initialize analyzer
    histone_analyzer = HistoneModificationAnalyzer()
    visualizer = EpigenomeVisualizer()
    
    print("Generating sample histone modification data...")
    histone_data = generate_sample_histone_data(n_regions=50, window_size=100)
    
    print(f"Dataset shape: {histone_data.shape}")
    print(f"Histone marks: {[col for col in histone_data.columns if col.startswith('H3')]}")
    
    # Histone mark information
    print_subheader("Histone Modification Marks")
    for mark_name, mark_info in HISTONE_MARKS.items():
        if mark_name in histone_data.columns:
            mean_signal = histone_data[mark_name].mean()
            print(f"{mark_name}: {mark_info['function']} (mean signal: {mean_signal:.2f})")
    
    # Chromatin state classification
    print_subheader("Chromatin State Classification")
    
    # Get data for first gene as example
    gene_data = histone_data[histone_data['gene'] == 'Gene_1'].copy()
    if len(gene_data) > 0:
        # Average across positions for each gene
        gene_summary = histone_data.groupby('gene')[['H3K4me3', 'H3K4me1', 'H3K27ac', 'H3K27me3', 'H3K9me3']].mean()
        
        # Classify chromatin states
        classified_data = histone_analyzer.classify_chromatin_states(gene_summary)
        
        # Count chromatin states
        state_counts = classified_data['chromatin_state'].value_counts()
        print("Chromatin state distribution:")
        for state, count in state_counts.items():
            percentage = (count / len(classified_data)) * 100
            print(f"  {state}: {count} genes ({percentage:.1f}%)")
    
    # TSS enrichment analysis (simulated)
    print_subheader("TSS Enrichment Analysis")
    
    # Create mock TSS positions
    tss_positions = pd.DataFrame({
        'gene_name': [f'Gene_{i}' for i in range(1, 11)],
        'chr': ['chr1'] * 10,
        'tss_position': np.random.randint(1000000, 2000000, 10)
    })
    
    # Mock ChIP-seq data
    chip_data = pd.DataFrame({
        'chr': ['chr1'] * 100,
        'start': np.random.randint(1000000, 2000000, 100),
        'end': np.random.randint(1000000, 2000000, 100),
        'signal': np.random.exponential(2, 100)
    })
    chip_data['end'] = chip_data['start'] + np.random.randint(200, 1000, 100)
    
    try:
        enrichment_results = histone_analyzer.calculate_tss_enrichment(
            chip_data, tss_positions, window_size=2000
        )
        
        print("TSS enrichment scores:")
        print(enrichment_results[['gene', 'enrichment_score', 'peak_count']].head())
        
    except Exception as e:
        print(f"TSS enrichment analysis error: {e}")
    
    # Visualization
    print_subheader("Creating Visualizations")
    
    try:
        # Histone modification profile
        gene_profile = histone_data[histone_data['gene'] == 'Gene_1']
        if len(gene_profile) > 0:
            fig1 = visualizer.plot_histone_modification_profile(
                gene_profile.set_index('position'), 
                gene_name='Gene_1'
            )
            plt.savefig('histone_profile.png', dpi=300, bbox_inches='tight')
            plt.close()
            print("✓ Saved histone modification profile")
        
        # Chromatin state pie chart
        if 'classified_data' in locals() and len(classified_data) > 0:
            fig2 = visualizer.plot_chromatin_state_pie(classified_data['chromatin_state'])
            plt.savefig('chromatin_states.png', dpi=300, bbox_inches='tight')
            plt.close()
            print("✓ Saved chromatin state distribution")
        
    except Exception as e:
        print(f"Visualization error: {e}")
    
    return histone_data


def demonstrate_gene_analysis():
    """Demonstrate gene-specific epigenetic analysis."""
    print_header("GENE-SPECIFIC EPIGENETIC ANALYSIS")
    
    # Analyze key developmental genes
    print_subheader("Pluripotency Genes")
    
    for gene_name, gene_info in DEMO_GENES['pluripotency'].items():
        print(f"\n{gene_info['name']}:")
        print(f"  Location: {gene_info['chr']}:{gene_info['start']}-{gene_info['end']}")
        print(f"  Function: {gene_info['function']}")
        
        # Simulate epigenetic analysis for this gene
        region_length = gene_info['end'] - gene_info['start']
        
        # Simulate methylation levels (pluripotency genes typically have low methylation in ESCs)
        esc_methylation = np.random.beta(1, 5)  # Low methylation
        differentiated_methylation = np.random.beta(3, 2)  # Higher methylation
        
        print(f"  Simulated methylation - ESC: {esc_methylation:.3f}, Differentiated: {differentiated_methylation:.3f}")
        
        # Simulate histone marks
        h3k4me3_level = np.random.uniform(2, 5)  # High active mark
        h3k27me3_level = np.random.uniform(0, 1)  # Low repressive mark
        
        print(f"  Simulated H3K4me3: {h3k4me3_level:.2f}, H3K27me3: {h3k27me3_level:.2f}")
    
    print_subheader("Cancer-Related Genes")
    
    for gene_name, gene_info in DEMO_GENES['cancer'].items():
        print(f"\n{gene_info['name']}:")
        print(f"  Location: {gene_info['chr']}:{gene_info['start']}-{gene_info['end']}")
        print(f"  Function: {gene_info['function']}")
        
        # Tumor suppressors often become hypermethylated in cancer
        normal_methylation = np.random.beta(1, 4)  # Low methylation
        cancer_methylation = np.random.beta(4, 1)  # High methylation
        
        print(f"  Simulated methylation - Normal: {normal_methylation:.3f}, Cancer: {cancer_methylation:.3f}")


def demonstrate_integration_analysis():
    """Demonstrate integrated epigenomic analysis."""
    print_header("INTEGRATED EPIGENOMIC ANALYSIS")
    
    print("Integrating DNA methylation and histone modification data...")
    
    # Generate coordinated data
    np.random.seed(42)
    n_regions = 200
    
    # Create regions with different epigenetic states
    regions = []
    for i in range(n_regions):
        region_id = f"Region_{i+1}"
        
        # Assign chromatin state
        state = np.random.choice([
            'Active_Promoter', 'Strong_Enhancer', 'Polycomb_Repressed', 
            'Heterochromatin', 'Quiescent'
        ], p=[0.15, 0.20, 0.15, 0.10, 0.40])
        
        # Generate consistent methylation and histone patterns
        if state == 'Active_Promoter':
            methylation = np.random.beta(1, 5)  # Low
            h3k4me3 = np.random.uniform(2, 5)   # High
            h3k27me3 = np.random.uniform(0, 0.5)  # Low
        elif state == 'Strong_Enhancer':
            methylation = np.random.beta(1, 4)  # Low
            h3k4me3 = np.random.uniform(0, 1)   # Low
            h3k27me3 = np.random.uniform(0, 0.5)  # Low
        elif state == 'Polycomb_Repressed':
            methylation = np.random.beta(2, 3)  # Variable
            h3k4me3 = np.random.uniform(0, 1)   # Low
            h3k27me3 = np.random.uniform(1, 3)  # High
        elif state == 'Heterochromatin':
            methylation = np.random.beta(4, 1)  # High
            h3k4me3 = np.random.uniform(0, 0.5) # Very low
            h3k27me3 = np.random.uniform(0, 1)  # Variable
        else:  # Quiescent
            methylation = np.random.beta(2, 2)  # Variable
            h3k4me3 = np.random.uniform(0, 1)   # Low
            h3k27me3 = np.random.uniform(0, 1)  # Low
        
        regions.append({
            'region_id': region_id,
            'chromatin_state': state,
            'methylation': methylation,
            'H3K4me3': h3k4me3,
            'H3K27me3': h3k27me3
        })
    
    integrated_data = pd.DataFrame(regions)
    
    print_subheader("Integrated Analysis Results")
    
    # Correlation analysis
    correlation_matrix = integrated_data[['methylation', 'H3K4me3', 'H3K27me3']].corr()
    print("Correlation matrix:")
    print(correlation_matrix.round(3))
    
    # State-specific analysis
    print("\nMean values by chromatin state:")
    state_summary = integrated_data.groupby('chromatin_state')[['methylation', 'H3K4me3', 'H3K27me3']].mean()
    print(state_summary.round(3))
    
    # Visualization
    print_subheader("Creating Integrated Visualizations")
    
    try:
        # Correlation heatmap
        plt.figure(figsize=(8, 6))
        sns.heatmap(correlation_matrix, annot=True, cmap='RdBu_r', center=0)
        plt.title('Epigenetic Mark Correlations')
        plt.tight_layout()
        plt.savefig('epigenetic_correlations.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("✓ Saved correlation heatmap")
        
        # Scatter plot of methylation vs H3K4me3
        plt.figure(figsize=(10, 6))
        for state in integrated_data['chromatin_state'].unique():
            state_data = integrated_data[integrated_data['chromatin_state'] == state]
            plt.scatter(state_data['methylation'], state_data['H3K4me3'], 
                       label=state, alpha=0.7)
        
        plt.xlabel('DNA Methylation Level')
        plt.ylabel('H3K4me3 Signal')
        plt.title('Relationship between DNA Methylation and H3K4me3')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig('methylation_vs_h3k4me3.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("✓ Saved methylation vs H3K4me3 scatter plot")
        
    except Exception as e:
        print(f"Visualization error: {e}")
    
    return integrated_data


def demonstrate_biological_insights():
    """Demonstrate biological insights from epigenetic analysis."""
    print_header("BIOLOGICAL INSIGHTS FROM EPIGENETIC ANALYSIS")
    
    print_subheader("Development and Differentiation")
    print("""
Key insights from epigenetic analysis in development:

1. PLURIPOTENCY MAINTENANCE:
   - OCT4, NANOG, SOX2 promoters have low DNA methylation in ESCs
   - High H3K4me3 marks at pluripotency gene promoters
   - Bivalent domains (H3K4me3 + H3K27me3) mark poised developmental genes

2. CELL DIFFERENTIATION:
   - Progressive DNA methylation of pluripotency gene promoters
   - Loss of H3K4me3 and gain of H3K27me3 at silenced genes
   - Activation of lineage-specific enhancers (H3K27ac gain)

3. CHROMATIN STATE TRANSITIONS:
   - Active → Poised → Repressed states during differentiation
   - Enhancer switching drives cell-type specific gene expression
   - DNA methylation provides stable gene silencing
    """)
    
    print_subheader("Disease and Cancer")
    print("""
Epigenetic changes in cancer:

1. TUMOR SUPPRESSOR SILENCING:
   - Hypermethylation of CDKN2A, MLH1, BRCA1 promoters
   - Loss of H3K4me3 at silenced tumor suppressor genes
   - CpG island methylator phenotype (CIMP) in certain cancers

2. ONCOGENE ACTIVATION:
   - Hypomethylation of oncogene promoters
   - Enhancer activation through H3K27ac deposition
   - Chromatin remodeling complex mutations

3. THERAPEUTIC TARGETS:
   - DNA methyltransferase inhibitors (5-azacytidine)
   - Histone deacetylase inhibitors (SAHA)
   - EZH2 inhibitors for H3K27me3-driven cancers
    """)
    
    print_subheader("Environmental Responses")
    print("""
Epigenetic responses to environment:

1. LIFESTYLE FACTORS:
   - Diet affects DNA methylation patterns
   - Exercise modulates histone modifications
   - Stress alters chromatin accessibility

2. AGING:
   - Progressive DNA methylation changes
   - Histone mark redistribution
   - Chromatin compaction increases

3. TRANSGENERATIONAL EFFECTS:
   - Some epigenetic marks are heritable
   - Environmental exposures can affect offspring
   - Imprinting defects cause developmental disorders
    """)


def main():
    """Main demonstration function."""
    print("="*60)
    print(" EPIGENETIC DATA ANALYSIS DEMONSTRATION")
    print(" Stanford Data Ocean - Data Science in Precision Medicine")
    print("="*60)
    
    print("""
This demonstration showcases the capabilities of the epigenetic analysis module,
including:

• DNA methylation pattern analysis
• Histone modification profiling  
• Chromatin state classification
• Integrated epigenomic analysis
• Biological interpretation

All analyses use simulated data based on real epigenomic patterns from the
NIH Roadmap Epigenomics Program.
    """)
    
    try:
        # Run demonstrations
        methylation_data = demonstrate_methylation_analysis()
        histone_data = demonstrate_histone_analysis()
        demonstrate_gene_analysis()
        integrated_data = demonstrate_integration_analysis()
        demonstrate_biological_insights()
        
        print_header("DEMONSTRATION COMPLETE")
        print("""
Summary of generated files:
• methylation_distribution.png - Distribution of methylation levels
• methylation_pca.png - PCA analysis of methylation data
• histone_profile.png - Histone modification profile around TSS
• chromatin_states.png - Chromatin state distribution
• epigenetic_correlations.png - Correlation between epigenetic marks
• methylation_vs_h3k4me3.png - Relationship between methylation and H3K4me3

Next steps:
1. Explore the Jupyter notebook for interactive analysis
2. Try the analysis with your own epigenomic datasets
3. Modify parameters to explore different biological questions

For more information, see the README.md file and documentation.
        """)
        
    except Exception as e:
        print(f"\nError during demonstration: {e}")
        print("Please check your Python environment and dependencies.")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main()) 