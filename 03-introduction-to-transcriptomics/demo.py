#!/usr/bin/env python3
"""
Transcriptomic Data Analysis Demonstration

This script demonstrates the key capabilities of the transcriptomic analysis module,
including differential expression analysis, pathway enrichment, and visualization.

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

from transcriptomic_utils import (
    TranscriptomeDataProcessor,
    DifferentialExpressionAnalyzer,
    PathwayEnrichmentAnalyzer,
    TranscriptomeVisualizer
)
from config import (
    SAMPLE_METADATA,
    PATHWAY_GENE_SETS,
    DEMO_GENES,
    TREATMENT_INFO,
    generate_sample_count_matrix,
    get_sample_metadata_df,
    get_pathway_gene_sets,
    get_treatment_contrasts,
    get_demo_gene_info
)


def print_header(title):
    """Print a formatted header."""
    print("\n" + "="*60)
    print(f" {title}")
    print("="*60)


def print_subheader(title):
    """Print a formatted subheader."""
    print(f"\n--- {title} ---")


def demonstrate_data_loading():
    """Demonstrate RNA-seq data loading and processing."""
    print_header("RNA-SEQ DATA LOADING AND PROCESSING")
    
    # Initialize processor
    processor = TranscriptomeDataProcessor()
    
    print("Generating sample RNA-seq count matrix...")
    count_matrix = generate_sample_count_matrix(n_genes=2000, n_samples=12)
    
    print(f"Count matrix shape: {count_matrix.shape}")
    print(f"Genes: {count_matrix.shape[0]:,}")
    print(f"Samples: {count_matrix.shape[1]}")
    
    # Load metadata
    metadata = get_sample_metadata_df()
    print(f"\nSample metadata shape: {metadata.shape}")
    print("Sample conditions:", metadata['condition'].value_counts().to_dict())
    
    # Basic quality metrics
    print_subheader("Quality Control Metrics")
    
    # Library sizes
    library_sizes = processor.calculate_library_sizes(count_matrix)
    print(f"Library sizes (total counts per sample):")
    print(f"  Mean: {library_sizes.mean():,.0f}")
    print(f"  Range: {library_sizes.min():,.0f} - {library_sizes.max():,.0f}")
    
    # Detected genes per sample
    detected_genes = (count_matrix > 0).sum(axis=0)
    print(f"\nDetected genes per sample:")
    print(f"  Mean: {detected_genes.mean():.0f}")
    print(f"  Range: {detected_genes.min()} - {detected_genes.max()}")
    
    # Filter low expression genes
    filtered_matrix = processor.filter_low_expression_genes(
        count_matrix, min_count=5, min_samples=3
    )
    
    # Normalize data
    normalized_data = processor.normalize_counts(filtered_matrix, method='log_cpm')
    print(f"Normalized data shape: {normalized_data.shape}")
    
    return count_matrix, filtered_matrix, normalized_data, metadata


def demonstrate_differential_expression():
    """Demonstrate differential expression analysis."""
    print_header("DIFFERENTIAL EXPRESSION ANALYSIS")
    
    # Generate data
    count_matrix = generate_sample_count_matrix(n_genes=2000, n_samples=12)
    metadata = get_sample_metadata_df()
    
    # Initialize analyzer
    de_analyzer = DifferentialExpressionAnalyzer()
    processor = TranscriptomeDataProcessor()
    
    # Filter and normalize
    filtered_data = processor.filter_low_expression_genes(count_matrix, min_count=5, min_samples=3)
    normalized_data = processor.normalize_counts(filtered_data, method='log_cpm')
    
    # Perform differential expression analysis
    contrasts = get_treatment_contrasts()
    
    for treatment, control in contrasts[:2]:  # Analyze first two contrasts
        print_subheader(f"Contrast: {treatment} vs {control}")
        
        # Get sample lists
        treatment_samples = [s for s, info in SAMPLE_METADATA.items() 
                           if info['condition'] == treatment and s in normalized_data.columns]
        control_samples = [s for s, info in SAMPLE_METADATA.items() 
                         if info['condition'] == control and s in normalized_data.columns]
        
        print(f"Treatment samples ({treatment}): {treatment_samples}")
        print(f"Control samples ({control}): {control_samples}")
        
        if len(treatment_samples) >= 2 and len(control_samples) >= 2:
            # Perform DESeq2-like analysis
            de_results = de_analyzer.perform_deseq2_like_analysis(
                filtered_data, treatment_samples, control_samples
            )
            
            if len(de_results) > 0:
                print(f"Analyzed {len(de_results)} genes")
                
                # Identify significant genes
                significant_genes = de_analyzer.identify_top_genes(
                    de_results, p_threshold=0.05, fc_threshold=1.5, n_top=20
                )
                
                print(f"Significant DEGs: {len(significant_genes)}")
                
                if len(significant_genes) > 0:
                    print("\nTop 10 differentially expressed genes:")
                    top_genes = significant_genes.head(10)
                    for _, gene_row in top_genes.iterrows():
                        gene_name = gene_row['gene']
                        log2fc = gene_row['log2FoldChange']
                        padj = gene_row['padj']
                        direction = "↑" if log2fc > 0 else "↓"
                        print(f"  {gene_name}: {direction} {abs(log2fc):.2f} (padj={padj:.2e})")
                
                # Analyze known hormone-responsive genes
                print_subheader("Analysis of Known Hormone-Responsive Genes")
                hormone_genes = list(DEMO_GENES['hormone_responsive'].keys())
                
                for gene in hormone_genes:
                    if gene in de_results['gene'].values:
                        gene_result = de_results[de_results['gene'] == gene].iloc[0]
                        gene_info = get_demo_gene_info(gene)
                        expected = gene_info['expected_response'].get(treatment, 'unknown')
                        
                        log2fc = gene_result['log2FoldChange']
                        padj = gene_result['padj']
                        observed = "up" if log2fc > 0.5 else "down" if log2fc < -0.5 else "stable"
                        
                        match = "✓" if expected == observed else "✗" if expected != 'unknown' else "?"
                        
                        print(f"  {gene}: {log2fc:+.2f} (padj={padj:.3f}) "
                              f"Expected: {expected}, Observed: {observed} {match}")
            
            else:
                print("No genes passed filtering criteria")
        else:
            print("Insufficient samples for analysis")
    
    return de_results if 'de_results' in locals() else None


def demonstrate_pathway_analysis():
    """Demonstrate pathway enrichment analysis."""
    print_header("PATHWAY ENRICHMENT ANALYSIS")
    
    # Initialize analyzer
    pathway_analyzer = PathwayEnrichmentAnalyzer()
    pathway_analyzer.load_gene_sets(get_pathway_gene_sets())
    
    # Generate some example differentially expressed genes
    print("Simulating differentially expressed genes...")
    
    # Create a realistic gene list based on hormone treatment
    upregulated_genes = [
        'ESR1', 'PGR', 'TFF1', 'GREB1',  # Hormone responsive
        'CCND1', 'MKI67', 'PCNA',  # Cell cycle
        'MYC', 'JUN', 'FOS'  # Transcription factors
    ]
    
    downregulated_genes = [
        'CDKN1A', 'TP53',  # Cell cycle inhibitors
        'BAX', 'CASP3',  # Apoptosis
        'IL1B', 'TNF'  # Inflammatory
    ]
    
    all_de_genes = upregulated_genes + downregulated_genes
    
    print(f"Analyzing {len(all_de_genes)} differentially expressed genes:")
    print(f"  Upregulated: {upregulated_genes}")
    print(f"  Downregulated: {downregulated_genes}")
    
    # Perform overrepresentation analysis
    print_subheader("Gene Set Overrepresentation Analysis")
    
    enrichment_results = pathway_analyzer.perform_overrepresentation_analysis(all_de_genes)
    
    if len(enrichment_results) > 0:
        print("Significantly enriched pathways:")
        significant_pathways = enrichment_results[enrichment_results['p_adj'] < 0.05]
        
        if len(significant_pathways) > 0:
            for _, pathway_row in significant_pathways.head(10).iterrows():
                pathway_name = pathway_row['pathway']
                p_value = pathway_row['p_value']
                p_adj = pathway_row['p_adj']
                overlap_size = pathway_row['overlap_size']
                pathway_size = pathway_row['pathway_size']
                overlap_genes = pathway_row['overlap_genes']
                
                print(f"\n  {pathway_name}:")
                print(f"    P-value: {p_value:.2e} (adjusted: {p_adj:.2e})")
                print(f"    Overlap: {overlap_size}/{pathway_size} genes")
                print(f"    Genes: {', '.join(overlap_genes)}")
        else:
            print("No pathways reached significance threshold (p_adj < 0.05)")
    else:
        print("No pathway enrichment results found")
    
    # Demonstrate GSEA-like analysis
    print_subheader("Gene Set Enrichment Analysis (GSEA)")
    
    # Create ranked gene list (simulate log2 fold changes)
    np.random.seed(42)
    all_genes = []
    for pathway_info in PATHWAY_GENE_SETS.values():
        all_genes.extend(pathway_info['genes'])
    all_genes = list(set(all_genes))  # Remove duplicates
    
    # Assign fold changes
    fold_changes = {}
    for gene in all_genes:
        if gene in upregulated_genes:
            fold_changes[gene] = np.random.uniform(1.0, 3.0)
        elif gene in downregulated_genes:
            fold_changes[gene] = np.random.uniform(-3.0, -1.0)
        else:
            fold_changes[gene] = np.random.normal(0, 0.5)
    
    # Create ranked series
    ranked_genes = pd.Series(fold_changes).sort_values(ascending=False)
    
    print(f"Performing GSEA on {len(ranked_genes)} ranked genes...")
    
    gsea_results = pathway_analyzer.perform_gsea(ranked_genes)
    
    if len(gsea_results) > 0:
        print("GSEA results (top pathways):")
        for _, gsea_row in gsea_results.head(5).iterrows():
            pathway_name = gsea_row['pathway']
            es_score = gsea_row['es_score']
            p_value = gsea_row['p_value']
            
            direction = "Enriched in upregulated" if es_score > 0 else "Enriched in downregulated"
            print(f"  {pathway_name}: ES={es_score:.3f}, p={p_value:.3f} ({direction})")
    
    return enrichment_results if 'enrichment_results' in locals() else None


def demonstrate_visualization():
    """Demonstrate transcriptomic data visualization."""
    print_header("TRANSCRIPTOMIC DATA VISUALIZATION")
    
    # Generate data for visualization
    count_matrix = generate_sample_count_matrix(n_genes=2000, n_samples=12)
    metadata = get_sample_metadata_df()
    
    processor = TranscriptomeDataProcessor()
    filtered_data = processor.filter_low_expression_genes(count_matrix, min_count=5, min_samples=3)
    normalized_data = processor.normalize_counts(filtered_data, method='log_cpm')
    
    # Initialize visualizer
    visualizer = TranscriptomeVisualizer()
    
    print_subheader("Creating Visualizations")
    
    try:
        # PCA plot
        print("Creating PCA plot...")
        fig1 = visualizer.plot_pca(normalized_data, metadata=metadata, color_by='condition')
        plt.savefig('transcriptomics_pca.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("✓ Saved PCA plot")
        
        # Expression heatmap
        print("Creating expression heatmap...")
        fig2 = visualizer.plot_heatmap(normalized_data, sample_metadata=metadata, top_n=50)
        plt.savefig('transcriptomics_heatmap.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("✓ Saved expression heatmap")
        
        # Generate mock DE results for volcano plot
        de_analyzer = DifferentialExpressionAnalyzer()
        treatment_samples = [s for s, info in SAMPLE_METADATA.items() 
                           if info['condition'] == 'OHT' and s in normalized_data.columns]
        control_samples = [s for s, info in SAMPLE_METADATA.items() 
                         if info['condition'] == 'Control' and s in normalized_data.columns]
        
        if len(treatment_samples) >= 2 and len(control_samples) >= 2:
            de_results = de_analyzer.perform_deseq2_like_analysis(
                filtered_data, treatment_samples, control_samples
            )
            
            if len(de_results) > 0:
                # Volcano plot
                print("Creating volcano plot...")
                fig3 = visualizer.plot_volcano(de_results, p_col='padj', fc_col='log2FoldChange')
                plt.savefig('transcriptomics_volcano.png', dpi=300, bbox_inches='tight')
                plt.close()
                print("✓ Saved volcano plot")
        
        # Gene expression boxplot
        print("Creating gene expression boxplots...")
        demo_genes = ['ESR1', 'PGR', 'MKI67']
        available_genes = [g for g in demo_genes if g in normalized_data.index]
        
        if available_genes:
            fig4 = visualizer.plot_expression_boxplot(
                normalized_data, available_genes, metadata=metadata, group_by='condition'
            )
            plt.savefig('transcriptomics_boxplot.png', dpi=300, bbox_inches='tight')
            plt.close()
            print("✓ Saved gene expression boxplots")
        
        # Pathway enrichment plot (mock data)
        print("Creating pathway enrichment plot...")
        mock_enrichment = pd.DataFrame({
            'pathway': ['Hormone Response', 'Cell Cycle', 'Apoptosis', 'Metabolism', 'DNA Repair'],
            'p_value': [0.001, 0.005, 0.01, 0.02, 0.03],
            'overlap_size': [8, 12, 6, 15, 9],
            'pathway_size': [20, 50, 25, 80, 30]
        })
        
        fig5 = visualizer.plot_pathway_enrichment(mock_enrichment)
        plt.savefig('transcriptomics_pathways.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("✓ Saved pathway enrichment plot")
        
    except Exception as e:
        print(f"Visualization error: {e}")
        print("Note: Some visualizations may require additional dependencies")


def demonstrate_biological_insights():
    """Demonstrate biological insights from transcriptomic analysis."""
    print_header("BIOLOGICAL INSIGHTS FROM TRANSCRIPTOMIC ANALYSIS")
    
    print_subheader("Hormone Response Mechanisms")
    print("""
Key insights from hormone treatment transcriptomic analysis:

1. ESTROGEN RECEPTOR SIGNALING:
   - ESR1 (Estrogen Receptor 1) shows upregulation with both OHT and DPN
   - Classical estrogen-responsive genes (TFF1, GREB1) are activated
   - Demonstrates preserved estrogen signaling pathway functionality

2. SELECTIVE ESTROGEN RECEPTOR MODULATION:
   - OHT acts as both agonist and antagonist depending on tissue context
   - DPN shows ERβ-selective activation patterns
   - Different receptor subtypes mediate distinct transcriptional programs

3. CELL CYCLE REGULATION:
   - Hormone treatments affect proliferation markers (MKI67, CCND1)
   - Cell cycle checkpoint genes (CDKN1A) show treatment-specific responses
   - Balance between growth promotion and growth inhibition
    """)
    
    print_subheader("Clinical Applications")
    print("""
Transcriptomic analysis applications in precision medicine:

1. BIOMARKER DISCOVERY:
   - Expression signatures predict treatment response
   - Pathway activity scores indicate drug mechanism
   - Resistance mechanisms identified through expression changes

2. DRUG DEVELOPMENT:
   - Target identification through differential expression
   - Mechanism of action studies via pathway analysis
   - Biomarker-driven patient stratification

3. PERSONALIZED THERAPY:
   - Expression-based treatment selection
   - Monitoring therapeutic response
   - Early detection of resistance development
    """)
    
    print_subheader("Pathway-Level Understanding")
    print("""
Systems-level insights from pathway enrichment:

1. HORMONAL PATHWAYS:
   - Estrogen signaling cascade activation
   - Cross-talk with growth factor pathways
   - Metabolic reprogramming in response to hormones

2. CELL FATE DECISIONS:
   - Proliferation vs. differentiation balance
   - Apoptosis pathway modulation
   - Stress response activation

3. THERAPEUTIC IMPLICATIONS:
   - Combination therapy opportunities
   - Resistance pathway identification
   - Biomarker development for patient selection
    """)


def main():
    """Main demonstration function."""
    print("="*60)
    print(" TRANSCRIPTOMIC DATA ANALYSIS DEMONSTRATION")
    print(" Stanford Data Ocean - Data Science in Precision Medicine")
    print("="*60)
    
    print("""
This demonstration showcases the capabilities of the transcriptomic analysis module,
including:

• RNA-seq data processing and quality control
• Differential expression analysis (DESeq2-like methods)
• Pathway enrichment and functional analysis
• Comprehensive visualization tools
• Biological interpretation and clinical insights

The analysis uses simulated data based on real hormone treatment studies,
demonstrating estrogen receptor modulator effects on gene expression.
    """)
    
    try:
        # Run demonstrations
        count_data, filtered_data, normalized_data, metadata = demonstrate_data_loading()
        de_results = demonstrate_differential_expression()
        enrichment_results = demonstrate_pathway_analysis()
        demonstrate_visualization()
        demonstrate_biological_insights()
        
        print_header("DEMONSTRATION COMPLETE")
        print("""
Summary of generated files:
• transcriptomics_pca.png - Principal component analysis of samples
• transcriptomics_heatmap.png - Expression heatmap of top variable genes
• transcriptomics_volcano.png - Volcano plot of differential expression
• transcriptomics_boxplot.png - Gene expression boxplots by condition
• transcriptomics_pathways.png - Pathway enrichment analysis

Key findings demonstrated:
• Quality control and normalization of RNA-seq data
• Identification of differentially expressed genes
• Pathway-level analysis of biological processes
• Visualization of complex transcriptomic patterns

Next steps:
1. Explore the Jupyter notebook for interactive analysis
2. Apply methods to your own RNA-seq datasets
3. Customize pathway databases for specific research questions
4. Integrate with other omics data types

For more information, see the README.md file and comprehensive documentation.
        """)
        
    except Exception as e:
        print(f"\nError during demonstration: {e}")
        print("Please check your Python environment and dependencies.")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main()) 