#!/usr/bin/env python3
"""
Proteomic Data Analysis Demonstration

This script demonstrates the key capabilities of the proteomic analysis module,
including protein abundance analysis, differential expression, pathway enrichment,
and clinical interpretation using iPOP study data.

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

from proteomic_utils import (
    ProteomeDataProcessor,
    DifferentialProteinAnalyzer,
    ProteinPathwayAnalyzer,
    ProteomicsVisualizer
)
from config import (
    SAMPLE_METADATA,
    PROTEIN_PATHWAYS,
    DEMO_PROTEINS,
    IPOP_STUDY_INFO,
    MS_PARAMETERS,
    generate_sample_abundance_matrix,
    get_sample_metadata_df,
    get_protein_pathways,
    get_analysis_contrasts,
    get_demo_protein_info
)


def print_header(title):
    """Print a formatted header."""
    print("\n" + "="*60)
    print(f" {title}")
    print("="*60)


def print_subheader(title):
    """Print a formatted subheader."""
    print(f"\n--- {title} ---")


def demonstrate_ipop_study_context():
    """Demonstrate the iPOP study context and proteomics background."""
    print_header("iPOP STUDY CONTEXT AND PROTEOMICS BACKGROUND")
    
    print("INTEGRATED PERSONAL OMICS PROFILING (iPOP) STUDY")
    print(f"Institution: {IPOP_STUDY_INFO['institution']}")
    print(f"Study Name: {IPOP_STUDY_INFO['name']}")
    print(f"Participants: {IPOP_STUDY_INFO['patient_count']} individuals")
    print(f"Focus: {IPOP_STUDY_INFO['focus']}")
    print(f"Description: {IPOP_STUDY_INFO['description']}")
    
    print_subheader("Mass Spectrometry Parameters")
    print(f"Instrument: {MS_PARAMETERS['instrument_type']}")
    print(f"Acquisition Mode: {MS_PARAMETERS['acquisition_mode']}")
    print(f"Mass Range: {MS_PARAMETERS['mass_range']}")
    print(f"Resolution: {MS_PARAMETERS['resolution']:,}")
    print(f"Quantification: {MS_PARAMETERS['quantification_method']}")
    
    print_subheader("Patient Z Case Study")
    print("""
Patient Z represents a detailed longitudinal case study from the iPOP cohort,
demonstrating how proteomics can track individual health status over time:

• BASELINE: Healthy state protein profile
• INFECTION: Acute phase response during viral infection
• RECOVERY: Gradual return to baseline protein levels
• CLINICAL CONTEXT: Real patient data with detailed monitoring

This case demonstrates the power of personalized proteomics for:
- Early disease detection
- Treatment monitoring
- Recovery assessment
- Precision medicine applications
    """)


def demonstrate_data_loading():
    """Demonstrate protein abundance data loading and processing."""
    print_header("PROTEIN ABUNDANCE DATA LOADING AND PROCESSING")
    
    # Initialize processor
    processor = ProteomeDataProcessor()
    
    print("Generating sample protein abundance matrix...")
    abundance_matrix = generate_sample_abundance_matrix(n_proteins=3000, n_samples=12)
    
    print(f"Abundance matrix shape: {abundance_matrix.shape}")
    print(f"Proteins: {abundance_matrix.shape[0]:,}")
    print(f"Samples: {abundance_matrix.shape[1]}")
    
    # Load metadata
    metadata = get_sample_metadata_df()
    print(f"\nSample metadata shape: {metadata.shape}")
    print("Sample conditions:", metadata['condition'].value_counts().to_dict())
    print("Infection status:", metadata['infection_status'].value_counts().to_dict())
    
    # Basic quality metrics
    print_subheader("Quality Control Metrics")
    
    # Missing value analysis
    missing_pct = abundance_matrix.isnull().sum().sum() / abundance_matrix.size
    print(f"Overall missing values: {missing_pct*100:.1f}%")
    
    # Protein statistics
    protein_stats = processor.calculate_protein_statistics(abundance_matrix)
    print(f"Protein abundance statistics:")
    print(f"  Mean abundance: {protein_stats['mean_abundance'].mean():.2f}")
    print(f"  Median CV: {protein_stats['cv'].median():.2f}")
    print(f"  Average detection frequency: {protein_stats['detection_freq'].mean():.2f}")
    
    # Handle missing values
    processed_data = processor.handle_missing_values(abundance_matrix, method='knn')
    print(f"Processed data shape: {processed_data.shape}")
    
    # Normalize data
    normalized_data = processor.normalize_abundance_data(processed_data, method='median_centering')
    print(f"Normalized data range: {normalized_data.min().min():.2f} to {normalized_data.max().max():.2f}")
    
    return abundance_matrix, processed_data, normalized_data, metadata


def demonstrate_differential_protein_analysis():
    """Demonstrate differential protein expression analysis."""
    print_header("DIFFERENTIAL PROTEIN EXPRESSION ANALYSIS")
    
    # Generate data
    abundance_matrix = generate_sample_abundance_matrix(n_proteins=3000, n_samples=12)
    metadata = get_sample_metadata_df()
    
    # Initialize analyzer
    de_analyzer = DifferentialProteinAnalyzer()
    processor = ProteomeDataProcessor()
    
    # Process data
    processed_data = processor.handle_missing_values(abundance_matrix, method='knn')
    normalized_data = processor.normalize_abundance_data(processed_data, method='median_centering')
    
    # Perform differential analysis
    contrasts = get_analysis_contrasts()
    
    for condition1, condition2 in contrasts[:2]:  # Analyze first two contrasts
        print_subheader(f"Contrast: {condition1} vs {condition2}")
        
        # Get sample lists
        group1_samples = [s for s, info in SAMPLE_METADATA.items() 
                         if info['condition'] == condition1 and s in normalized_data.columns]
        group2_samples = [s for s, info in SAMPLE_METADATA.items() 
                         if info['condition'] == condition2 and s in normalized_data.columns]
        
        print(f"Group 1 samples ({condition1}): {group1_samples}")
        print(f"Group 2 samples ({condition2}): {group2_samples}")
        
        if len(group1_samples) >= 2 and len(group2_samples) >= 2:
            # Perform t-test analysis
            de_results = de_analyzer.perform_ttest_analysis(
                normalized_data, group1_samples, group2_samples
            )
            
            if len(de_results) > 0:
                print(f"Analyzed {len(de_results)} proteins")
                
                # Identify significant proteins
                significant_proteins = de_analyzer.identify_significant_proteins(
                    de_results, p_threshold=0.05, fc_threshold=1.5
                )
                
                print(f"Significant proteins: {len(significant_proteins)}")
                
                if len(significant_proteins) > 0:
                    print("\nTop 10 differentially expressed proteins:")
                    top_proteins = significant_proteins.head(10)
                    for _, protein_row in top_proteins.iterrows():
                        protein_name = protein_row['protein']
                        log2fc = protein_row['log2_fold_change']
                        p_adj = protein_row['p_adj']
                        direction = "↑" if log2fc > 0 else "↓"
                        print(f"  {protein_name}: {direction} {abs(log2fc):.2f} (p_adj={p_adj:.2e})")
                
                # Analyze known acute phase proteins
                print_subheader("Analysis of Known Acute Phase Proteins")
                acute_phase_proteins = list(DEMO_PROTEINS['acute_phase'].keys())
                
                for protein in acute_phase_proteins:
                    if protein in de_results['protein'].values:
                        protein_result = de_results[de_results['protein'] == protein].iloc[0]
                        protein_info = get_demo_protein_info(protein)
                        expected = protein_info['expected_response'].get(condition1, 'unknown')
                        
                        log2fc = protein_result['log2_fold_change']
                        p_adj = protein_result['p_adj']
                        observed = "up" if log2fc > 0.5 else "down" if log2fc < -0.5 else "stable"
                        
                        match = "✓" if expected == observed else "✗" if expected != 'unknown' else "?"
                        
                        print(f"  {protein}: {log2fc:+.2f} (p_adj={p_adj:.3f}) "
                              f"Expected: {expected}, Observed: {observed} {match}")
            
            else:
                print("No proteins passed filtering criteria")
        else:
            print("Insufficient samples for analysis")
    
    return de_results if 'de_results' in locals() else None


def demonstrate_pathway_analysis():
    """Demonstrate protein pathway enrichment analysis."""
    print_header("PROTEIN PATHWAY ENRICHMENT ANALYSIS")
    
    # Initialize analyzer
    pathway_analyzer = ProteinPathwayAnalyzer()
    pathway_analyzer.load_pathway_databases(get_protein_pathways())
    
    # Generate example differentially expressed proteins based on infection response
    print("Simulating infection-responsive proteins...")
    
    # Create realistic protein lists based on acute phase response
    upregulated_proteins = [
        'CRP', 'SAA1', 'HP', 'FGA', 'FGB',  # Acute phase proteins
        'IGHG1', 'IGKC',  # Immunoglobulins
        'C3', 'C4A',  # Complement proteins
        'HSPA1A', 'HSP90AA1'  # Heat shock proteins
    ]
    
    downregulated_proteins = [
        'ALB', 'APOA1', 'APOB',  # Negative acute phase
        'TTR',  # Transthyretin
        'CFD'  # Complement factor D (initially)
    ]
    
    all_de_proteins = upregulated_proteins + downregulated_proteins
    
    print(f"Analyzing {len(all_de_proteins)} differentially expressed proteins:")
    print(f"  Upregulated: {upregulated_proteins}")
    print(f"  Downregulated: {downregulated_proteins}")
    
    # Perform pathway enrichment analysis
    print_subheader("Protein Pathway Overrepresentation Analysis")
    
    enrichment_results = pathway_analyzer.perform_pathway_enrichment(all_de_proteins)
    
    if len(enrichment_results) > 0:
        print("Significantly enriched pathways:")
        significant_pathways = enrichment_results[enrichment_results['p_adj'] < 0.05]
        
        if len(significant_pathways) > 0:
            for _, pathway_row in significant_pathways.head(8).iterrows():
                pathway_name = pathway_row['pathway']
                p_value = pathway_row['p_value']
                p_adj = pathway_row['p_adj']
                overlap_size = pathway_row['overlap_size']
                pathway_size = pathway_row['pathway_size']
                enrichment_ratio = pathway_row['enrichment_ratio']
                overlap_proteins = pathway_row['overlap_proteins']
                
                print(f"\n  {pathway_name}:")
                print(f"    P-value: {p_value:.2e} (adjusted: {p_adj:.2e})")
                print(f"    Overlap: {overlap_size}/{pathway_size} proteins")
                print(f"    Enrichment ratio: {enrichment_ratio:.2f}")
                print(f"    Proteins: {', '.join(overlap_proteins)}")
        else:
            print("No pathways reached significance threshold (p_adj < 0.05)")
    else:
        print("No pathway enrichment results found")
    
    # Analyze protein functions
    print_subheader("Functional Analysis of Differentially Expressed Proteins")
    
    function_analysis = pathway_analyzer.analyze_protein_functions(all_de_proteins)
    
    if len(function_analysis) > 0:
        print("Protein functional categories:")
        for _, func_row in function_analysis.iterrows():
            protein = func_row['protein']
            n_pathways = func_row['n_pathways']
            pathways = func_row['pathways']
            
            if n_pathways > 0:
                print(f"  {protein}: {n_pathways} pathways - {pathways}")
            else:
                print(f"  {protein}: No known pathways")
    
    return enrichment_results if 'enrichment_results' in locals() else None


def demonstrate_visualization():
    """Demonstrate proteomic data visualization."""
    print_header("PROTEOMIC DATA VISUALIZATION")
    
    # Generate data for visualization
    abundance_matrix = generate_sample_abundance_matrix(n_proteins=3000, n_samples=12)
    metadata = get_sample_metadata_df()
    
    processor = ProteomeDataProcessor()
    processed_data = processor.handle_missing_values(abundance_matrix, method='knn')
    normalized_data = processor.normalize_abundance_data(processed_data, method='median_centering')
    
    # Initialize visualizer
    visualizer = ProteomicsVisualizer()
    
    print_subheader("Creating Visualizations")
    
    try:
        # PCA plot
        print("Creating PCA plot...")
        fig1 = visualizer.plot_pca(normalized_data, metadata=metadata, color_by='condition')
        plt.savefig('proteomics_pca.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("✓ Saved PCA plot")
        
        # Abundance heatmap
        print("Creating protein abundance heatmap...")
        fig2 = visualizer.plot_abundance_heatmap(normalized_data, sample_metadata=metadata, top_n=50)
        plt.savefig('proteomics_heatmap.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("✓ Saved abundance heatmap")
        
        # Generate mock DE results for volcano plot
        de_analyzer = DifferentialProteinAnalyzer()
        infection_samples = [s for s, info in SAMPLE_METADATA.items() 
                           if info['condition'] == 'Infection' and s in normalized_data.columns]
        healthy_samples = [s for s, info in SAMPLE_METADATA.items() 
                         if info['condition'] == 'Healthy' and s in normalized_data.columns]
        
        if len(infection_samples) >= 2 and len(healthy_samples) >= 2:
            de_results = de_analyzer.perform_ttest_analysis(
                normalized_data, infection_samples, healthy_samples
            )
            
            if len(de_results) > 0:
                # Volcano plot
                print("Creating volcano plot...")
                fig3 = visualizer.plot_volcano(de_results, p_col='p_adj', fc_col='log2_fold_change')
                plt.savefig('proteomics_volcano.png', dpi=300, bbox_inches='tight')
                plt.close()
                print("✓ Saved volcano plot")
        
        # Protein abundance comparison
        print("Creating protein abundance comparison...")
        demo_proteins = ['CRP', 'ALB', 'HP', 'APOA1']
        available_proteins = [p for p in demo_proteins if p in normalized_data.index]
        
        if available_proteins:
            fig4 = visualizer.plot_protein_abundance_comparison(
                normalized_data, available_proteins, metadata=metadata, group_by='condition'
            )
            plt.savefig('proteomics_protein_comparison.png', dpi=300, bbox_inches='tight')
            plt.close()
            print("✓ Saved protein abundance comparison")
        
        # Abundance distribution
        print("Creating abundance distribution plots...")
        fig5 = visualizer.plot_abundance_distribution(
            normalized_data, sample_metadata=metadata, group_by='condition'
        )
        plt.savefig('proteomics_abundance_distribution.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("✓ Saved abundance distribution plots")
        
        # Missing value visualization
        print("Creating missing value pattern visualization...")
        fig6 = visualizer.plot_missing_values(abundance_matrix)
        plt.savefig('proteomics_missing_values.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("✓ Saved missing value visualization")
        
        # Pathway enrichment plot (mock data)
        print("Creating pathway enrichment plot...")
        mock_enrichment = pd.DataFrame({
            'pathway': ['Acute Phase Response', 'Complement System', 'Blood Coagulation', 
                       'Immunoglobulins', 'Lipid Transport'],
            'p_value': [0.001, 0.003, 0.008, 0.02, 0.04],
            'enrichment_ratio': [3.2, 2.8, 2.1, 1.9, 1.6],
            'overlap_size': [8, 6, 5, 4, 3],
            'pathway_size': [15, 12, 18, 16, 14]
        })
        
        fig7 = visualizer.plot_pathway_enrichment(mock_enrichment)
        plt.savefig('proteomics_pathways.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("✓ Saved pathway enrichment plot")
        
    except Exception as e:
        print(f"Visualization error: {e}")
        print("Note: Some visualizations may require additional dependencies")


def demonstrate_clinical_insights():
    """Demonstrate clinical insights from proteomic analysis."""
    print_header("CLINICAL INSIGHTS FROM PROTEOMIC ANALYSIS")
    
    print_subheader("Acute Phase Response in Infection")
    print("""
Key insights from Patient Z's infection and recovery proteomics:

1. POSITIVE ACUTE PHASE PROTEINS (Increase during infection):
   - CRP (C-Reactive Protein): Dramatic increase (>10-fold)
   - SAA1 (Serum Amyloid A): Major inflammatory marker
   - Haptoglobin (HP): Binds free hemoglobin, prevents oxidative damage
   - Fibrinogen (FGA/FGB): Coagulation and inflammatory response

2. NEGATIVE ACUTE PHASE PROTEINS (Decrease during infection):
   - Albumin (ALB): Reduced synthesis, increased vascular permeability
   - Apolipoproteins (APOA1/APOB): Altered lipid metabolism
   - Transferrin: Iron sequestration as antimicrobial strategy

3. COMPLEMENT SYSTEM MODULATION:
   - C3, C4A: Central complement components
   - CFD: Initially decreased, recovers during healing
   - Regulatory proteins maintain immune balance
    """)
    
    print_subheader("Precision Medicine Applications")
    print("""
Proteomic analysis enables personalized medicine approaches:

1. BIOMARKER DISCOVERY:
   - Protein signatures predict infection severity
   - Recovery markers guide treatment decisions
   - Individual baseline profiles enable personalized monitoring

2. THERAPEUTIC MONITORING:
   - Real-time assessment of treatment response
   - Early detection of complications
   - Optimization of intervention timing

3. DRUG DEVELOPMENT:
   - Target identification through protein pathway analysis
   - Mechanism of action studies
   - Biomarker-driven clinical trials

4. CLINICAL DECISION SUPPORT:
   - Risk stratification based on protein profiles
   - Personalized treatment recommendations
   - Prediction of treatment outcomes
    """)
    
    print_subheader("Systems Biology Insights")
    print("""
Pathway-level understanding from proteomic analysis:

1. INFLAMMATORY NETWORKS:
   - Acute phase response coordination
   - Complement cascade activation
   - Cytokine signaling pathways

2. METABOLIC REPROGRAMMING:
   - Energy metabolism shifts during infection
   - Lipid transport alterations
   - Protein synthesis redirection

3. HOMEOSTATIC MECHANISMS:
   - Coagulation system activation
   - Antioxidant response upregulation
   - Tissue repair protein expression

4. RECOVERY PROCESSES:
   - Gradual normalization of protein levels
   - Resolution of inflammatory responses
   - Restoration of metabolic balance
    """)
    
    print_subheader("Future Directions")
    print("""
Emerging applications of proteomics in precision medicine:

1. MULTI-OMICS INTEGRATION:
   - Combining proteomics with genomics and metabolomics
   - Systems-level understanding of disease mechanisms
   - Comprehensive biomarker panels

2. REAL-TIME MONITORING:
   - Point-of-care protein analysis
   - Continuous health monitoring
   - Early disease detection

3. THERAPEUTIC DEVELOPMENT:
   - Protein-based drug targets
   - Companion diagnostics
   - Personalized treatment strategies

4. POPULATION HEALTH:
   - Protein-based risk assessment
   - Public health surveillance
   - Preventive medicine applications
    """)


def main():
    """Main demonstration function."""
    print("="*60)
    print(" PROTEOMIC DATA ANALYSIS DEMONSTRATION")
    print(" Stanford Data Ocean - Data Science in Precision Medicine")
    print("="*60)
    
    print("""
This demonstration showcases the capabilities of the proteomic analysis module
using data inspired by the Stanford iPOP study, including:

• Mass spectrometry-based protein quantification analysis
• Differential protein expression analysis
• Pathway enrichment and functional interpretation
• Clinical insights and precision medicine applications
• Comprehensive visualization tools

The analysis focuses on Patient Z's infection and recovery, demonstrating
how proteomics enables personalized medicine through detailed protein profiling.
    """)
    
    try:
        # Run demonstrations
        demonstrate_ipop_study_context()
        abundance_data, processed_data, normalized_data, metadata = demonstrate_data_loading()
        de_results = demonstrate_differential_protein_analysis()
        enrichment_results = demonstrate_pathway_analysis()
        demonstrate_visualization()
        demonstrate_clinical_insights()
        
        print_header("DEMONSTRATION COMPLETE")
        print("""
Summary of generated files:
• proteomics_pca.png - Principal component analysis of samples
• proteomics_heatmap.png - Protein abundance heatmap
• proteomics_volcano.png - Volcano plot of differential expression
• proteomics_protein_comparison.png - Key protein abundance comparison
• proteomics_abundance_distribution.png - Overall abundance distributions
• proteomics_missing_values.png - Missing value pattern analysis
• proteomics_pathways.png - Pathway enrichment visualization

Key findings demonstrated:
• Quality control and normalization of MS-based protein data
• Identification of differentially expressed proteins during infection
• Pathway-level analysis of acute phase response
• Clinical interpretation of proteomic changes
• Personalized medicine applications

Clinical insights:
• Acute phase response proteins as infection biomarkers
• Complement system modulation during immune response
• Metabolic reprogramming reflected in protein changes
• Recovery monitoring through protein normalization

Next steps:
1. Explore the Jupyter notebook for interactive analysis
2. Apply methods to your own proteomics datasets
3. Integrate with clinical data for biomarker discovery
4. Develop personalized protein signatures

For more information, see the README.md file and comprehensive documentation.
        """)
        
    except Exception as e:
        print(f"\nError during demonstration: {e}")
        print("Please check your Python environment and dependencies.")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main()) 