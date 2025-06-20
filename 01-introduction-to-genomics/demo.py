#!/usr/bin/env python3
"""
Genomic Data Analysis Demonstration
===================================

A simple demonstration script that showcases the key functionality
of the genomic data analysis module.

Author: Emediong Francis
Part of: Stanford Data Ocean - Data Science in Precision Medicine

Usage:
    python demo.py
"""

import sys
import os
sys.path.append('scripts')

# Import our custom utilities
from genomic_utils import GenomicDataProcessor, VariantAnnotator, GenomicVisualizer
from config import DATABASE_TABLES, DEMO_SNPS, get_demo_snp_info, get_aws_config

import pandas as pd
import matplotlib.pyplot as plt

def demonstrate_population_analysis():
    """Demonstrate population genetics analysis."""
    print("=" * 60)
    print("🧬 GENOMIC DATA ANALYSIS DEMONSTRATION")
    print("=" * 60)
    
    # Initialize processor
    processor = GenomicDataProcessor()
    
    # Demonstrate with sample data (simulated info field)
    print("\n👁️ EYE COLOR GENETICS ANALYSIS")
    print("-" * 40)
    
    # Sample VCF info field for rs12913832 (eye color SNP)
    sample_info = "AC=2;AN=5008;AF=0.000399361;EAS_AF=0.002;AMR_AF=0.2017;AFR_AF=0.0023;EUR_AF=0.6362;SAS_AF=0.028"
    
    snp_info = get_demo_snp_info("rs12913832")
    print(f"SNP: rs12913832 - {snp_info['name']}")
    print(f"Description: {snp_info['description']}")
    
    # Extract population frequencies
    pop_frequencies = processor.extract_population_frequencies(sample_info)
    
    print("\n🌍 Population Frequencies:")
    for population, frequency in pop_frequencies.items():
        print(f"   {population:20s}: {frequency:6.2f}%")
    
    # Create visualization
    print("\n📊 Creating population frequency visualization...")
    processor.plot_population_frequencies(
        pop_frequencies,
        title="Eye Color SNP (rs12913832) - Population Frequencies"
    )
    
    return pop_frequencies

def demonstrate_variant_annotation():
    """Demonstrate variant annotation functionality."""
    print("\n🔬 VARIANT ANNOTATION ANALYSIS")
    print("-" * 40)
    
    # Initialize annotator
    annotator = VariantAnnotator()
    
    # Sample variants for demonstration
    sample_variants = [
        ("A", "G"),  # Transition
        ("C", "T"),  # Transition  
        ("A", "T"),  # Transversion
        ("G", "C"),  # Transversion
        ("ATG", "A"),  # Deletion
        ("A", "ATG"),  # Insertion
    ]
    
    print("Variant Type Analysis:")
    transitions = 0
    transversions = 0
    
    for ref, alt in sample_variants:
        variant_type = annotator.annotate_variant_type(ref, alt)
        is_transition = annotator.is_transition(ref, alt) if len(ref) == 1 and len(alt) == 1 else False
        
        print(f"   {ref:>5s} → {alt:<5s}: {variant_type:12s} {'(Transition)' if is_transition else '(Transversion)' if variant_type == 'SNV' else ''}")
        
        if variant_type == "SNV":
            if is_transition:
                transitions += 1
            else:
                transversions += 1
    
    # Calculate Ti/Tv ratio
    ti_tv_ratio = transitions / transversions if transversions > 0 else float('inf')
    print(f"\n📈 Transition/Transversion Ratio: {ti_tv_ratio:.2f}")
    print(f"   Transitions: {transitions}")
    print(f"   Transversions: {transversions}")

def demonstrate_quality_analysis():
    """Demonstrate data quality analysis."""
    print("\n📊 DATA QUALITY ANALYSIS")
    print("-" * 40)
    
    # Create sample data for demonstration
    sample_data = pd.DataFrame({
        'chrm': ['1', '1', '2', '2', '3', '3', '4', '4', '5', '5'] * 5,
        'start_position': range(1000000, 1000050),
        'qual': [45.2, 23.1, 67.8, 12.3, 89.1, 34.5, 56.7, 78.9, 41.2, 29.8] * 5,
        'reference_bases': ['A', 'C', 'G', 'T'] * 12 + ['A', 'C'],
        'alternate_bases': ['G', 'T', 'A', 'C'] * 12 + ['T', 'G']
    })
    
    processor = GenomicDataProcessor()
    
    # Analyze quality
    quality_stats = processor.analyze_variant_quality(sample_data)
    
    print("Quality Statistics:")
    for stat_name, value in quality_stats.items():
        print(f"   {stat_name.replace('_', ' ').title():15s}: {value:.2f}")
    
    # Filter high-quality variants
    print(f"\nFiltering variants with quality >= 30...")
    high_qual_variants = processor.filter_high_quality_variants(sample_data, min_quality=30.0)
    
    # Get chromosome summary
    chr_summary = processor.get_chromosome_summary(sample_data)
    print(f"\nChromosome Summary:")
    print(chr_summary.to_string(index=False))

def demonstrate_visualization():
    """Demonstrate visualization capabilities."""
    print("\n📈 VISUALIZATION DEMONSTRATION")
    print("-" * 40)
    
    # Create sample chromosome data
    chr_data = pd.DataFrame({
        'chrm': ['1', '2', '3', '4', '5'],
        'variant_count': [1250, 980, 1100, 890, 1050]
    })
    
    visualizer = GenomicVisualizer()
    
    print("Creating chromosome distribution plot...")
    visualizer.plot_chromosome_distribution(chr_data)
    
    # Create sample quality data
    quality_data = pd.DataFrame({
        'qual': [45.2, 23.1, 67.8, 12.3, 89.1, 34.5, 56.7, 78.9, 41.2, 29.8] * 10
    })
    
    print("Creating quality distribution plot...")
    visualizer.plot_quality_distribution(quality_data)

def show_configuration():
    """Display configuration information."""
    print("\n⚙️ CONFIGURATION INFORMATION")
    print("-" * 40)
    
    print("Available Databases:")
    for key, table in DATABASE_TABLES.items():
        print(f"   {key:15s}: {table}")
    
    print(f"\nDemo SNPs Available:")
    for rsid, info in list(DEMO_SNPS.items())[:3]:  # Show first 3
        print(f"   {rsid}: {info['name']}")
        print(f"      → {info['description']}")
    
    print(f"\nAWS Configuration:")
    aws_config = get_aws_config()
    for key, value in aws_config.items():
        display_value = value if value else "Not configured"
        print(f"   {key:15s}: {display_value}")

def main():
    """Main demonstration function."""
    try:
        # Show configuration
        show_configuration()
        
        # Run demonstrations
        pop_frequencies = demonstrate_population_analysis()
        demonstrate_variant_annotation()
        demonstrate_quality_analysis()
        demonstrate_visualization()
        
        print("\n" + "=" * 60)
        print("✅ DEMONSTRATION COMPLETED SUCCESSFULLY!")
        print("=" * 60)
        print("\n🎯 Key Features Demonstrated:")
        print("   • Population genetics analysis")
        print("   • Variant type annotation")
        print("   • Data quality assessment")
        print("   • Genomic data visualization")
        print("\n📚 This demonstrates the capabilities of the genomic")
        print("   analysis module for precision medicine research.")
        
        # Keep plots open
        print("\n💡 Close the plot windows to exit the program.")
        plt.show()
        
    except Exception as e:
        print(f"\n❌ Error during demonstration: {e}")
        print("💡 Make sure all dependencies are installed:")
        print("   pip install -r requirements.txt")

if __name__ == "__main__":
    main() 