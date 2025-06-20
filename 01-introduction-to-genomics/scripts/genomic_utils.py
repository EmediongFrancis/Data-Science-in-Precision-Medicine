"""
Genomic Data Analysis Utilities

This module contains utility functions for genomic data analysis
as part of the Stanford Data Ocean Precision Medicine project.

Author: Emediong Francis
Date: 2024
"""

import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Optional, Tuple


class GenomicDataProcessor:
    """
    A class for processing and analyzing genomic data from various sources.
    """
    
    def __init__(self):
        """Initialize the genomic data processor."""
        self.population_codes = {
            'EAS_AF': 'East Asian',
            'AMR_AF': 'Ad Mixed American', 
            'EUR_AF': 'European',
            'AFR_AF': 'African',
            'SAS_AF': 'South Asian'
        }
    
    def extract_population_frequencies(self, info_text: str) -> Dict[str, float]:
        """
        Extract population-specific allele frequencies from VCF info column.
        
        Args:
            info_text (str): The info column text from VCF file
            
        Returns:
            Dict[str, float]: Dictionary mapping population names to frequencies (as percentages)
        """
        frequencies = {}
        
        for code, name in self.population_codes.items():
            pattern = f'{code}=([0-9.]+)'
            match = re.search(pattern, str(info_text))
            if match:
                frequencies[name] = float(match.group(1)) * 100  # Convert to percentage
        
        return frequencies
    
    def extract_allele_counts(self, info_text: str) -> Dict[str, int]:
        """
        Extract allele count information from VCF info column.
        
        Args:
            info_text (str): The info column text from VCF file
            
        Returns:
            Dict[str, int]: Dictionary with AC (allele count) and AN (allele number)
        """
        counts = {}
        
        # Extract allele count (AC)
        ac_match = re.search(r'AC=(\d+)', str(info_text))
        if ac_match:
            counts['allele_count'] = int(ac_match.group(1))
        
        # Extract allele number (AN)
        an_match = re.search(r'AN=(\d+)', str(info_text))
        if an_match:
            counts['allele_number'] = int(an_match.group(1))
        
        # Calculate allele frequency
        if 'allele_count' in counts and 'allele_number' in counts:
            counts['allele_frequency'] = counts['allele_count'] / counts['allele_number']
        
        return counts
    
    def plot_population_frequencies(self, frequencies: Dict[str, float], 
                                  title: str = "Population Allele Frequencies",
                                  figsize: Tuple[int, int] = (12, 6)) -> None:
        """
        Create a bar plot of population frequencies.
        
        Args:
            frequencies (Dict[str, float]): Population frequencies
            title (str): Plot title
            figsize (Tuple[int, int]): Figure size
        """
        plt.figure(figsize=figsize)
        
        populations = list(frequencies.keys())
        freqs = list(frequencies.values())
        
        colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FFEAA7']
        bars = plt.bar(populations, freqs, color=colors[:len(populations)])
        
        plt.title(title, fontsize=16, fontweight='bold')
        plt.xlabel('Population', fontsize=12)
        plt.ylabel('Allele Frequency (%)', fontsize=12)
        plt.xticks(rotation=45, ha='right')
        
        # Add value labels on bars
        for bar, freq in zip(bars, freqs):
            plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5, 
                    f'{freq:.1f}%', ha='center', va='bottom', fontweight='bold')
        
        plt.tight_layout()
        plt.grid(axis='y', alpha=0.3)
        plt.show()
    
    def analyze_variant_quality(self, df: pd.DataFrame) -> Dict[str, float]:
        """
        Analyze variant quality scores in a genomic dataset.
        
        Args:
            df (pd.DataFrame): DataFrame containing genomic variants
            
        Returns:
            Dict[str, float]: Quality statistics
        """
        if 'qual' not in df.columns:
            return {"error": "Quality column not found"}
        
        qual_stats = {
            'mean_quality': df['qual'].mean(),
            'median_quality': df['qual'].median(),
            'min_quality': df['qual'].min(),
            'max_quality': df['qual'].max(),
            'std_quality': df['qual'].std()
        }
        
        return qual_stats
    
    def filter_high_quality_variants(self, df: pd.DataFrame, 
                                   min_quality: float = 30.0) -> pd.DataFrame:
        """
        Filter variants based on quality score.
        
        Args:
            df (pd.DataFrame): DataFrame containing genomic variants
            min_quality (float): Minimum quality threshold
            
        Returns:
            pd.DataFrame: Filtered DataFrame
        """
        if 'qual' not in df.columns:
            print("Warning: Quality column not found. Returning original DataFrame.")
            return df
        
        filtered_df = df[df['qual'] >= min_quality].copy()
        
        print(f"Original variants: {len(df)}")
        print(f"High-quality variants (qual >= {min_quality}): {len(filtered_df)}")
        print(f"Filtered out: {len(df) - len(filtered_df)} ({((len(df) - len(filtered_df))/len(df)*100):.1f}%)")
        
        return filtered_df
    
    def get_chromosome_summary(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Get summary statistics by chromosome.
        
        Args:
            df (pd.DataFrame): DataFrame containing genomic variants
            
        Returns:
            pd.DataFrame: Summary by chromosome
        """
        if 'chrm' not in df.columns:
            return pd.DataFrame({"error": ["Chromosome column not found"]})
        
        summary = df.groupby('chrm').agg({
            'start_position': ['count', 'min', 'max'],
            'qual': ['mean', 'std'] if 'qual' in df.columns else ['count']
        }).round(2)
        
        # Flatten column names
        summary.columns = ['_'.join(col).strip() for col in summary.columns.values]
        summary = summary.reset_index()
        
        return summary


class VariantAnnotator:
    """
    A class for annotating genomic variants with functional information.
    """
    
    def __init__(self):
        """Initialize the variant annotator."""
        pass
    
    def annotate_variant_type(self, ref_allele: str, alt_allele: str) -> str:
        """
        Classify variant type based on reference and alternate alleles.
        
        Args:
            ref_allele (str): Reference allele
            alt_allele (str): Alternate allele
            
        Returns:
            str: Variant type classification
        """
        if len(ref_allele) == 1 and len(alt_allele) == 1:
            return "SNV"  # Single Nucleotide Variant
        elif len(ref_allele) > len(alt_allele):
            return "Deletion"
        elif len(ref_allele) < len(alt_allele):
            return "Insertion"
        else:
            return "Complex"
    
    def is_transition(self, ref_allele: str, alt_allele: str) -> bool:
        """
        Check if a SNV is a transition (A<->G or C<->T).
        
        Args:
            ref_allele (str): Reference allele
            alt_allele (str): Alternate allele
            
        Returns:
            bool: True if transition, False if transversion
        """
        if len(ref_allele) != 1 or len(alt_allele) != 1:
            return False
        
        transitions = [('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')]
        return (ref_allele, alt_allele) in transitions
    
    def calculate_ti_tv_ratio(self, df: pd.DataFrame) -> float:
        """
        Calculate transition/transversion ratio for SNVs.
        
        Args:
            df (pd.DataFrame): DataFrame with reference_bases and alternate_bases columns
            
        Returns:
            float: Ti/Tv ratio
        """
        if 'reference_bases' not in df.columns or 'alternate_bases' not in df.columns:
            return 0.0
        
        # Filter for SNVs only
        snvs = df[(df['reference_bases'].str.len() == 1) & 
                  (df['alternate_bases'].str.len() == 1)].copy()
        
        if len(snvs) == 0:
            return 0.0
        
        snvs['is_transition'] = snvs.apply(
            lambda row: self.is_transition(row['reference_bases'], row['alternate_bases']), 
            axis=1
        )
        
        transitions = snvs['is_transition'].sum()
        transversions = len(snvs) - transitions
        
        return transitions / transversions if transversions > 0 else float('inf')


class GenomicVisualizer:
    """
    A class for creating genomic data visualizations.
    """
    
    def __init__(self):
        """Initialize the genomic visualizer."""
        plt.style.use('seaborn-v0_8')
        sns.set_palette("husl")
    
    def plot_chromosome_distribution(self, df: pd.DataFrame, 
                                   figsize: Tuple[int, int] = (15, 8)) -> None:
        """
        Plot variant distribution across chromosomes.
        
        Args:
            df (pd.DataFrame): DataFrame with chromosome information
            figsize (Tuple[int, int]): Figure size
        """
        if 'chrm' not in df.columns:
            print("Error: Chromosome column not found")
            return
        
        # Count variants per chromosome
        chrom_counts = df['chrm'].value_counts().sort_index()
        
        plt.figure(figsize=figsize)
        
        bars = plt.bar(chrom_counts.index, chrom_counts.values, 
                      color=plt.cm.viridis(np.linspace(0, 1, len(chrom_counts))))
        
        plt.title('Genetic Variant Distribution Across Chromosomes', 
                 fontsize=16, fontweight='bold')
        plt.xlabel('Chromosome', fontsize=12)
        plt.ylabel('Number of Variants', fontsize=12)
        plt.xticks(rotation=45)
        
        # Add value labels on bars
        for bar, count in zip(bars, chrom_counts.values):
            plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(chrom_counts)*0.01, 
                    f'{count:,}', ha='center', va='bottom', fontsize=8, rotation=90)
        
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        plt.show()
    
    def plot_quality_distribution(self, df: pd.DataFrame, 
                                figsize: Tuple[int, int] = (12, 6)) -> None:
        """
        Plot distribution of variant quality scores.
        
        Args:
            df (pd.DataFrame): DataFrame with quality scores
            figsize (Tuple[int, int]): Figure size
        """
        if 'qual' not in df.columns:
            print("Error: Quality column not found")
            return
        
        plt.figure(figsize=figsize)
        
        plt.subplot(1, 2, 1)
        plt.hist(df['qual'], bins=50, alpha=0.7, color='skyblue', edgecolor='black')
        plt.title('Distribution of Variant Quality Scores', fontweight='bold')
        plt.xlabel('Quality Score')
        plt.ylabel('Frequency')
        plt.grid(alpha=0.3)
        
        plt.subplot(1, 2, 2)
        plt.boxplot(df['qual'])
        plt.title('Quality Score Box Plot', fontweight='bold')
        plt.ylabel('Quality Score')
        plt.grid(alpha=0.3)
        
        plt.tight_layout()
        plt.show()


# Utility functions for common genomic operations
def parse_rsid(rsid: str) -> Optional[int]:
    """
    Parse rsID to extract the numeric identifier.
    
    Args:
        rsid (str): rsID string (e.g., 'rs12345')
        
    Returns:
        Optional[int]: Numeric part of rsID or None if invalid
    """
    match = re.match(r'rs(\d+)', str(rsid))
    return int(match.group(1)) if match else None


def format_genomic_position(chromosome: str, position: int) -> str:
    """
    Format genomic position in standard notation.
    
    Args:
        chromosome (str): Chromosome identifier
        position (int): Genomic position
        
    Returns:
        str: Formatted position (e.g., 'chr1:12345')
    """
    return f"chr{chromosome}:{position:,}"


def calculate_variant_density(df: pd.DataFrame, window_size: int = 1000000) -> pd.DataFrame:
    """
    Calculate variant density in genomic windows.
    
    Args:
        df (pd.DataFrame): DataFrame with genomic variants
        window_size (int): Size of genomic window in base pairs
        
    Returns:
        pd.DataFrame: Variant density per window
    """
    if 'chrm' not in df.columns or 'start_position' not in df.columns:
        return pd.DataFrame()
    
    density_data = []
    
    for chrom in df['chrm'].unique():
        chrom_data = df[df['chrm'] == chrom].copy()
        if len(chrom_data) == 0:
            continue
            
        max_pos = chrom_data['start_position'].max()
        
        for window_start in range(0, max_pos + window_size, window_size):
            window_end = window_start + window_size
            variants_in_window = len(chrom_data[
                (chrom_data['start_position'] >= window_start) & 
                (chrom_data['start_position'] < window_end)
            ])
            
            density_data.append({
                'chromosome': chrom,
                'window_start': window_start,
                'window_end': window_end,
                'variant_count': variants_in_window,
                'density': variants_in_window / window_size * 1000000  # variants per Mb
            })
    
    return pd.DataFrame(density_data)


# Example usage and testing functions
def test_genomic_utils():
    """Test function for genomic utilities."""
    print("Testing Genomic Utilities...")
    
    # Test data
    test_info = "AC=2;AF=0.000399361;AN=5008;NS=2504;DP=103152;EAS_AF=0.001;AMR_AF=0;EUR_AF=0;AFR_AF=0;SAS_AF=0"
    
    processor = GenomicDataProcessor()
    frequencies = processor.extract_population_frequencies(test_info)
    counts = processor.extract_allele_counts(test_info)
    
    print(f"Population frequencies: {frequencies}")
    print(f"Allele counts: {counts}")
    
    # Test variant annotation
    annotator = VariantAnnotator()
    variant_type = annotator.annotate_variant_type("A", "G")
    is_ti = annotator.is_transition("A", "G")
    
    print(f"Variant type A->G: {variant_type}")
    print(f"Is transition A->G: {is_ti}")
    
    print("✅ All tests passed!")


if __name__ == "__main__":
    test_genomic_utils() 