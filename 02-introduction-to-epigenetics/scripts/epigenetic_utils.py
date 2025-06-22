"""
Epigenetic Data Analysis Utilities

This module provides comprehensive tools for analyzing epigenomic data,
including DNA methylation patterns and histone modifications.

Classes:
    - EpigenomeDataProcessor: Core data processing and analysis
    - MethylationAnalyzer: Specialized DNA methylation analysis
    - HistoneModificationAnalyzer: ChIP-seq and histone modification analysis
    - EpigenomeVisualizer: Visualization tools for epigenomic data

Author: Stanford Data Ocean - Data Science in Precision Medicine
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')


class EpigenomeDataProcessor:
    """
    Core class for processing and analyzing epigenomic datasets.
    
    This class handles loading, processing, and basic analysis of epigenomic
    data from the NIH Roadmap Epigenomics Program and similar sources.
    """
    
    def __init__(self):
        """Initialize the epigenome data processor."""
        self.methylation_data = None
        self.histone_data = None
        self.sample_metadata = None
        
    def load_methylation_data(self, filepath, sample_limit=None):
        """
        Load DNA methylation data from WGBS DMR files.
        
        Args:
            filepath (str): Path to the methylation data file
            sample_limit (int): Limit number of samples for memory management
            
        Returns:
            pd.DataFrame: Processed methylation data
        """
        try:
            data = pd.read_csv(filepath, sep='\t', low_memory=False)
            
            if sample_limit:
                # Keep genomic coordinates and limited samples
                coord_cols = ['chr', 'start', 'end']
                sample_cols = [col for col in data.columns if col not in coord_cols][:sample_limit]
                data = data[coord_cols + sample_cols]
            
            # Create genomic region identifiers
            data['region_id'] = data['chr'].astype(str) + ':' + \
                               data['start'].astype(str) + '-' + \
                               data['end'].astype(str)
            
            self.methylation_data = data
            print(f"Loaded methylation data: {data.shape[0]} regions, {data.shape[1]-4} samples")
            return data
            
        except Exception as e:
            print(f"Error loading methylation data: {e}")
            return None
    
    def load_sample_metadata(self, filepath):
        """
        Load sample metadata for epigenome annotations.
        
        Args:
            filepath (str): Path to metadata file
            
        Returns:
            pd.DataFrame: Sample metadata
        """
        try:
            if filepath.endswith('.xlsx'):
                metadata = pd.read_excel(filepath)
            else:
                metadata = pd.read_csv(filepath)
            
            self.sample_metadata = metadata
            print(f"Loaded metadata for {len(metadata)} samples")
            return metadata
            
        except Exception as e:
            print(f"Error loading metadata: {e}")
            return None
    
    def filter_high_variance_regions(self, data, top_n=5000):
        """
        Filter for most variable methylation regions across samples.
        
        Args:
            data (pd.DataFrame): Methylation data
            top_n (int): Number of top variable regions to keep
            
        Returns:
            pd.DataFrame: Filtered data with most variable regions
        """
        # Get numeric columns (samples)
        sample_cols = [col for col in data.columns 
                      if col not in ['chr', 'start', 'end', 'region_id']]
        
        # Calculate variance for each region
        region_variance = data[sample_cols].var(axis=1, skipna=True)
        
        # Get top variable regions
        top_regions = region_variance.nlargest(top_n).index
        filtered_data = data.loc[top_regions].copy()
        
        print(f"Filtered to {len(filtered_data)} most variable regions")
        return filtered_data
    
    def calculate_methylation_statistics(self, data):
        """
        Calculate basic statistics for methylation data.
        
        Args:
            data (pd.DataFrame): Methylation data
            
        Returns:
            dict: Dictionary containing various statistics
        """
        sample_cols = [col for col in data.columns 
                      if col not in ['chr', 'start', 'end', 'region_id']]
        
        stats_dict = {
            'total_regions': len(data),
            'total_samples': len(sample_cols),
            'mean_methylation_per_sample': data[sample_cols].mean().to_dict(),
            'median_methylation_per_sample': data[sample_cols].median().to_dict(),
            'missing_data_per_sample': data[sample_cols].isnull().sum().to_dict(),
            'chromosomes_covered': sorted(data['chr'].unique()),
            'methylation_range': {
                'min': data[sample_cols].min().min(),
                'max': data[sample_cols].max().max()
            }
        }
        
        return stats_dict
    
    def identify_differential_regions(self, data, group1_samples, group2_samples, 
                                    min_diff=0.2, p_threshold=0.05):
        """
        Identify differentially methylated regions between two groups.
        
        Args:
            data (pd.DataFrame): Methylation data
            group1_samples (list): Sample names for group 1
            group2_samples (list): Sample names for group 2
            min_diff (float): Minimum methylation difference threshold
            p_threshold (float): P-value threshold for significance
            
        Returns:
            pd.DataFrame: Differentially methylated regions
        """
        # Filter for samples that exist in the data
        available_cols = data.columns
        group1_available = [s for s in group1_samples if s in available_cols]
        group2_available = [s for s in group2_samples if s in available_cols]
        
        if not group1_available or not group2_available:
            print("Warning: Some samples not found in data")
            return pd.DataFrame()
        
        # Calculate group means
        group1_mean = data[group1_available].mean(axis=1)
        group2_mean = data[group2_available].mean(axis=1)
        
        # Calculate difference
        methylation_diff = group1_mean - group2_mean
        
        # Perform t-tests
        p_values = []
        for idx in data.index:
            try:
                group1_vals = data.loc[idx, group1_available].dropna()
                group2_vals = data.loc[idx, group2_available].dropna()
                
                if len(group1_vals) >= 2 and len(group2_vals) >= 2:
                    _, p_val = stats.ttest_ind(group1_vals, group2_vals)
                    p_values.append(p_val)
                else:
                    p_values.append(1.0)
            except:
                p_values.append(1.0)
        
        # Create results dataframe
        results = data[['chr', 'start', 'end', 'region_id']].copy()
        results['group1_mean'] = group1_mean
        results['group2_mean'] = group2_mean
        results['methylation_diff'] = methylation_diff
        results['p_value'] = p_values
        
        # Filter for significant differences
        significant = results[
            (abs(results['methylation_diff']) >= min_diff) & 
            (results['p_value'] <= p_threshold)
        ].copy()
        
        print(f"Found {len(significant)} differentially methylated regions")
        return significant.sort_values('p_value')


class MethylationAnalyzer:
    """
    Specialized class for DNA methylation analysis.
    
    Provides advanced methods for analyzing DNA methylation patterns,
    including CpG island analysis and methylation clustering.
    """
    
    def __init__(self):
        """Initialize methylation analyzer."""
        pass
    
    def classify_methylation_levels(self, methylation_values):
        """
        Classify methylation levels into categories.
        
        Args:
            methylation_values (array-like): Methylation values (0-1)
            
        Returns:
            list: Classification labels
        """
        classifications = []
        for val in methylation_values:
            if pd.isna(val):
                classifications.append('Unknown')
            elif val < 0.2:
                classifications.append('Unmethylated')
            elif val < 0.8:
                classifications.append('Partially Methylated')
            else:
                classifications.append('Highly Methylated')
        
        return classifications
    
    def analyze_cpg_content(self, sequence):
        """
        Analyze CpG content in a DNA sequence.
        
        Args:
            sequence (str): DNA sequence
            
        Returns:
            dict: CpG analysis results
        """
        sequence = sequence.upper()
        
        # Count CpG dinucleotides
        cpg_count = sequence.count('CG')
        
        # Calculate GC content
        gc_count = sequence.count('G') + sequence.count('C')
        gc_content = gc_count / len(sequence) if len(sequence) > 0 else 0
        
        # Calculate CpG observed/expected ratio
        c_count = sequence.count('C')
        g_count = sequence.count('G')
        expected_cpg = (c_count * g_count) / len(sequence) if len(sequence) > 0 else 0
        obs_exp_ratio = cpg_count / expected_cpg if expected_cpg > 0 else 0
        
        return {
            'cpg_count': cpg_count,
            'gc_content': gc_content,
            'cpg_density': cpg_count / len(sequence) if len(sequence) > 0 else 0,
            'obs_exp_ratio': obs_exp_ratio,
            'is_cpg_island': obs_exp_ratio > 0.6 and gc_content > 0.5
        }
    
    def perform_methylation_clustering(self, data, method='ward'):
        """
        Perform hierarchical clustering on methylation data.
        
        Args:
            data (pd.DataFrame): Methylation data (samples as columns)
            method (str): Clustering method
            
        Returns:
            dict: Clustering results
        """
        # Get numeric data and handle missing values
        numeric_cols = data.select_dtypes(include=[np.number]).columns
        cluster_data = data[numeric_cols].fillna(data[numeric_cols].mean())
        
        # Perform clustering
        linkage_matrix = linkage(cluster_data.T, method=method)
        
        return {
            'linkage_matrix': linkage_matrix,
            'sample_order': cluster_data.columns,
            'method': method
        }


class HistoneModificationAnalyzer:
    """
    Class for analyzing histone modification ChIP-seq data.
    
    Handles processing and analysis of histone modification patterns
    around genes and regulatory elements.
    """
    
    def __init__(self):
        """Initialize histone modification analyzer."""
        self.histone_marks = {
            'H3K4me3': 'Active promoter mark',
            'H3K4me1': 'Enhancer mark',
            'H3K27ac': 'Active enhancer mark',
            'H3K27me3': 'Repressive mark',
            'H3K9me3': 'Heterochromatin mark'
        }
    
    def load_chip_seq_data(self, filepath):
        """
        Load ChIP-seq peak data.
        
        Args:
            filepath (str): Path to ChIP-seq data file
            
        Returns:
            pd.DataFrame: ChIP-seq data
        """
        try:
            data = pd.read_csv(filepath, sep='\t')
            print(f"Loaded ChIP-seq data: {data.shape}")
            return data
        except Exception as e:
            print(f"Error loading ChIP-seq data: {e}")
            return None
    
    def calculate_tss_enrichment(self, data, tss_positions, window_size=2000):
        """
        Calculate histone modification enrichment around transcription start sites.
        
        Args:
            data (pd.DataFrame): ChIP-seq data with genomic coordinates
            tss_positions (pd.DataFrame): TSS positions
            window_size (int): Window size around TSS
            
        Returns:
            pd.DataFrame: Enrichment scores around TSS
        """
        enrichment_profiles = []
        
        for _, tss in tss_positions.iterrows():
            # Define window around TSS
            start_pos = tss['tss_position'] - window_size // 2
            end_pos = tss['tss_position'] + window_size // 2
            
            # Find overlapping ChIP-seq peaks
            overlapping = data[
                (data['chr'] == tss['chr']) &
                (data['start'] <= end_pos) &
                (data['end'] >= start_pos)
            ]
            
            # Calculate enrichment score (sum of signal values)
            enrichment_score = overlapping['signal'].sum() if 'signal' in data.columns else len(overlapping)
            
            enrichment_profiles.append({
                'gene': tss['gene_name'],
                'chr': tss['chr'],
                'tss_position': tss['tss_position'],
                'enrichment_score': enrichment_score,
                'peak_count': len(overlapping)
            })
        
        return pd.DataFrame(enrichment_profiles)
    
    def classify_chromatin_states(self, histone_data, thresholds=None):
        """
        Classify chromatin states based on histone modification patterns.
        
        Args:
            histone_data (pd.DataFrame): Data with histone modification signals
            thresholds (dict): Thresholds for each mark
            
        Returns:
            pd.DataFrame: Chromatin state classifications
        """
        if thresholds is None:
            thresholds = {
                'H3K4me3': 2.0,
                'H3K4me1': 1.5,
                'H3K27ac': 1.5,
                'H3K27me3': 1.0,
                'H3K9me3': 1.0
            }
        
        states = []
        
        for _, row in histone_data.iterrows():
            # Check for active promoter (H3K4me3 high)
            if row.get('H3K4me3', 0) > thresholds['H3K4me3']:
                state = 'Active Promoter'
            
            # Check for active enhancer (H3K4me1 + H3K27ac)
            elif (row.get('H3K4me1', 0) > thresholds['H3K4me1'] and 
                  row.get('H3K27ac', 0) > thresholds['H3K27ac']):
                state = 'Active Enhancer'
            
            # Check for poised enhancer (H3K4me1 only)
            elif row.get('H3K4me1', 0) > thresholds['H3K4me1']:
                state = 'Poised Enhancer'
            
            # Check for repressed regions
            elif row.get('H3K27me3', 0) > thresholds['H3K27me3']:
                state = 'Polycomb Repressed'
            
            elif row.get('H3K9me3', 0) > thresholds['H3K9me3']:
                state = 'Heterochromatin'
            
            else:
                state = 'Low Signal'
            
            states.append(state)
        
        histone_data_copy = histone_data.copy()
        histone_data_copy['chromatin_state'] = states
        
        return histone_data_copy


class EpigenomeVisualizer:
    """
    Visualization tools for epigenomic data.
    
    Provides methods for creating publication-quality plots of
    methylation patterns and histone modifications.
    """
    
    def __init__(self):
        """Initialize visualizer with custom styling."""
        plt.style.use('seaborn-v0_8')
        self.colors = {
            'methylated': '#d62728',
            'unmethylated': '#2ca02c',
            'partial': '#ff7f0e',
            'active': '#1f77b4',
            'repressive': '#9467bd'
        }
    
    def plot_methylation_heatmap(self, data, sample_metadata=None, 
                               figsize=(12, 8), cluster=True):
        """
        Create a heatmap of methylation levels across samples and regions.
        
        Args:
            data (pd.DataFrame): Methylation data
            sample_metadata (pd.DataFrame): Sample annotations
            figsize (tuple): Figure size
            cluster (bool): Whether to cluster samples and regions
            
        Returns:
            matplotlib.figure.Figure: The created figure
        """
        # Get numeric columns (samples)
        sample_cols = [col for col in data.columns 
                      if col not in ['chr', 'start', 'end', 'region_id']]
        
        # Prepare data for heatmap
        plot_data = data[sample_cols].fillna(0)
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Create heatmap
        if cluster:
            sns.clustermap(plot_data, 
                         cmap='RdYlBu_r', 
                         center=0.5,
                         figsize=figsize,
                         cbar_kws={'label': 'Methylation Level'})
        else:
            sns.heatmap(plot_data, 
                       cmap='RdYlBu_r', 
                       center=0.5,
                       cbar_kws={'label': 'Methylation Level'},
                       ax=ax)
        
        plt.title('DNA Methylation Levels Across Samples')
        plt.xlabel('Samples')
        plt.ylabel('Genomic Regions')
        
        return fig
    
    def plot_methylation_distribution(self, data, sample_name=None, figsize=(10, 6)):
        """
        Plot distribution of methylation levels.
        
        Args:
            data (pd.DataFrame or pd.Series): Methylation data
            sample_name (str): Name for the sample/dataset
            figsize (tuple): Figure size
            
        Returns:
            matplotlib.figure.Figure: The created figure
        """
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
        
        # Get methylation values
        if isinstance(data, pd.DataFrame):
            sample_cols = [col for col in data.columns 
                          if col not in ['chr', 'start', 'end', 'region_id']]
            meth_values = data[sample_cols].values.flatten()
        else:
            meth_values = data.values
        
        # Remove NaN values
        meth_values = meth_values[~pd.isna(meth_values)]
        
        # Histogram
        ax1.hist(meth_values, bins=50, alpha=0.7, color=self.colors['active'])
        ax1.set_xlabel('Methylation Level')
        ax1.set_ylabel('Frequency')
        ax1.set_title('Methylation Level Distribution')
        
        # Box plot by methylation category
        categories = []
        values = []
        for val in meth_values:
            if val < 0.2:
                categories.append('Unmethylated')
                values.append(val)
            elif val < 0.8:
                categories.append('Partial')
                values.append(val)
            else:
                categories.append('Methylated')
                values.append(val)
        
        cat_df = pd.DataFrame({'Category': categories, 'Methylation': values})
        sns.boxplot(data=cat_df, x='Category', y='Methylation', ax=ax2)
        ax2.set_title('Methylation by Category')
        
        plt.tight_layout()
        return fig
    
    def plot_histone_modification_profile(self, data, gene_name=None, 
                                        window_size=5000, figsize=(12, 6)):
        """
        Plot histone modification profile around a gene or region.
        
        Args:
            data (pd.DataFrame): Histone modification data
            gene_name (str): Gene name for title
            window_size (int): Window size for the profile
            figsize (tuple): Figure size
            
        Returns:
            matplotlib.figure.Figure: The created figure
        """
        fig, ax = plt.subplots(figsize=figsize)
        
        # Create position array
        positions = np.linspace(-window_size//2, window_size//2, len(data))
        
        # Plot each histone mark
        histone_marks = ['H3K4me3', 'H3K4me1', 'H3K27ac', 'H3K27me3', 'H3K9me3']
        colors = ['red', 'orange', 'green', 'blue', 'purple']
        
        for mark, color in zip(histone_marks, colors):
            if mark in data.columns:
                ax.plot(positions, data[mark], label=mark, color=color, linewidth=2)
        
        ax.axvline(x=0, color='black', linestyle='--', alpha=0.5, label='TSS')
        ax.set_xlabel('Distance from TSS (bp)')
        ax.set_ylabel('ChIP-seq Signal')
        ax.set_title(f'Histone Modification Profile{" - " + gene_name if gene_name else ""}')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        return fig
    
    def plot_chromatin_state_pie(self, chromatin_states, figsize=(8, 8)):
        """
        Create a pie chart of chromatin state distribution.
        
        Args:
            chromatin_states (pd.Series or list): Chromatin state classifications
            figsize (tuple): Figure size
            
        Returns:
            matplotlib.figure.Figure: The created figure
        """
        fig, ax = plt.subplots(figsize=figsize)
        
        # Count chromatin states
        if isinstance(chromatin_states, pd.Series):
            state_counts = chromatin_states.value_counts()
        else:
            state_counts = pd.Series(chromatin_states).value_counts()
        
        # Define colors for different states
        state_colors = {
            'Active Promoter': '#d62728',
            'Active Enhancer': '#2ca02c',
            'Poised Enhancer': '#ff7f0e',
            'Polycomb Repressed': '#9467bd',
            'Heterochromatin': '#8c564b',
            'Low Signal': '#e377c2'
        }
        
        colors = [state_colors.get(state, '#1f77b4') for state in state_counts.index]
        
        # Create pie chart
        wedges, texts, autotexts = ax.pie(state_counts.values, 
                                         labels=state_counts.index,
                                         colors=colors,
                                         autopct='%1.1f%%',
                                         startangle=90)
        
        ax.set_title('Chromatin State Distribution')
        
        return fig
    
    def plot_pca_analysis(self, data, metadata=None, color_by=None, figsize=(10, 8)):
        """
        Perform and plot PCA analysis of epigenomic data.
        
        Args:
            data (pd.DataFrame): Epigenomic data (samples as columns)
            metadata (pd.DataFrame): Sample metadata
            color_by (str): Column in metadata to color samples by
            figsize (tuple): Figure size
            
        Returns:
            matplotlib.figure.Figure: The created figure
        """
        # Prepare data
        sample_cols = [col for col in data.columns 
                      if col not in ['chr', 'start', 'end', 'region_id']]
        
        # Transpose and standardize
        pca_data = data[sample_cols].T.fillna(0)
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(pca_data)
        
        # Perform PCA
        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(scaled_data)
        
        # Create plot
        fig, ax = plt.subplots(figsize=figsize)
        
        if metadata is not None and color_by is not None:
            # Color by metadata category
            for category in metadata[color_by].unique():
                mask = metadata[color_by] == category
                ax.scatter(pca_result[mask, 0], pca_result[mask, 1], 
                          label=category, alpha=0.7, s=50)
        else:
            ax.scatter(pca_result[:, 0], pca_result[:, 1], alpha=0.7, s=50)
        
        ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
        ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
        ax.set_title('PCA of Epigenomic Data')
        
        if metadata is not None and color_by is not None:
            ax.legend()
        
        ax.grid(True, alpha=0.3)
        
        return fig 