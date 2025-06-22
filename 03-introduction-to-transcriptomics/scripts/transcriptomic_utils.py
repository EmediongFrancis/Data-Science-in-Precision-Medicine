"""
Transcriptomic Data Analysis Utilities

This module provides comprehensive tools for analyzing RNA-seq data,
including bulk and single-cell transcriptomics analysis.

Classes:
    - TranscriptomeDataProcessor: Core data processing and analysis
    - DifferentialExpressionAnalyzer: Statistical analysis of gene expression
    - PathwayEnrichmentAnalyzer: Functional analysis and pathway enrichment
    - TranscriptomeVisualizer: Visualization tools for transcriptomic data

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
from sklearn.manifold import TSNE
import warnings
warnings.filterwarnings('ignore')


class TranscriptomeDataProcessor:
    """
    Core class for processing and analyzing transcriptomic datasets.
    
    This class handles loading, processing, and basic analysis of RNA-seq
    data including count matrices, sample metadata, and quality control.
    """
    
    def __init__(self):
        """Initialize the transcriptome data processor."""
        self.count_data = None
        self.sample_metadata = None
        self.gene_annotations = None
        
    def load_count_matrix(self, filepath, gene_id_col=0, sample_start_col=1):
        """
        Load RNA-seq count matrix from file.
        
        Args:
            filepath (str): Path to the count matrix file
            gene_id_col (int): Column index for gene IDs
            sample_start_col (int): Starting column index for sample data
            
        Returns:
            pd.DataFrame: Processed count matrix
        """
        try:
            data = pd.read_csv(filepath, sep='\t', index_col=gene_id_col)
            
            # Remove non-sample columns if present
            if 'transcript_id.s.' in data.columns:
                data = data.drop('transcript_id.s.', axis=1)
            
            # Round to integers (counts should be whole numbers)
            data = data.round(0).astype(int)
            
            # Remove genes with no expression across all samples
            data = data[data.sum(axis=1) > 0]
            
            self.count_data = data
            print(f"Loaded count matrix: {data.shape[0]} genes, {data.shape[1]} samples")
            return data
            
        except Exception as e:
            print(f"Error loading count matrix: {e}")
            return None
    
    def load_sample_metadata(self, filepath):
        """
        Load sample metadata for experimental design.
        
        Args:
            filepath (str): Path to metadata file
            
        Returns:
            pd.DataFrame: Sample metadata
        """
        try:
            metadata = pd.read_csv(filepath, sep='\t', index_col=0)
            
            # Convert categorical columns to factors
            for col in metadata.columns:
                if metadata[col].dtype == 'object':
                    metadata[col] = metadata[col].astype('category')
            
            self.sample_metadata = metadata
            print(f"Loaded metadata for {len(metadata)} samples")
            return metadata
            
        except Exception as e:
            print(f"Error loading sample metadata: {e}")
            return None
    
    def calculate_library_sizes(self, data=None):
        """
        Calculate library sizes (total counts per sample).
        
        Args:
            data (pd.DataFrame): Count matrix (uses self.count_data if None)
            
        Returns:
            pd.Series: Library sizes for each sample
        """
        if data is None:
            data = self.count_data
        
        if data is None:
            raise ValueError("No count data available")
        
        library_sizes = data.sum(axis=0)
        return library_sizes
    
    def filter_low_expression_genes(self, data=None, min_count=10, min_samples=3):
        """
        Filter out genes with low expression across samples.
        
        Args:
            data (pd.DataFrame): Count matrix
            min_count (int): Minimum count threshold
            min_samples (int): Minimum number of samples above threshold
            
        Returns:
            pd.DataFrame: Filtered count matrix
        """
        if data is None:
            data = self.count_data
        
        if data is None:
            raise ValueError("No count data available")
        
        # Keep genes with at least min_count in at least min_samples
        keep_genes = (data >= min_count).sum(axis=1) >= min_samples
        filtered_data = data[keep_genes]
        
        print(f"Filtered genes: {data.shape[0]} -> {filtered_data.shape[0]} "
              f"({filtered_data.shape[0]/data.shape[0]*100:.1f}% retained)")
        
        return filtered_data
    
    def normalize_counts(self, data=None, method='cpm'):
        """
        Normalize count data using various methods.
        
        Args:
            data (pd.DataFrame): Count matrix
            method (str): Normalization method ('cpm', 'tpm', 'log_cpm')
            
        Returns:
            pd.DataFrame: Normalized expression data
        """
        if data is None:
            data = self.count_data
        
        if data is None:
            raise ValueError("No count data available")
        
        if method == 'cpm':
            # Counts per million
            library_sizes = data.sum(axis=0)
            normalized = data.div(library_sizes, axis=1) * 1e6
            
        elif method == 'log_cpm':
            # Log2 counts per million
            library_sizes = data.sum(axis=0)
            cpm = data.div(library_sizes, axis=1) * 1e6
            normalized = np.log2(cpm + 1)
            
        elif method == 'tpm':
            # Transcripts per million (simplified - assumes uniform gene length)
            library_sizes = data.sum(axis=0)
            normalized = data.div(library_sizes, axis=1) * 1e6
            
        else:
            raise ValueError(f"Unknown normalization method: {method}")
        
        print(f"Normalized data using {method} method")
        return normalized
    
    def calculate_quality_metrics(self, data=None, metadata=None):
        """
        Calculate quality control metrics for samples.
        
        Args:
            data (pd.DataFrame): Count matrix
            metadata (pd.DataFrame): Sample metadata
            
        Returns:
            pd.DataFrame: Quality metrics
        """
        if data is None:
            data = self.count_data
        if metadata is None:
            metadata = self.sample_metadata
        
        if data is None:
            raise ValueError("No count data available")
        
        metrics = pd.DataFrame(index=data.columns)
        
        # Library size
        metrics['library_size'] = data.sum(axis=0)
        
        # Number of detected genes
        metrics['detected_genes'] = (data > 0).sum(axis=0)
        
        # Percentage of top 100 genes
        top100_counts = data.apply(lambda x: x.nlargest(100).sum(), axis=0)
        metrics['pct_top100'] = (top100_counts / metrics['library_size']) * 100
        
        # Add metadata if available
        if metadata is not None:
            metrics = metrics.join(metadata, how='left')
        
        return metrics


class DifferentialExpressionAnalyzer:
    """
    Class for statistical analysis of differential gene expression.
    
    Provides methods for identifying differentially expressed genes
    between conditions using various statistical approaches.
    """
    
    def __init__(self):
        """Initialize differential expression analyzer."""
        pass
    
    def perform_ttest(self, data, group1_samples, group2_samples, 
                     equal_var=False, min_expression=1):
        """
        Perform t-test for differential expression.
        
        Args:
            data (pd.DataFrame): Expression data (genes x samples)
            group1_samples (list): Sample names for group 1
            group2_samples (list): Sample names for group 2
            equal_var (bool): Assume equal variances
            min_expression (float): Minimum expression threshold
            
        Returns:
            pd.DataFrame: Results with statistics and p-values
        """
        results = []
        
        for gene in data.index:
            group1_vals = data.loc[gene, group1_samples].dropna()
            group2_vals = data.loc[gene, group2_samples].dropna()
            
            # Skip genes with insufficient data or low expression
            if (len(group1_vals) < 2 or len(group2_vals) < 2 or 
                group1_vals.mean() < min_expression and group2_vals.mean() < min_expression):
                continue
            
            # Perform t-test
            try:
                t_stat, p_val = stats.ttest_ind(group1_vals, group2_vals, 
                                              equal_var=equal_var)
                
                # Calculate fold change and log2 fold change
                mean1 = group1_vals.mean()
                mean2 = group2_vals.mean()
                fold_change = mean1 / mean2 if mean2 > 0 else float('inf')
                log2_fc = np.log2(fold_change) if fold_change > 0 and fold_change != float('inf') else np.nan
                
                results.append({
                    'gene': gene,
                    'group1_mean': mean1,
                    'group2_mean': mean2,
                    'fold_change': fold_change,
                    'log2_fold_change': log2_fc,
                    't_statistic': t_stat,
                    'p_value': p_val
                })
                
            except Exception as e:
                continue
        
        results_df = pd.DataFrame(results)
        
        if len(results_df) > 0:
            # Multiple testing correction (Benjamini-Hochberg)
            from statsmodels.stats.multitest import multipletests
            _, results_df['p_adj'], _, _ = multipletests(results_df['p_value'], 
                                                       method='fdr_bh')
            
            # Sort by adjusted p-value
            results_df = results_df.sort_values('p_adj')
        
        print(f"Performed differential expression analysis on {len(results_df)} genes")
        return results_df
    
    def perform_deseq2_like_analysis(self, data, group1_samples, group2_samples):
        """
        Perform DESeq2-like analysis using negative binomial model approximation.
        
        Args:
            data (pd.DataFrame): Count data
            group1_samples (list): Sample names for group 1
            group2_samples (list): Sample names for group 2
            
        Returns:
            pd.DataFrame: Results with DESeq2-like statistics
        """
        results = []
        
        for gene in data.index:
            group1_counts = data.loc[gene, group1_samples].dropna()
            group2_counts = data.loc[gene, group2_samples].dropna()
            
            if len(group1_counts) < 2 or len(group2_counts) < 2:
                continue
            
            # Calculate means and dispersions
            mean1 = group1_counts.mean()
            mean2 = group2_counts.mean()
            
            # Skip genes with very low expression
            if mean1 < 1 and mean2 < 1:
                continue
            
            # Estimate dispersion (simplified)
            var1 = group1_counts.var()
            var2 = group2_counts.var()
            
            # Negative binomial test approximation
            try:
                # Use Wald test approximation
                pooled_mean = (mean1 + mean2) / 2
                pooled_var = (var1 + var2) / 2
                
                # Calculate log2 fold change
                log2_fc = np.log2((mean1 + 0.1) / (mean2 + 0.1))
                
                # Standard error estimation
                se = np.sqrt(1/len(group1_counts) + 1/len(group2_counts)) * np.sqrt(pooled_var)
                
                # Wald statistic
                wald_stat = log2_fc / se if se > 0 else 0
                p_val = 2 * (1 - stats.norm.cdf(abs(wald_stat)))
                
                results.append({
                    'gene': gene,
                    'baseMean': pooled_mean,
                    'log2FoldChange': log2_fc,
                    'lfcSE': se,
                    'stat': wald_stat,
                    'pvalue': p_val,
                    'group1_mean': mean1,
                    'group2_mean': mean2
                })
                
            except Exception as e:
                continue
        
        results_df = pd.DataFrame(results)
        
        if len(results_df) > 0:
            # Multiple testing correction
            from statsmodels.stats.multitest import multipletests
            _, results_df['padj'], _, _ = multipletests(results_df['pvalue'], 
                                                      method='fdr_bh')
            
            # Sort by adjusted p-value
            results_df = results_df.sort_values('padj')
        
        print(f"Performed DESeq2-like analysis on {len(results_df)} genes")
        return results_df
    
    def identify_top_genes(self, results, p_threshold=0.05, fc_threshold=1.5, 
                          n_top=None):
        """
        Identify top differentially expressed genes.
        
        Args:
            results (pd.DataFrame): Results from differential expression analysis
            p_threshold (float): P-value threshold
            fc_threshold (float): Fold change threshold
            n_top (int): Number of top genes to return
            
        Returns:
            pd.DataFrame: Top differentially expressed genes
        """
        # Filter by significance and fold change
        if 'padj' in results.columns:
            p_col = 'padj'
        else:
            p_col = 'p_value'
        
        if 'log2FoldChange' in results.columns:
            fc_col = 'log2FoldChange'
        else:
            fc_col = 'log2_fold_change'
        
        significant = results[
            (results[p_col] < p_threshold) & 
            (abs(results[fc_col]) > np.log2(fc_threshold))
        ].copy()
        
        if n_top is not None:
            significant = significant.head(n_top)
        
        print(f"Identified {len(significant)} significantly differentially expressed genes")
        return significant


class PathwayEnrichmentAnalyzer:
    """
    Class for functional analysis and pathway enrichment.
    
    Provides methods for Gene Ontology analysis, pathway enrichment,
    and gene set enrichment analysis.
    """
    
    def __init__(self):
        """Initialize pathway enrichment analyzer."""
        self.gene_sets = {}
        self.go_terms = {}
    
    def load_gene_sets(self, gene_sets_dict=None):
        """
        Load gene sets for enrichment analysis.
        
        Args:
            gene_sets_dict (dict): Dictionary of gene sets
        """
        if gene_sets_dict is None:
            # Load default gene sets (simplified for demo)
            self.gene_sets = {
                'Cell Cycle': ['CCNA2', 'CCNB1', 'CDK1', 'CDK2', 'PCNA', 'MKI67'],
                'Apoptosis': ['TP53', 'BAX', 'BCL2', 'CASP3', 'CASP9', 'FAS'],
                'Immune Response': ['IL1B', 'TNF', 'IFNG', 'IL6', 'CD8A', 'CD4'],
                'Metabolism': ['GAPDH', 'PKM', 'LDHA', 'HK1', 'PFKL', 'ALDOA'],
                'DNA Repair': ['BRCA1', 'BRCA2', 'ATM', 'ATR', 'TP53BP1', 'RAD51'],
                'Transcription': ['MYC', 'JUN', 'FOS', 'ETS1', 'SP1', 'STAT1'],
                'Signal Transduction': ['EGFR', 'MAPK1', 'AKT1', 'PIK3CA', 'MTOR', 'RAS']
            }
        else:
            self.gene_sets = gene_sets_dict
        
        print(f"Loaded {len(self.gene_sets)} gene sets")
    
    def perform_overrepresentation_analysis(self, gene_list, background_genes=None):
        """
        Perform gene set overrepresentation analysis.
        
        Args:
            gene_list (list): List of genes of interest
            background_genes (list): Background gene list
            
        Returns:
            pd.DataFrame: Enrichment results
        """
        if not self.gene_sets:
            self.load_gene_sets()
        
        if background_genes is None:
            # Use all genes in gene sets as background
            background_genes = set()
            for genes in self.gene_sets.values():
                background_genes.update(genes)
            background_genes = list(background_genes)
        
        gene_set = set(gene_list)
        background_set = set(background_genes)
        
        results = []
        
        for pathway_name, pathway_genes in self.gene_sets.items():
            pathway_set = set(pathway_genes) & background_set
            
            # Calculate overlap
            overlap = gene_set & pathway_set
            
            # Fisher's exact test
            a = len(overlap)  # genes in both lists
            b = len(pathway_set) - a  # genes in pathway but not in gene list
            c = len(gene_set) - a  # genes in gene list but not in pathway
            d = len(background_set) - len(gene_set) - b  # genes in neither
            
            if a > 0:
                # Fisher's exact test
                oddsratio, p_value = stats.fisher_exact([[a, b], [c, d]])
                
                results.append({
                    'pathway': pathway_name,
                    'pathway_size': len(pathway_set),
                    'overlap_size': a,
                    'gene_list_size': len(gene_set),
                    'background_size': len(background_set),
                    'odds_ratio': oddsratio,
                    'p_value': p_value,
                    'overlap_genes': list(overlap)
                })
        
        results_df = pd.DataFrame(results)
        
        if len(results_df) > 0:
            # Multiple testing correction
            from statsmodels.stats.multitest import multipletests
            _, results_df['p_adj'], _, _ = multipletests(results_df['p_value'], 
                                                       method='fdr_bh')
            
            # Sort by p-value
            results_df = results_df.sort_values('p_value')
        
        print(f"Performed overrepresentation analysis on {len(results_df)} pathways")
        return results_df
    
    def perform_gsea(self, ranked_genes, min_set_size=5, max_set_size=500):
        """
        Perform Gene Set Enrichment Analysis (GSEA).
        
        Args:
            ranked_genes (pd.Series): Genes ranked by metric (e.g., log2FC)
            min_set_size (int): Minimum gene set size
            max_set_size (int): Maximum gene set size
            
        Returns:
            pd.DataFrame: GSEA results
        """
        if not self.gene_sets:
            self.load_gene_sets()
        
        results = []
        
        for pathway_name, pathway_genes in self.gene_sets.items():
            # Filter gene set
            pathway_set = set(pathway_genes) & set(ranked_genes.index)
            
            if len(pathway_set) < min_set_size or len(pathway_set) > max_set_size:
                continue
            
            # Calculate enrichment score
            es_score = self._calculate_enrichment_score(ranked_genes, pathway_set)
            
            # Simple p-value estimation (for demo purposes)
            # In practice, would use permutation testing
            p_value = self._estimate_gsea_pvalue(ranked_genes, pathway_set, es_score)
            
            results.append({
                'pathway': pathway_name,
                'es_score': es_score,
                'p_value': p_value,
                'gene_set_size': len(pathway_set)
            })
        
        results_df = pd.DataFrame(results)
        
        if len(results_df) > 0:
            results_df = results_df.sort_values('p_value')
        
        print(f"Performed GSEA on {len(results_df)} pathways")
        return results_df
    
    def _calculate_enrichment_score(self, ranked_genes, gene_set):
        """Calculate GSEA enrichment score."""
        # Simplified enrichment score calculation
        ranks = []
        for gene in gene_set:
            if gene in ranked_genes.index:
                rank = list(ranked_genes.index).index(gene)
                ranks.append(rank)
        
        if not ranks:
            return 0
        
        # Calculate mean rank difference from expected
        n_genes = len(ranked_genes)
        expected_rank = n_genes / 2
        mean_rank = np.mean(ranks)
        
        # Normalize by total genes
        es_score = (expected_rank - mean_rank) / n_genes
        
        return es_score
    
    def _estimate_gsea_pvalue(self, ranked_genes, gene_set, observed_es):
        """Estimate GSEA p-value using simplified method."""
        # Simplified p-value estimation
        # In practice, would use proper permutation testing
        n_permutations = 100
        null_scores = []
        
        for _ in range(n_permutations):
            # Random gene set of same size
            random_genes = np.random.choice(ranked_genes.index, 
                                          size=len(gene_set), 
                                          replace=False)
            random_set = set(random_genes)
            null_score = self._calculate_enrichment_score(ranked_genes, random_set)
            null_scores.append(null_score)
        
        # Calculate p-value
        if observed_es >= 0:
            p_value = np.mean([score >= observed_es for score in null_scores])
        else:
            p_value = np.mean([score <= observed_es for score in null_scores])
        
        return max(p_value, 1/n_permutations)  # Minimum p-value


class TranscriptomeVisualizer:
    """
    Visualization tools for transcriptomic data.
    
    Provides methods for creating publication-quality plots of
    gene expression data and analysis results.
    """
    
    def __init__(self):
        """Initialize visualizer with custom styling."""
        plt.style.use('seaborn-v0_8')
        self.colors = {
            'upregulated': '#d62728',
            'downregulated': '#2ca02c',
            'not_significant': '#808080',
            'significant': '#1f77b4'
        }
    
    def plot_volcano(self, results, p_col='padj', fc_col='log2FoldChange',
                    p_threshold=0.05, fc_threshold=1.0, figsize=(10, 8)):
        """
        Create volcano plot for differential expression results.
        
        Args:
            results (pd.DataFrame): DE results
            p_col (str): P-value column name
            fc_col (str): Fold change column name
            p_threshold (float): P-value threshold
            fc_threshold (float): Log2 fold change threshold
            figsize (tuple): Figure size
            
        Returns:
            matplotlib.figure.Figure: The created figure
        """
        fig, ax = plt.subplots(figsize=figsize)
        
        # Prepare data
        x = results[fc_col]
        y = -np.log10(results[p_col])
        
        # Color points based on significance
        colors = []
        for i, row in results.iterrows():
            if row[p_col] < p_threshold and abs(row[fc_col]) > fc_threshold:
                if row[fc_col] > 0:
                    colors.append(self.colors['upregulated'])
                else:
                    colors.append(self.colors['downregulated'])
            else:
                colors.append(self.colors['not_significant'])
        
        # Create scatter plot
        scatter = ax.scatter(x, y, c=colors, alpha=0.6, s=20)
        
        # Add threshold lines
        ax.axhline(y=-np.log10(p_threshold), color='black', linestyle='--', alpha=0.5)
        ax.axvline(x=fc_threshold, color='black', linestyle='--', alpha=0.5)
        ax.axvline(x=-fc_threshold, color='black', linestyle='--', alpha=0.5)
        
        # Labels and title
        ax.set_xlabel(f'Log2 Fold Change ({fc_col})')
        ax.set_ylabel(f'-Log10 P-value ({p_col})')
        ax.set_title('Volcano Plot - Differential Gene Expression')
        
        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor=self.colors['upregulated'], label='Upregulated'),
            Patch(facecolor=self.colors['downregulated'], label='Downregulated'),
            Patch(facecolor=self.colors['not_significant'], label='Not Significant')
        ]
        ax.legend(handles=legend_elements)
        
        plt.tight_layout()
        return fig
    
    def plot_heatmap(self, data, sample_metadata=None, top_n=50, 
                    cluster_samples=True, cluster_genes=True, 
                    figsize=(12, 10)):
        """
        Create expression heatmap.
        
        Args:
            data (pd.DataFrame): Expression data (genes x samples)
            sample_metadata (pd.DataFrame): Sample annotations
            top_n (int): Number of top variable genes to plot
            cluster_samples (bool): Whether to cluster samples
            cluster_genes (bool): Whether to cluster genes
            figsize (tuple): Figure size
            
        Returns:
            matplotlib.figure.Figure: The created figure
        """
        # Select top variable genes
        gene_var = data.var(axis=1)
        top_genes = gene_var.nlargest(top_n).index
        plot_data = data.loc[top_genes]
        
        # Z-score normalization
        plot_data = plot_data.sub(plot_data.mean(axis=1), axis=0).div(plot_data.std(axis=1), axis=0)
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Create heatmap
        sns.heatmap(plot_data, 
                   cmap='RdBu_r', 
                   center=0,
                   cbar_kws={'label': 'Z-score'},
                   ax=ax,
                   xticklabels=True,
                   yticklabels=True)
        
        ax.set_title(f'Expression Heatmap - Top {top_n} Variable Genes')
        ax.set_xlabel('Samples')
        ax.set_ylabel('Genes')
        
        plt.tight_layout()
        return fig
    
    def plot_pca(self, data, metadata=None, color_by=None, figsize=(10, 8)):
        """
        Create PCA plot of samples.
        
        Args:
            data (pd.DataFrame): Expression data (genes x samples)
            metadata (pd.DataFrame): Sample metadata
            color_by (str): Column in metadata to color by
            figsize (tuple): Figure size
            
        Returns:
            matplotlib.figure.Figure: The created figure
        """
        # Prepare data (samples as rows)
        pca_data = data.T.fillna(0)
        
        # Standardize
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
            ax.legend()
        else:
            ax.scatter(pca_result[:, 0], pca_result[:, 1], alpha=0.7, s=50)
        
        ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
        ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
        ax.set_title('PCA of Gene Expression Data')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        return fig
    
    def plot_pathway_enrichment(self, enrichment_results, top_n=15, figsize=(10, 8)):
        """
        Create pathway enrichment plot.
        
        Args:
            enrichment_results (pd.DataFrame): Pathway enrichment results
            top_n (int): Number of top pathways to plot
            figsize (tuple): Figure size
            
        Returns:
            matplotlib.figure.Figure: The created figure
        """
        # Select top pathways
        plot_data = enrichment_results.head(top_n).copy()
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Create horizontal bar plot
        y_pos = np.arange(len(plot_data))
        bars = ax.barh(y_pos, -np.log10(plot_data['p_value']), 
                      color=self.colors['significant'], alpha=0.7)
        
        # Customize plot
        ax.set_yticks(y_pos)
        ax.set_yticklabels(plot_data['pathway'])
        ax.set_xlabel('-Log10 P-value')
        ax.set_title('Pathway Enrichment Analysis')
        ax.grid(True, alpha=0.3, axis='x')
        
        # Add significance threshold line
        ax.axvline(x=-np.log10(0.05), color='red', linestyle='--', alpha=0.5)
        
        plt.tight_layout()
        return fig
    
    def plot_expression_boxplot(self, data, genes, metadata=None, 
                              group_by=None, figsize=(12, 6)):
        """
        Create boxplot of gene expression across groups.
        
        Args:
            data (pd.DataFrame): Expression data
            genes (list): List of genes to plot
            metadata (pd.DataFrame): Sample metadata
            group_by (str): Column to group samples by
            figsize (tuple): Figure size
            
        Returns:
            matplotlib.figure.Figure: The created figure
        """
        n_genes = len(genes)
        fig, axes = plt.subplots(1, n_genes, figsize=figsize, sharey=True)
        
        if n_genes == 1:
            axes = [axes]
        
        for i, gene in enumerate(genes):
            if gene not in data.index:
                axes[i].text(0.5, 0.5, f'Gene {gene}\nnot found', 
                           ha='center', va='center', transform=axes[i].transAxes)
                continue
            
            gene_data = data.loc[gene]
            
            if metadata is not None and group_by is not None:
                # Create DataFrame for plotting
                plot_df = pd.DataFrame({
                    'expression': gene_data,
                    'group': metadata[group_by]
                })
                
                sns.boxplot(data=plot_df, x='group', y='expression', ax=axes[i])
            else:
                axes[i].boxplot(gene_data.dropna())
            
            axes[i].set_title(gene)
            axes[i].set_ylabel('Expression' if i == 0 else '')
            axes[i].tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        return fig 