"""
Proteomic Data Analysis Utilities

This module provides comprehensive tools for analyzing mass spectrometry-based proteomics data,
including protein quantification, differential expression analysis, and functional interpretation.

Classes:
    - ProteomeDataProcessor: Core data processing and protein quantification
    - DifferentialProteinAnalyzer: Statistical analysis of protein expression changes
    - ProteinPathwayAnalyzer: Functional analysis and pathway enrichment
    - ProteomicsVisualizer: Visualization tools for proteomics data

Author: Stanford Data Ocean - Data Science in Precision Medicine
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
import warnings
warnings.filterwarnings('ignore')


class ProteomeDataProcessor:
    """
    Core class for processing and analyzing proteomics datasets.
    
    This class handles loading, processing, and quality control of mass spectrometry
    data including protein abundance matrices, sample metadata, and protein annotations.
    """
    
    def __init__(self):
        """Initialize the proteome data processor."""
        self.abundance_data = None
        self.sample_metadata = None
        self.protein_annotations = None
        self.missing_value_threshold = 0.5
        
    def load_abundance_matrix(self, filepath, protein_id_col=0, sample_start_col=1):
        """
        Load protein abundance matrix from mass spectrometry data.
        
        Args:
            filepath (str): Path to the abundance matrix file
            protein_id_col (int): Column index for protein IDs
            sample_start_col (int): Starting column index for sample data
            
        Returns:
            pd.DataFrame: Processed abundance matrix
        """
        try:
            data = pd.read_csv(filepath, sep='\t', index_col=protein_id_col)
            
            # Remove non-sample columns if present
            non_sample_cols = ['Gene.names', 'Protein.names', 'Fasta.headers']
            for col in non_sample_cols:
                if col in data.columns:
                    data = data.drop(col, axis=1)
            
            # Convert to numeric, handling missing values
            data = data.apply(pd.to_numeric, errors='coerce')
            
            # Log transform if values appear to be raw intensities
            if data.max().max() > 100:
                data = np.log2(data + 1)
                print("Applied log2 transformation to abundance data")
            
            self.abundance_data = data
            print(f"Loaded abundance matrix: {data.shape[0]} proteins, {data.shape[1]} samples")
            return data
            
        except Exception as e:
            print(f"Error loading abundance matrix: {e}")
            return None
    
    def load_sample_metadata(self, filepath=None, metadata_dict=None):
        """
        Load sample metadata for experimental design.
        
        Args:
            filepath (str): Path to metadata file
            metadata_dict (dict): Dictionary of sample metadata
            
        Returns:
            pd.DataFrame: Sample metadata
        """
        try:
            if filepath:
                metadata = pd.read_csv(filepath, sep='\t', index_col=0)
            elif metadata_dict:
                metadata = pd.DataFrame.from_dict(metadata_dict, orient='index')
            else:
                raise ValueError("Either filepath or metadata_dict must be provided")
            
            # Convert categorical columns
            for col in metadata.columns:
                if metadata[col].dtype == 'object':
                    metadata[col] = metadata[col].astype('category')
            
            self.sample_metadata = metadata
            print(f"Loaded metadata for {len(metadata)} samples")
            return metadata
            
        except Exception as e:
            print(f"Error loading sample metadata: {e}")
            return None
    
    def handle_missing_values(self, data=None, method='knn', missing_threshold=0.5):
        """
        Handle missing values in proteomics data.
        
        Args:
            data (pd.DataFrame): Abundance matrix
            method (str): Method for handling missing values ('knn', 'min', 'median', 'drop')
            missing_threshold (float): Threshold for dropping proteins with too many missing values
            
        Returns:
            pd.DataFrame: Data with missing values handled
        """
        if data is None:
            data = self.abundance_data
        
        if data is None:
            raise ValueError("No abundance data available")
        
        # Calculate missing value percentage per protein
        missing_pct = data.isnull().sum(axis=1) / data.shape[1]
        
        # Drop proteins with too many missing values
        keep_proteins = missing_pct <= missing_threshold
        filtered_data = data[keep_proteins]
        
        print(f"Removed {(~keep_proteins).sum()} proteins with >{missing_threshold*100}% missing values")
        
        if method == 'knn':
            # Simple KNN imputation using row means
            for protein in filtered_data.index:
                missing_mask = filtered_data.loc[protein].isnull()
                if missing_mask.any():
                    # Use median of available values for imputation
                    impute_value = filtered_data.loc[protein].median()
                    if pd.isna(impute_value):
                        # If all values missing, use global median
                        impute_value = filtered_data.median().median()
                    filtered_data.loc[protein, missing_mask] = impute_value
                    
        elif method == 'min':
            # Replace with minimum value (common for missing at random)
            global_min = filtered_data.min().min()
            filtered_data = filtered_data.fillna(global_min)
            
        elif method == 'median':
            # Replace with median value
            filtered_data = filtered_data.fillna(filtered_data.median())
            
        elif method == 'drop':
            # Drop any remaining missing values
            filtered_data = filtered_data.dropna()
        
        print(f"Processed missing values using {method} method")
        return filtered_data
    
    def normalize_abundance_data(self, data=None, method='median_centering'):
        """
        Normalize protein abundance data.
        
        Args:
            data (pd.DataFrame): Abundance matrix
            method (str): Normalization method ('median_centering', 'quantile', 'zscore')
            
        Returns:
            pd.DataFrame: Normalized abundance data
        """
        if data is None:
            data = self.abundance_data
        
        if data is None:
            raise ValueError("No abundance data available")
        
        if method == 'median_centering':
            # Subtract median from each sample
            sample_medians = data.median(axis=0)
            normalized = data.sub(sample_medians, axis=1)
            
        elif method == 'quantile':
            # Quantile normalization
            rank_mean = data.stack().groupby(data.rank(method='first').stack().astype(int)).mean()
            normalized = data.rank(method='min').stack().astype(int).map(rank_mean).unstack()
            
        elif method == 'zscore':
            # Z-score normalization per sample
            normalized = data.apply(lambda x: (x - x.mean()) / x.std(), axis=0)
            
        else:
            raise ValueError(f"Unknown normalization method: {method}")
        
        print(f"Normalized data using {method} method")
        return normalized
    
    def calculate_protein_statistics(self, data=None):
        """
        Calculate basic statistics for proteins.
        
        Args:
            data (pd.DataFrame): Abundance matrix
            
        Returns:
            pd.DataFrame: Protein statistics
        """
        if data is None:
            data = self.abundance_data
        
        if data is None:
            raise ValueError("No abundance data available")
        
        stats_df = pd.DataFrame(index=data.index)
        
        # Basic statistics
        stats_df['mean_abundance'] = data.mean(axis=1)
        stats_df['median_abundance'] = data.median(axis=1)
        stats_df['std_abundance'] = data.std(axis=1)
        stats_df['cv'] = stats_df['std_abundance'] / stats_df['mean_abundance']
        
        # Missing value statistics
        stats_df['missing_values'] = data.isnull().sum(axis=1)
        stats_df['missing_pct'] = stats_df['missing_values'] / data.shape[1]
        
        # Detection frequency
        stats_df['detection_freq'] = (data > 0).sum(axis=1) / data.shape[1]
        
        return stats_df
    
    def filter_proteins(self, data=None, min_detection_freq=0.5, max_cv=1.0):
        """
        Filter proteins based on quality criteria.
        
        Args:
            data (pd.DataFrame): Abundance matrix
            min_detection_freq (float): Minimum detection frequency
            max_cv (float): Maximum coefficient of variation
            
        Returns:
            pd.DataFrame: Filtered protein data
        """
        if data is None:
            data = self.abundance_data
        
        stats = self.calculate_protein_statistics(data)
        
        # Apply filters
        keep_proteins = (
            (stats['detection_freq'] >= min_detection_freq) &
            (stats['cv'] <= max_cv) &
            (~stats['cv'].isna())
        )
        
        filtered_data = data[keep_proteins]
        
        print(f"Protein filtering: {data.shape[0]} -> {filtered_data.shape[0]} proteins")
        print(f"  Detection frequency >= {min_detection_freq}: {(stats['detection_freq'] >= min_detection_freq).sum()} proteins")
        print(f"  CV <= {max_cv}: {(stats['cv'] <= max_cv).sum()} proteins")
        
        return filtered_data


class DifferentialProteinAnalyzer:
    """
    Class for statistical analysis of differential protein expression.
    
    Provides methods for identifying differentially expressed proteins
    between conditions using various statistical approaches.
    """
    
    def __init__(self):
        """Initialize differential protein expression analyzer."""
        pass
    
    def perform_ttest_analysis(self, data, group1_samples, group2_samples, 
                              paired=False, equal_var=False):
        """
        Perform t-test for differential protein expression.
        
        Args:
            data (pd.DataFrame): Abundance data (proteins x samples)
            group1_samples (list): Sample names for group 1
            group2_samples (list): Sample names for group 2
            paired (bool): Whether to perform paired t-test
            equal_var (bool): Assume equal variances
            
        Returns:
            pd.DataFrame: Results with statistics and p-values
        """
        results = []
        
        for protein in data.index:
            group1_vals = data.loc[protein, group1_samples].dropna()
            group2_vals = data.loc[protein, group2_samples].dropna()
            
            # Skip proteins with insufficient data
            if len(group1_vals) < 2 or len(group2_vals) < 2:
                continue
            
            try:
                if paired:
                    # Paired t-test
                    if len(group1_vals) == len(group2_vals):
                        t_stat, p_val = stats.ttest_rel(group1_vals, group2_vals)
                    else:
                        continue  # Skip if unequal sample sizes for paired test
                else:
                    # Independent t-test
                    t_stat, p_val = stats.ttest_ind(group1_vals, group2_vals, 
                                                  equal_var=equal_var)
                
                # Calculate effect size and fold change
                mean1 = group1_vals.mean()
                mean2 = group2_vals.mean()
                
                # Calculate fold change (avoiding log of negative values)
                if mean1 > 0 and mean2 > 0:
                    fold_change = mean1 / mean2
                    log2_fc = np.log2(fold_change)
                else:
                    fold_change = np.nan
                    log2_fc = mean1 - mean2  # Use difference for log-transformed data
                
                # Cohen's d (effect size)
                pooled_std = np.sqrt(((len(group1_vals)-1)*group1_vals.var() + 
                                    (len(group2_vals)-1)*group2_vals.var()) / 
                                   (len(group1_vals) + len(group2_vals) - 2))
                cohens_d = (mean1 - mean2) / pooled_std if pooled_std > 0 else 0
                
                results.append({
                    'protein': protein,
                    'group1_mean': mean1,
                    'group2_mean': mean2,
                    'fold_change': fold_change,
                    'log2_fold_change': log2_fc,
                    't_statistic': t_stat,
                    'p_value': p_val,
                    'cohens_d': cohens_d,
                    'group1_n': len(group1_vals),
                    'group2_n': len(group2_vals)
                })
                
            except Exception as e:
                continue
        
        results_df = pd.DataFrame(results)
        
        if len(results_df) > 0:
            # Multiple testing correction
            from statsmodels.stats.multitest import multipletests
            _, results_df['p_adj'], _, _ = multipletests(results_df['p_value'], 
                                                       method='fdr_bh')
            
            # Sort by adjusted p-value
            results_df = results_df.sort_values('p_adj')
        
        print(f"Performed differential protein analysis on {len(results_df)} proteins")
        return results_df
    
    def perform_limma_like_analysis(self, data, design_matrix, contrasts):
        """
        Perform limma-like analysis for complex experimental designs.
        
        Args:
            data (pd.DataFrame): Abundance data
            design_matrix (pd.DataFrame): Design matrix for linear modeling
            contrasts (dict): Contrast definitions
            
        Returns:
            dict: Results for each contrast
        """
        from sklearn.linear_model import LinearRegression
        from scipy.stats import t
        
        results = {}
        
        for contrast_name, contrast_vector in contrasts.items():
            contrast_results = []
            
            for protein in data.index:
                protein_data = data.loc[protein].dropna()
                
                # Align design matrix with available data
                available_samples = protein_data.index
                aligned_design = design_matrix.loc[available_samples]
                aligned_data = protein_data.loc[available_samples]
                
                if len(aligned_data) < 3:  # Need minimum samples
                    continue
                
                try:
                    # Fit linear model
                    model = LinearRegression()
                    model.fit(aligned_design, aligned_data)
                    
                    # Calculate residuals and standard error
                    predictions = model.predict(aligned_design)
                    residuals = aligned_data - predictions
                    mse = np.sum(residuals**2) / (len(aligned_data) - aligned_design.shape[1])
                    
                    # Calculate contrast
                    contrast_coef = np.dot(contrast_vector, model.coef_)
                    
                    # Standard error of contrast
                    contrast_var = mse * np.dot(contrast_vector, 
                                              np.dot(np.linalg.pinv(np.dot(aligned_design.T, aligned_design)), 
                                                   contrast_vector))
                    contrast_se = np.sqrt(contrast_var)
                    
                    # T-statistic and p-value
                    t_stat = contrast_coef / contrast_se if contrast_se > 0 else 0
                    df = len(aligned_data) - aligned_design.shape[1]
                    p_val = 2 * (1 - t.cdf(abs(t_stat), df))
                    
                    contrast_results.append({
                        'protein': protein,
                        'log2_fold_change': contrast_coef,
                        't_statistic': t_stat,
                        'p_value': p_val,
                        'standard_error': contrast_se,
                        'degrees_freedom': df
                    })
                    
                except Exception as e:
                    continue
            
            contrast_df = pd.DataFrame(contrast_results)
            
            if len(contrast_df) > 0:
                # Multiple testing correction
                from statsmodels.stats.multitest import multipletests
                _, contrast_df['p_adj'], _, _ = multipletests(contrast_df['p_value'], 
                                                            method='fdr_bh')
                contrast_df = contrast_df.sort_values('p_adj')
            
            results[contrast_name] = contrast_df
            print(f"Analyzed {len(contrast_df)} proteins for contrast {contrast_name}")
        
        return results
    
    def identify_significant_proteins(self, results, p_threshold=0.05, 
                                    fc_threshold=1.5, effect_size_threshold=0.5):
        """
        Identify significantly differentially expressed proteins.
        
        Args:
            results (pd.DataFrame): Results from differential analysis
            p_threshold (float): P-value threshold
            fc_threshold (float): Fold change threshold
            effect_size_threshold (float): Effect size threshold
            
        Returns:
            pd.DataFrame: Significant proteins
        """
        # Determine column names
        p_col = 'p_adj' if 'p_adj' in results.columns else 'p_value'
        fc_col = 'log2_fold_change' if 'log2_fold_change' in results.columns else 'log2_fold_change'
        
        # Apply filters
        significant = results[
            (results[p_col] < p_threshold) & 
            (abs(results[fc_col]) > np.log2(fc_threshold))
        ].copy()
        
        # Add effect size filter if available
        if 'cohens_d' in results.columns:
            significant = significant[abs(significant['cohens_d']) > effect_size_threshold]
        
        # Classify as up/down regulated
        significant['regulation'] = significant[fc_col].apply(
            lambda x: 'upregulated' if x > 0 else 'downregulated'
        )
        
        print(f"Identified {len(significant)} significantly differentially expressed proteins")
        if len(significant) > 0:
            up_count = (significant['regulation'] == 'upregulated').sum()
            down_count = (significant['regulation'] == 'downregulated').sum()
            print(f"  Upregulated: {up_count}, Downregulated: {down_count}")
        
        return significant


class ProteinPathwayAnalyzer:
    """
    Class for functional analysis and pathway enrichment of proteins.
    
    Provides methods for Gene Ontology analysis, pathway enrichment,
    and protein-protein interaction analysis.
    """
    
    def __init__(self):
        """Initialize protein pathway analyzer."""
        self.protein_pathways = {}
        self.go_annotations = {}
        
    def load_pathway_databases(self, pathway_dict=None):
        """
        Load pathway databases for enrichment analysis.
        
        Args:
            pathway_dict (dict): Dictionary of pathway definitions
        """
        if pathway_dict is None:
            # Load default pathways relevant to proteomics
            self.protein_pathways = {
                'Protein Folding': ['HSP90AA1', 'HSP70', 'HSPA1A', 'CALR', 'PDIA3', 'CANX'],
                'Protein Degradation': ['PSMA1', 'PSMA2', 'PSMB1', 'PSMB2', 'UBE2I', 'UBC'],
                'Metabolic Enzymes': ['GAPDH', 'PKM', 'LDHA', 'ENO1', 'PGK1', 'ALDOA'],
                'Immune Response': ['C3', 'C4A', 'CFD', 'HP', 'SAA1', 'CRP'],
                'Blood Coagulation': ['FGA', 'FGB', 'FGG', 'F2', 'F9', 'SERPINC1'],
                'Lipid Metabolism': ['APOA1', 'APOA2', 'APOB', 'APOE', 'LCAT', 'CETP'],
                'Oxidative Stress': ['SOD1', 'SOD2', 'CAT', 'GPX1', 'PRDX1', 'TXNRD1'],
                'Cell Adhesion': ['ITGA1', 'ITGB1', 'CDH1', 'VCAM1', 'ICAM1', 'PECAM1'],
                'Signal Transduction': ['EGFR', 'MAPK1', 'AKT1', 'PIK3CA', 'PTEN', 'RAF1'],
                'Cytoskeleton': ['ACTB', 'TUBB', 'VIM', 'KRT1', 'KRT10', 'ACTN1']
            }
        else:
            self.protein_pathways = pathway_dict
        
        print(f"Loaded {len(self.protein_pathways)} protein pathway databases")
    
    def perform_pathway_enrichment(self, protein_list, background_proteins=None):
        """
        Perform pathway enrichment analysis.
        
        Args:
            protein_list (list): List of proteins of interest
            background_proteins (list): Background protein list
            
        Returns:
            pd.DataFrame: Enrichment results
        """
        if not self.protein_pathways:
            self.load_pathway_databases()
        
        if background_proteins is None:
            # Use all proteins in pathways as background
            background_proteins = set()
            for proteins in self.protein_pathways.values():
                background_proteins.update(proteins)
            background_proteins = list(background_proteins)
        
        protein_set = set(protein_list)
        background_set = set(background_proteins)
        
        results = []
        
        for pathway_name, pathway_proteins in self.protein_pathways.items():
            pathway_set = set(pathway_proteins) & background_set
            
            if len(pathway_set) == 0:
                continue
            
            # Calculate overlap
            overlap = protein_set & pathway_set
            
            # Fisher's exact test
            a = len(overlap)  # proteins in both lists
            b = len(pathway_set) - a  # proteins in pathway but not in protein list
            c = len(protein_set) - a  # proteins in protein list but not in pathway
            d = len(background_set) - len(protein_set) - b  # proteins in neither
            
            if a > 0:
                # Fisher's exact test
                oddsratio, p_value = stats.fisher_exact([[a, b], [c, d]])
                
                # Calculate enrichment ratio
                expected = len(protein_set) * len(pathway_set) / len(background_set)
                enrichment_ratio = a / expected if expected > 0 else 0
                
                results.append({
                    'pathway': pathway_name,
                    'pathway_size': len(pathway_set),
                    'overlap_size': a,
                    'protein_list_size': len(protein_set),
                    'background_size': len(background_set),
                    'odds_ratio': oddsratio,
                    'enrichment_ratio': enrichment_ratio,
                    'p_value': p_value,
                    'overlap_proteins': list(overlap)
                })
        
        results_df = pd.DataFrame(results)
        
        if len(results_df) > 0:
            # Multiple testing correction
            from statsmodels.stats.multitest import multipletests
            _, results_df['p_adj'], _, _ = multipletests(results_df['p_value'], 
                                                       method='fdr_bh')
            
            # Sort by p-value
            results_df = results_df.sort_values('p_value')
        
        print(f"Performed pathway enrichment analysis on {len(results_df)} pathways")
        return results_df
    
    def analyze_protein_functions(self, protein_list):
        """
        Analyze functional categories of proteins.
        
        Args:
            protein_list (list): List of proteins to analyze
            
        Returns:
            pd.DataFrame: Functional analysis results
        """
        if not self.protein_pathways:
            self.load_pathway_databases()
        
        protein_functions = {}
        
        for protein in protein_list:
            protein_functions[protein] = []
            
            # Find pathways containing this protein
            for pathway_name, pathway_proteins in self.protein_pathways.items():
                if protein in pathway_proteins:
                    protein_functions[protein].append(pathway_name)
        
        # Create summary
        function_summary = []
        for protein, functions in protein_functions.items():
            function_summary.append({
                'protein': protein,
                'n_pathways': len(functions),
                'pathways': '; '.join(functions) if functions else 'Unknown'
            })
        
        return pd.DataFrame(function_summary)
    
    def create_protein_network(self, protein_list, interaction_threshold=0.7):
        """
        Create protein-protein interaction network (simplified).
        
        Args:
            protein_list (list): List of proteins
            interaction_threshold (float): Threshold for interactions
            
        Returns:
            dict: Network representation
        """
        # Simplified network based on pathway co-membership
        network = {'nodes': [], 'edges': []}
        
        # Add nodes
        for protein in protein_list:
            network['nodes'].append({
                'id': protein,
                'label': protein
            })
        
        # Add edges based on pathway co-membership
        for i, protein1 in enumerate(protein_list):
            for j, protein2 in enumerate(protein_list[i+1:], i+1):
                
                # Count shared pathways
                pathways1 = set()
                pathways2 = set()
                
                for pathway_name, pathway_proteins in self.protein_pathways.items():
                    if protein1 in pathway_proteins:
                        pathways1.add(pathway_name)
                    if protein2 in pathway_proteins:
                        pathways2.add(pathway_name)
                
                shared_pathways = pathways1 & pathways2
                
                if len(shared_pathways) > 0:
                    # Calculate interaction strength
                    strength = len(shared_pathways) / max(len(pathways1), len(pathways2))
                    
                    if strength >= interaction_threshold:
                        network['edges'].append({
                            'source': protein1,
                            'target': protein2,
                            'weight': strength,
                            'shared_pathways': list(shared_pathways)
                        })
        
        print(f"Created protein network with {len(network['nodes'])} nodes and {len(network['edges'])} edges")
        return network


class ProteomicsVisualizer:
    """
    Visualization tools for proteomics data.
    
    Provides methods for creating publication-quality plots of
    protein abundance data and analysis results.
    """
    
    def __init__(self):
        """Initialize visualizer with custom styling."""
        plt.style.use('seaborn-v0_8')
        self.colors = {
            'upregulated': '#d62728',
            'downregulated': '#2ca02c',
            'not_significant': '#808080',
            'significant': '#1f77b4',
            'missing': '#ffcc99'
        }
    
    def plot_volcano(self, results, p_col='p_adj', fc_col='log2_fold_change',
                    p_threshold=0.05, fc_threshold=1.0, figsize=(10, 8)):
        """
        Create volcano plot for differential protein expression.
        
        Args:
            results (pd.DataFrame): Differential expression results
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
        scatter = ax.scatter(x, y, c=colors, alpha=0.6, s=30)
        
        # Add threshold lines
        ax.axhline(y=-np.log10(p_threshold), color='black', linestyle='--', alpha=0.5)
        ax.axvline(x=fc_threshold, color='black', linestyle='--', alpha=0.5)
        ax.axvline(x=-fc_threshold, color='black', linestyle='--', alpha=0.5)
        
        # Labels and title
        ax.set_xlabel(f'Log2 Fold Change ({fc_col})')
        ax.set_ylabel(f'-Log10 P-value ({p_col})')
        ax.set_title('Volcano Plot - Differential Protein Expression')
        
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
    
    def plot_abundance_heatmap(self, data, sample_metadata=None, top_n=50, 
                              cluster_samples=True, cluster_proteins=True, 
                              figsize=(14, 10)):
        """
        Create protein abundance heatmap.
        
        Args:
            data (pd.DataFrame): Abundance data (proteins x samples)
            sample_metadata (pd.DataFrame): Sample annotations
            top_n (int): Number of top variable proteins to plot
            cluster_samples (bool): Whether to cluster samples
            cluster_proteins (bool): Whether to cluster proteins
            figsize (tuple): Figure size
            
        Returns:
            matplotlib.figure.Figure: The created figure
        """
        # Select top variable proteins
        protein_var = data.var(axis=1)
        top_proteins = protein_var.nlargest(top_n).index
        plot_data = data.loc[top_proteins]
        
        # Z-score normalization for visualization
        plot_data = plot_data.sub(plot_data.mean(axis=1), axis=0).div(plot_data.std(axis=1), axis=0)
        plot_data = plot_data.fillna(0)  # Handle any remaining NaN values
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Create heatmap
        sns.heatmap(plot_data, 
                   cmap='RdBu_r', 
                   center=0,
                   cbar_kws={'label': 'Z-score (abundance)'},
                   ax=ax,
                   xticklabels=True,
                   yticklabels=True)
        
        ax.set_title(f'Protein Abundance Heatmap - Top {top_n} Variable Proteins')
        ax.set_xlabel('Samples')
        ax.set_ylabel('Proteins')
        
        plt.tight_layout()
        return fig
    
    def plot_pca(self, data, metadata=None, color_by=None, figsize=(10, 8)):
        """
        Create PCA plot of samples based on protein abundance.
        
        Args:
            data (pd.DataFrame): Abundance data (proteins x samples)
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
                          label=category, alpha=0.7, s=60)
            ax.legend()
        else:
            ax.scatter(pca_result[:, 0], pca_result[:, 1], alpha=0.7, s=60)
        
        ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
        ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
        ax.set_title('PCA of Protein Abundance Data')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        return fig
    
    def plot_pathway_enrichment(self, enrichment_results, top_n=15, figsize=(12, 8)):
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
        
        # Color bars by enrichment ratio
        if 'enrichment_ratio' in plot_data.columns:
            for i, bar in enumerate(bars):
                ratio = plot_data.iloc[i]['enrichment_ratio']
                if ratio > 2:
                    bar.set_color(self.colors['upregulated'])
                elif ratio > 1.5:
                    bar.set_color(self.colors['significant'])
        
        # Customize plot
        ax.set_yticks(y_pos)
        ax.set_yticklabels(plot_data['pathway'])
        ax.set_xlabel('-Log10 P-value')
        ax.set_title('Protein Pathway Enrichment Analysis')
        ax.grid(True, alpha=0.3, axis='x')
        
        # Add significance threshold line
        ax.axvline(x=-np.log10(0.05), color='red', linestyle='--', alpha=0.5)
        
        plt.tight_layout()
        return fig
    
    def plot_abundance_distribution(self, data, sample_metadata=None, 
                                  group_by=None, figsize=(12, 6)):
        """
        Plot distribution of protein abundances.
        
        Args:
            data (pd.DataFrame): Abundance data
            sample_metadata (pd.DataFrame): Sample metadata
            group_by (str): Column to group samples by
            figsize (tuple): Figure size
            
        Returns:
            matplotlib.figure.Figure: The created figure
        """
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
        
        # Overall distribution
        all_values = data.values.flatten()
        all_values = all_values[~np.isnan(all_values)]
        
        ax1.hist(all_values, bins=50, alpha=0.7, color=self.colors['significant'])
        ax1.set_xlabel('Log2 Abundance')
        ax1.set_ylabel('Frequency')
        ax1.set_title('Overall Abundance Distribution')
        ax1.grid(True, alpha=0.3)
        
        # Sample-wise distributions
        if sample_metadata is not None and group_by is not None:
            for category in sample_metadata[group_by].unique():
                samples = sample_metadata[sample_metadata[group_by] == category].index
                category_data = data[samples].values.flatten()
                category_data = category_data[~np.isnan(category_data)]
                
                ax2.hist(category_data, bins=30, alpha=0.6, label=category)
            
            ax2.legend()
            ax2.set_title(f'Abundance Distribution by {group_by}')
        else:
            # Box plot of samples
            data_for_box = [data[col].dropna() for col in data.columns]
            ax2.boxplot(data_for_box)
            ax2.set_title('Abundance Distribution by Sample')
            ax2.set_xticklabels(data.columns, rotation=45)
        
        ax2.set_xlabel('Log2 Abundance')
        ax2.set_ylabel('Frequency')
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        return fig
    
    def plot_missing_values(self, data, figsize=(12, 8)):
        """
        Visualize missing value patterns in proteomics data.
        
        Args:
            data (pd.DataFrame): Abundance data
            figsize (tuple): Figure size
            
        Returns:
            matplotlib.figure.Figure: The created figure
        """
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize)
        
        # Missing values heatmap
        missing_data = data.isnull()
        
        # Sample subset for visualization if too many proteins
        if len(data) > 100:
            # Show proteins with most missing values
            missing_counts = missing_data.sum(axis=1)
            top_missing = missing_counts.nlargest(50).index
            plot_missing = missing_data.loc[top_missing]
        else:
            plot_missing = missing_data
        
        sns.heatmap(plot_missing, cmap='RdYlBu_r', cbar_kws={'label': 'Missing'}, 
                   ax=ax1, xticklabels=True, yticklabels=True)
        ax1.set_title('Missing Value Pattern (Top Missing Proteins)')
        ax1.set_xlabel('Samples')
        ax1.set_ylabel('Proteins')
        
        # Missing value statistics
        missing_per_sample = missing_data.sum(axis=0)
        missing_per_protein = missing_data.sum(axis=1)
        
        ax2_twin = ax2.twinx()
        
        # Bar plot for samples
        bars1 = ax2.bar(range(len(missing_per_sample)), missing_per_sample, 
                       alpha=0.7, color=self.colors['missing'], label='Missing per sample')
        ax2.set_xlabel('Sample Index')
        ax2.set_ylabel('Missing Values per Sample', color=self.colors['missing'])
        ax2.tick_params(axis='y', labelcolor=self.colors['missing'])
        
        # Histogram for proteins
        ax2_twin.hist(missing_per_protein, bins=20, alpha=0.5, 
                     color=self.colors['significant'], label='Distribution of missing per protein')
        ax2_twin.set_ylabel('Frequency', color=self.colors['significant'])
        ax2_twin.tick_params(axis='y', labelcolor=self.colors['significant'])
        
        ax2.set_title('Missing Value Statistics')
        
        plt.tight_layout()
        return fig
    
    def plot_protein_abundance_comparison(self, data, proteins, metadata=None, 
                                        group_by=None, figsize=(15, 8)):
        """
        Compare abundance of specific proteins across conditions.
        
        Args:
            data (pd.DataFrame): Abundance data
            proteins (list): List of proteins to compare
            metadata (pd.DataFrame): Sample metadata
            group_by (str): Column to group samples by
            figsize (tuple): Figure size
            
        Returns:
            matplotlib.figure.Figure: The created figure
        """
        n_proteins = len(proteins)
        fig, axes = plt.subplots(1, n_proteins, figsize=figsize, sharey=True)
        
        if n_proteins == 1:
            axes = [axes]
        
        for i, protein in enumerate(proteins):
            if protein not in data.index:
                axes[i].text(0.5, 0.5, f'Protein {protein}\nnot found', 
                           ha='center', va='center', transform=axes[i].transAxes)
                continue
            
            protein_data = data.loc[protein]
            
            if metadata is not None and group_by is not None:
                # Create DataFrame for plotting
                plot_df = pd.DataFrame({
                    'abundance': protein_data,
                    'group': metadata[group_by]
                })
                
                sns.boxplot(data=plot_df, x='group', y='abundance', ax=axes[i])
                sns.swarmplot(data=plot_df, x='group', y='abundance', ax=axes[i], 
                            color='black', alpha=0.6, size=4)
            else:
                axes[i].boxplot(protein_data.dropna())
            
            axes[i].set_title(protein)
            axes[i].set_ylabel('Log2 Abundance' if i == 0 else '')
            axes[i].tick_params(axis='x', rotation=45)
            axes[i].grid(True, alpha=0.3)
        
        plt.tight_layout()
        return fig 