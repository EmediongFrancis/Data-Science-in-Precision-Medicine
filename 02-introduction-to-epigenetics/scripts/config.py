"""
Configuration file for epigenetic data analysis.

Contains constants, sample data, and configuration parameters
for the Introduction to Epigenetics module.
"""

import pandas as pd
import numpy as np

# Sample metadata for NIH Roadmap Epigenomics samples
SAMPLE_METADATA = {
    'E001': {'cell_type': 'ES-I3', 'tissue': 'ESC', 'category': 'ESC'},
    'E002': {'cell_type': 'ES-WA7', 'tissue': 'ESC', 'category': 'ESC'},
    'E003': {'cell_type': 'H1', 'tissue': 'ESC', 'category': 'ESC'},
    'E004': {'cell_type': 'H1_BMP4', 'tissue': 'ESC_Derived', 'category': 'ESC_Derived'},
    'E005': {'cell_type': 'H1_Trophoblast', 'tissue': 'ESC_Derived', 'category': 'ESC_Derived'},
    'E006': {'cell_type': 'H1_Mesenchymal', 'tissue': 'ESC_Derived', 'category': 'ESC_Derived'},
    'E007': {'cell_type': 'H1_Neuronal', 'tissue': 'ESC_Derived', 'category': 'ESC_Derived'},
    'E008': {'cell_type': 'H9', 'tissue': 'ESC', 'category': 'ESC'},
    'E009': {'cell_type': 'H9_Neuronal', 'tissue': 'ESC_Derived', 'category': 'ESC_Derived'},
    'E010': {'cell_type': 'H9_Mesenchymal', 'tissue': 'ESC_Derived', 'category': 'ESC_Derived'},
    'E011': {'cell_type': 'hESC_Derived_CD184+', 'tissue': 'ESC_Derived', 'category': 'ESC_Derived'},
    'E012': {'cell_type': 'hESC_Derived_CD56+', 'tissue': 'ESC_Derived', 'category': 'ESC_Derived'},
    'E013': {'cell_type': 'hESC_Derived_CD326+', 'tissue': 'ESC_Derived', 'category': 'ESC_Derived'},
    'E014': {'cell_type': 'HUES48', 'tissue': 'ESC', 'category': 'ESC'},
    'E015': {'cell_type': 'HUES6', 'tissue': 'ESC', 'category': 'ESC'},
    'E016': {'cell_type': 'HUES64', 'tissue': 'ESC', 'category': 'ESC'},
    'E017': {'cell_type': 'IMR90', 'tissue': 'Lung', 'category': 'Other'},
    'E018': {'cell_type': 'iPS-15b', 'tissue': 'iPSC', 'category': 'iPSC'},
    'E019': {'cell_type': 'iPS-18', 'tissue': 'iPSC', 'category': 'iPSC'},
    'E020': {'cell_type': 'iPS-20b', 'tissue': 'iPSC', 'category': 'iPSC'},
    'E021': {'cell_type': 'iPS_DF_6.9', 'tissue': 'iPSC', 'category': 'iPSC'},
    'E022': {'cell_type': 'iPS_DF_19.11', 'tissue': 'iPSC', 'category': 'iPSC'},
    'E023': {'cell_type': 'Mesenchymal_Stem_Cell', 'tissue': 'Mesenchymal', 'category': 'Stem_Cell'},
    'E024': {'cell_type': 'ES-UCSF4', 'tissue': 'ESC', 'category': 'ESC'},
    'E025': {'cell_type': 'Adipose_Derived_MSC', 'tissue': 'Adipose', 'category': 'Stem_Cell'}
}

# Histone modification marks and their biological significance
HISTONE_MARKS = {
    'H3K4me3': {
        'name': 'H3K4me3',
        'description': 'Trimethylation of lysine 4 on histone H3',
        'function': 'Active promoter mark',
        'location': 'Promoters of active genes',
        'enzyme': 'SET1/MLL complexes'
    },
    'H3K4me1': {
        'name': 'H3K4me1',
        'description': 'Monomethylation of lysine 4 on histone H3',
        'function': 'Enhancer mark',
        'location': 'Active and poised enhancers',
        'enzyme': 'SET1/MLL complexes'
    },
    'H3K27ac': {
        'name': 'H3K27ac',
        'description': 'Acetylation of lysine 27 on histone H3',
        'function': 'Active enhancer/promoter mark',
        'location': 'Active regulatory elements',
        'enzyme': 'p300/CBP'
    },
    'H3K27me3': {
        'name': 'H3K27me3',
        'description': 'Trimethylation of lysine 27 on histone H3',
        'function': 'Repressive mark',
        'location': 'Polycomb-repressed regions',
        'enzyme': 'EZH2 (PRC2 complex)'
    },
    'H3K9me3': {
        'name': 'H3K9me3',
        'description': 'Trimethylation of lysine 9 on histone H3',
        'function': 'Heterochromatin mark',
        'location': 'Constitutive heterochromatin',
        'enzyme': 'SUV39H1/H2, SETDB1'
    },
    'H3K36me3': {
        'name': 'H3K36me3',
        'description': 'Trimethylation of lysine 36 on histone H3',
        'function': 'Gene body mark',
        'location': 'Actively transcribed gene bodies',
        'enzyme': 'SETD2'
    }
}

# Chromatin states and their characteristics
CHROMATIN_STATES = {
    'Active_Promoter': {
        'name': 'Active Promoter',
        'description': 'Actively transcribed promoter regions',
        'histone_signature': ['H3K4me3', 'H3K27ac'],
        'color': '#ff0000',
        'methylation': 'Low'
    },
    'Weak_Promoter': {
        'name': 'Weak Promoter',
        'description': 'Weakly active promoter regions',
        'histone_signature': ['H3K4me3'],
        'color': '#ff6969',
        'methylation': 'Low-Medium'
    },
    'Poised_Promoter': {
        'name': 'Poised Promoter',
        'description': 'Bivalent promoter with both active and repressive marks',
        'histone_signature': ['H3K4me3', 'H3K27me3'],
        'color': '#cf0cf0',
        'methylation': 'Variable'
    },
    'Strong_Enhancer': {
        'name': 'Strong Enhancer',
        'description': 'Strongly active enhancer regions',
        'histone_signature': ['H3K4me1', 'H3K27ac'],
        'color': '#ffa500',
        'methylation': 'Low'
    },
    'Weak_Enhancer': {
        'name': 'Weak Enhancer',
        'description': 'Weakly active enhancer regions',
        'histone_signature': ['H3K4me1'],
        'color': '#ffff00',
        'methylation': 'Low-Medium'
    },
    'Polycomb_Repressed': {
        'name': 'Polycomb Repressed',
        'description': 'Polycomb-mediated repression',
        'histone_signature': ['H3K27me3'],
        'color': '#c0c0c0',
        'methylation': 'Variable'
    },
    'Heterochromatin': {
        'name': 'Heterochromatin',
        'description': 'Constitutive heterochromatin',
        'histone_signature': ['H3K9me3'],
        'color': '#969696',
        'methylation': 'High'
    },
    'Quiescent': {
        'name': 'Quiescent/Low Signal',
        'description': 'Regions with low histone modification signal',
        'histone_signature': [],
        'color': '#ffffff',
        'methylation': 'Variable'
    }
}

# Key genes for epigenetic analysis
DEMO_GENES = {
    'pluripotency': {
        'POU5F1': {
            'name': 'POU5F1 (OCT4)',
            'chr': 'chr6',
            'start': 31128153,
            'end': 31142413,
            'description': 'Master pluripotency transcription factor',
            'function': 'Maintains embryonic stem cell pluripotency'
        },
        'NANOG': {
            'name': 'NANOG',
            'chr': 'chr12',
            'start': 7940000,
            'end': 7950000,
            'description': 'Core pluripotency factor',
            'function': 'Prevents differentiation of embryonic stem cells'
        },
        'SOX2': {
            'name': 'SOX2',
            'chr': 'chr3',
            'start': 181711925,
            'end': 181714436,
            'description': 'SRY-box transcription factor 2',
            'function': 'Maintains pluripotency and neural development'
        }
    },
    'neural': {
        'PAX6': {
            'name': 'PAX6',
            'chr': 'chr11',
            'start': 31784717,
            'end': 31818243,
            'description': 'Paired box protein 6',
            'function': 'Master regulator of eye and neural development'
        },
        'SOX1': {
            'name': 'SOX1',
            'chr': 'chr13',
            'start': 112721000,
            'end': 112724000,
            'description': 'SRY-box transcription factor 1',
            'function': 'Neural tube development'
        },
        'FOXG1': {
            'name': 'FOXG1',
            'chr': 'chr14',
            'start': 29236000,
            'end': 29239000,
            'description': 'Forkhead box protein G1',
            'function': 'Forebrain development'
        }
    },
    'cancer': {
        'CDKN2A': {
            'name': 'CDKN2A (p16)',
            'chr': 'chr9',
            'start': 21967751,
            'end': 21995300,
            'description': 'Cyclin-dependent kinase inhibitor 2A',
            'function': 'Tumor suppressor, cell cycle regulation'
        },
        'MLH1': {
            'name': 'MLH1',
            'chr': 'chr3',
            'start': 37034840,
            'end': 37107177,
            'description': 'MutL homolog 1',
            'function': 'DNA mismatch repair, tumor suppressor'
        },
        'BRCA1': {
            'name': 'BRCA1',
            'chr': 'chr17',
            'start': 41196312,
            'end': 41277500,
            'description': 'Breast cancer 1',
            'function': 'DNA repair, tumor suppressor'
        }
    }
}

# Sample methylation data structure
def generate_sample_methylation_data(n_regions=1000, n_samples=25):
    """
    Generate sample methylation data for demonstration.
    
    Args:
        n_regions (int): Number of genomic regions
        n_samples (int): Number of samples
        
    Returns:
        pd.DataFrame: Sample methylation data
    """
    np.random.seed(42)  # For reproducible results
    
    # Generate genomic coordinates
    chromosomes = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
    
    data = []
    for i in range(n_regions):
        chr_name = np.random.choice(chromosomes)
        start = np.random.randint(1000000, 200000000)
        end = start + np.random.randint(200, 2000)
        
        region_data = {
            'chr': chr_name,
            'start': start,
            'end': end,
            'region_id': f"{chr_name}:{start}-{end}"
        }
        
        # Generate methylation values for each sample
        # Different patterns for different region types
        region_type = np.random.choice(['promoter', 'enhancer', 'gene_body', 'intergenic'], 
                                     p=[0.1, 0.15, 0.25, 0.5])
        
        for j, sample_id in enumerate(list(SAMPLE_METADATA.keys())[:n_samples]):
            if region_type == 'promoter':
                # Promoters typically have low methylation
                methylation = np.random.beta(1, 4)
            elif region_type == 'enhancer':
                # Enhancers have variable methylation
                methylation = np.random.beta(2, 2)
            elif region_type == 'gene_body':
                # Gene bodies often have higher methylation
                methylation = np.random.beta(3, 2)
            else:  # intergenic
                # Intergenic regions have variable methylation
                methylation = np.random.beta(2, 3)
            
            # Add some noise and sample-specific effects
            sample_category = SAMPLE_METADATA[sample_id]['category']
            if sample_category == 'ESC':
                # ESCs tend to have lower methylation at regulatory regions
                if region_type in ['promoter', 'enhancer']:
                    methylation *= 0.7
            
            # Add noise
            methylation += np.random.normal(0, 0.05)
            methylation = np.clip(methylation, 0, 1)
            
            region_data[sample_id] = methylation
        
        data.append(region_data)
    
    return pd.DataFrame(data)

# Sample histone modification data
def generate_sample_histone_data(n_regions=500, window_size=100):
    """
    Generate sample histone modification data around TSS.
    
    Args:
        n_regions (int): Number of gene regions
        window_size (int): Number of bins around TSS
        
    Returns:
        pd.DataFrame: Sample histone modification data
    """
    np.random.seed(42)
    
    data = []
    for i in range(n_regions):
        gene_name = f"Gene_{i+1}"
        
        # Generate histone modification profiles
        # Create realistic patterns around TSS
        positions = np.arange(-window_size//2, window_size//2)
        
        # H3K4me3 - sharp peak at TSS
        h3k4me3 = np.exp(-0.5 * (positions/10)**2) * np.random.uniform(0.5, 3.0)
        h3k4me3 += np.random.normal(0, 0.1, len(positions))
        
        # H3K4me1 - broader peak, slightly downstream
        h3k4me1 = np.exp(-0.5 * ((positions-5)/20)**2) * np.random.uniform(0.3, 2.0)
        h3k4me1 += np.random.normal(0, 0.1, len(positions))
        
        # H3K27ac - similar to H3K4me1 but more variable
        h3k27ac = np.exp(-0.5 * ((positions-3)/15)**2) * np.random.uniform(0.2, 2.5)
        h3k27ac += np.random.normal(0, 0.15, len(positions))
        
        # H3K27me3 - anti-correlated with active marks
        h3k27me3 = np.random.uniform(0, 1) * (1 - h3k4me3/np.max(h3k4me3))
        h3k27me3 += np.random.normal(0, 0.1, len(positions))
        
        # H3K9me3 - generally low at active promoters
        h3k9me3 = np.random.uniform(0, 0.5, len(positions))
        h3k9me3 += np.random.normal(0, 0.05, len(positions))
        
        # Ensure non-negative values
        h3k4me3 = np.maximum(h3k4me3, 0)
        h3k4me1 = np.maximum(h3k4me1, 0)
        h3k27ac = np.maximum(h3k27ac, 0)
        h3k27me3 = np.maximum(h3k27me3, 0)
        h3k9me3 = np.maximum(h3k9me3, 0)
        
        for j, pos in enumerate(positions):
            data.append({
                'gene': gene_name,
                'position': pos,
                'H3K4me3': h3k4me3[j],
                'H3K4me1': h3k4me1[j],
                'H3K27ac': h3k27ac[j],
                'H3K27me3': h3k27me3[j],
                'H3K9me3': h3k9me3[j]
            })
    
    return pd.DataFrame(data)

# Analysis parameters
ANALYSIS_PARAMS = {
    'methylation_thresholds': {
        'unmethylated': 0.2,
        'partially_methylated': 0.8
    },
    'histone_thresholds': {
        'H3K4me3': 2.0,
        'H3K4me1': 1.5,
        'H3K27ac': 1.5,
        'H3K27me3': 1.0,
        'H3K9me3': 1.0
    },
    'differential_analysis': {
        'min_methylation_diff': 0.2,
        'p_value_threshold': 0.05,
        'min_samples_per_group': 3
    },
    'visualization': {
        'default_figsize': (10, 6),
        'heatmap_figsize': (12, 8),
        'color_palette': 'RdYlBu_r',
        'dpi': 300
    }
}

# File paths and URLs
DATA_SOURCES = {
    'roadmap_portal': 'https://egg2.wustl.edu/roadmap/web_portal/index.html',
    'roadmap_paper': 'https://doi.org/10.1038/nature14248',
    'ucsc_browser': 'https://genome.ucsc.edu/',
    'encode_portal': 'https://www.encodeproject.org/'
}

def get_gene_info(gene_name, category=None):
    """
    Get information about a specific gene.
    
    Args:
        gene_name (str): Name of the gene
        category (str): Gene category (pluripotency, neural, cancer)
        
    Returns:
        dict: Gene information
    """
    if category:
        return DEMO_GENES.get(category, {}).get(gene_name, None)
    else:
        # Search across all categories
        for cat_genes in DEMO_GENES.values():
            if gene_name in cat_genes:
                return cat_genes[gene_name]
        return None

def get_histone_mark_info(mark_name):
    """
    Get information about a histone modification mark.
    
    Args:
        mark_name (str): Name of the histone mark
        
    Returns:
        dict: Histone mark information
    """
    return HISTONE_MARKS.get(mark_name, None)

def get_chromatin_state_info(state_name):
    """
    Get information about a chromatin state.
    
    Args:
        state_name (str): Name of the chromatin state
        
    Returns:
        dict: Chromatin state information
    """
    return CHROMATIN_STATES.get(state_name, None) 