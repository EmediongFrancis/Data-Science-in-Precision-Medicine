"""
Configuration file for transcriptomic data analysis.

Contains constants, sample data, and configuration parameters
for the Introduction to Transcriptomics module.
"""

import pandas as pd
import numpy as np

# Sample metadata for hormone treatment study (based on original R module)
SAMPLE_METADATA = {
    'Control_1': {'condition': 'Control', 'patient': 'Patient_1', 'batch': 'Batch_1'},
    'Control_2': {'condition': 'Control', 'patient': 'Patient_2', 'batch': 'Batch_1'},
    'Control_3': {'condition': 'Control', 'patient': 'Patient_3', 'batch': 'Batch_2'},
    'Control_4': {'condition': 'Control', 'patient': 'Patient_4', 'batch': 'Batch_2'},
    'OHT_1': {'condition': 'OHT', 'patient': 'Patient_1', 'batch': 'Batch_1'},
    'OHT_2': {'condition': 'OHT', 'patient': 'Patient_2', 'batch': 'Batch_1'},
    'OHT_3': {'condition': 'OHT', 'patient': 'Patient_3', 'batch': 'Batch_2'},
    'OHT_4': {'condition': 'OHT', 'patient': 'Patient_4', 'batch': 'Batch_2'},
    'DPN_1': {'condition': 'DPN', 'patient': 'Patient_1', 'batch': 'Batch_1'},
    'DPN_2': {'condition': 'DPN', 'patient': 'Patient_2', 'batch': 'Batch_1'},
    'DPN_3': {'condition': 'DPN', 'patient': 'Patient_3', 'batch': 'Batch_2'},
    'DPN_4': {'condition': 'DPN', 'patient': 'Patient_4', 'batch': 'Batch_2'}
}

# Gene sets for pathway analysis (curated from major biological pathways)
PATHWAY_GENE_SETS = {
    'Cell Cycle': {
        'description': 'Genes involved in cell cycle regulation and progression',
        'genes': ['CCNA2', 'CCNB1', 'CCND1', 'CCNE1', 'CDK1', 'CDK2', 'CDK4', 'CDK6', 
                 'CDKN1A', 'CDKN1B', 'CDKN2A', 'PCNA', 'MKI67', 'TOP2A', 'BUB1', 'AURKA'],
        'source': 'KEGG'
    },
    'Apoptosis': {
        'description': 'Programmed cell death pathway genes',
        'genes': ['TP53', 'BAX', 'BCL2', 'BCL2L1', 'CASP3', 'CASP8', 'CASP9', 'FAS', 
                 'FASLG', 'APAF1', 'CYCS', 'BID', 'BAK1', 'MCL1', 'PUMA', 'NOXA'],
        'source': 'KEGG'
    },
    'Immune Response': {
        'description': 'Immune system activation and inflammatory response',
        'genes': ['IL1B', 'TNF', 'IFNG', 'IL6', 'IL10', 'CD8A', 'CD4', 'CXCL10', 
                 'CCL2', 'TLR4', 'NFKB1', 'STAT1', 'IRF1', 'HLA-DRA', 'CD68', 'CD86'],
        'source': 'GO'
    },
    'Metabolism': {
        'description': 'Metabolic pathways and energy production',
        'genes': ['GAPDH', 'PKM', 'LDHA', 'HK1', 'HK2', 'PFKL', 'ALDOA', 'ENO1', 
                 'PGK1', 'TPI1', 'PDHB', 'CS', 'IDH1', 'OGDH', 'SUCLA2', 'MDH2'],
        'source': 'KEGG'
    },
    'DNA Repair': {
        'description': 'DNA damage response and repair mechanisms',
        'genes': ['BRCA1', 'BRCA2', 'ATM', 'ATR', 'TP53BP1', 'RAD51', 'RAD52', 'XRCC1', 
                 'PARP1', 'MSH2', 'MLH1', 'ERCC1', 'XPA', 'PCNA', 'RFC1', 'LIG1'],
        'source': 'GO'
    },
    'Transcription Regulation': {
        'description': 'Transcriptional control and gene expression regulation',
        'genes': ['MYC', 'JUN', 'FOS', 'ETS1', 'SP1', 'STAT1', 'STAT3', 'NF1', 
                 'E2F1', 'RB1', 'CREB1', 'ATF1', 'ELK1', 'SRF', 'TCF7L2', 'LEF1'],
        'source': 'GO'
    },
    'Signal Transduction': {
        'description': 'Cellular signaling pathways and communication',
        'genes': ['EGFR', 'MAPK1', 'MAPK3', 'AKT1', 'PIK3CA', 'MTOR', 'KRAS', 'RAF1', 
                 'MAP2K1', 'PTEN', 'GSK3B', 'CTNNB1', 'WNT3A', 'SMAD2', 'TGFB1', 'NOTCH1'],
        'source': 'KEGG'
    },
    'Hormone Response': {
        'description': 'Hormone signaling and endocrine response',
        'genes': ['ESR1', 'ESR2', 'PGR', 'AR', 'NR3C1', 'PPARA', 'PPARG', 'RXRA', 
                 'VDR', 'RARA', 'THRB', 'NCOR1', 'CREBBP', 'EP300', 'GREB1', 'TFF1'],
        'source': 'Custom'
    },
    'Oxidative Stress': {
        'description': 'Cellular response to oxidative stress and ROS',
        'genes': ['SOD1', 'SOD2', 'CAT', 'GPX1', 'PRDX1', 'TXNRD1', 'NRF2', 'KEAP1', 
                 'HMOX1', 'NQO1', 'GCLC', 'GCLM', 'GSR', 'GSTP1', 'MT1A', 'MT2A'],
        'source': 'GO'
    },
    'Angiogenesis': {
        'description': 'Blood vessel formation and vascular development',
        'genes': ['VEGFA', 'VEGFR2', 'ANGPT1', 'ANGPT2', 'TIE1', 'PDGFA', 'PDGFB', 'FGF2', 
                 'FGFR1', 'HIF1A', 'EPAS1', 'PECAM1', 'CDH5', 'TEK', 'NOS3', 'THBS1'],
        'source': 'GO'
    }
}

# Demo genes for specific analyses
DEMO_GENES = {
    'hormone_responsive': {
        'ESR1': {
            'name': 'Estrogen Receptor 1',
            'description': 'Primary estrogen receptor mediating estrogen signaling',
            'chromosome': 'chr6',
            'expected_response': {'OHT': 'up', 'DPN': 'up'}
        },
        'PGR': {
            'name': 'Progesterone Receptor',
            'description': 'Nuclear receptor for progesterone',
            'chromosome': 'chr11',
            'expected_response': {'OHT': 'up', 'DPN': 'down'}
        },
        'TFF1': {
            'name': 'Trefoil Factor 1',
            'description': 'Estrogen-responsive gene, breast cancer marker',
            'chromosome': 'chr21',
            'expected_response': {'OHT': 'up', 'DPN': 'up'}
        },
        'GREB1': {
            'name': 'Growth Regulation by Estrogen in Breast Cancer 1',
            'description': 'Estrogen-regulated gene in breast tissue',
            'chromosome': 'chr2',
            'expected_response': {'OHT': 'up', 'DPN': 'up'}
        }
    },
    'cell_cycle': {
        'MKI67': {
            'name': 'Marker of Proliferation Ki-67',
            'description': 'Cell proliferation marker, expressed in cycling cells',
            'chromosome': 'chr10',
            'expected_response': {'OHT': 'up', 'DPN': 'variable'}
        },
        'CCND1': {
            'name': 'Cyclin D1',
            'description': 'G1/S-specific cyclin, cell cycle progression',
            'chromosome': 'chr11',
            'expected_response': {'OHT': 'up', 'DPN': 'up'}
        },
        'CDKN1A': {
            'name': 'Cyclin Dependent Kinase Inhibitor 1A (p21)',
            'description': 'Cell cycle inhibitor, p53 target',
            'chromosome': 'chr6',
            'expected_response': {'OHT': 'down', 'DPN': 'variable'}
        }
    },
    'metabolism': {
        'GAPDH': {
            'name': 'Glyceraldehyde-3-Phosphate Dehydrogenase',
            'description': 'Key glycolytic enzyme, commonly used housekeeping gene',
            'chromosome': 'chr12',
            'expected_response': {'OHT': 'stable', 'DPN': 'stable'}
        },
        'PKM': {
            'name': 'Pyruvate Kinase M1/2',
            'description': 'Final step of glycolysis, metabolic regulation',
            'chromosome': 'chr15',
            'expected_response': {'OHT': 'up', 'DPN': 'up'}
        }
    },
    'stress_response': {
        'TP53': {
            'name': 'Tumor Protein p53',
            'description': 'Guardian of the genome, DNA damage response',
            'chromosome': 'chr17',
            'expected_response': {'OHT': 'variable', 'DPN': 'variable'}
        },
        'HSP90AA1': {
            'name': 'Heat Shock Protein 90 Alpha Family Class A Member 1',
            'description': 'Molecular chaperone, stress response',
            'chromosome': 'chr14',
            'expected_response': {'OHT': 'up', 'DPN': 'up'}
        }
    }
}

# Experimental design parameters
EXPERIMENTAL_DESIGN = {
    'study_design': 'paired',
    'factors': ['condition', 'patient'],
    'conditions': ['Control', 'OHT', 'DPN'],
    'n_patients': 4,
    'n_replicates': 1,
    'blocking_factor': 'patient'
}

# Treatment information
TREATMENT_INFO = {
    'Control': {
        'full_name': 'Vehicle Control',
        'description': 'Ethanol vehicle control treatment',
        'duration': '24 hours',
        'concentration': 'Vehicle only'
    },
    'OHT': {
        'full_name': '4-Hydroxytamoxifen',
        'description': 'Active metabolite of tamoxifen, selective estrogen receptor modulator',
        'duration': '24 hours',
        'concentration': '100 nM',
        'mechanism': 'Estrogen receptor antagonist/agonist'
    },
    'DPN': {
        'full_name': 'Diarylpropionitrile',
        'description': 'Selective estrogen receptor beta agonist',
        'duration': '24 hours',
        'concentration': '10 nM',
        'mechanism': 'ERβ-selective agonist'
    }
}

# Analysis parameters
ANALYSIS_PARAMS = {
    'differential_expression': {
        'p_value_threshold': 0.05,
        'fold_change_threshold': 1.5,
        'min_expression_level': 10,
        'min_samples_expressed': 3,
        'multiple_testing_method': 'fdr_bh'
    },
    'pathway_enrichment': {
        'p_value_threshold': 0.05,
        'min_gene_set_size': 5,
        'max_gene_set_size': 500,
        'background_correction': True
    },
    'quality_control': {
        'min_library_size': 1000000,
        'max_pct_top100': 50,
        'min_detected_genes': 5000,
        'outlier_threshold': 3
    },
    'normalization': {
        'default_method': 'log_cpm',
        'variance_stabilization': True,
        'batch_correction': False
    }
}

# File format specifications
FILE_FORMATS = {
    'count_matrix': {
        'format': 'TSV',
        'separator': '\t',
        'gene_id_column': 0,
        'sample_start_column': 1,
        'expected_columns': ['gene_id', 'transcript_id.s.', 'sample_columns...']
    },
    'sample_metadata': {
        'format': 'TSV',
        'separator': '\t',
        'sample_id_column': 0,
        'required_columns': ['condition', 'patient'],
        'optional_columns': ['batch', 'sequencing_date', 'library_prep']
    }
}

def generate_sample_count_matrix(n_genes=15000, n_samples=12, seed=42):
    """
    Generate sample RNA-seq count matrix for demonstration.
    
    Args:
        n_genes (int): Number of genes to simulate
        n_samples (int): Number of samples to simulate
        seed (int): Random seed for reproducibility
        
    Returns:
        pd.DataFrame: Simulated count matrix
    """
    np.random.seed(seed)
    
    # Generate gene names
    gene_names = [f"GENE_{i+1:05d}" for i in range(n_genes)]
    
    # Add some real gene names for demonstration
    real_genes = []
    for category in DEMO_GENES.values():
        real_genes.extend(list(category.keys()))
    
    # Replace first genes with real gene names
    for i, gene in enumerate(real_genes[:min(len(real_genes), n_genes)]):
        gene_names[i] = gene
    
    # Generate sample names
    sample_names = list(SAMPLE_METADATA.keys())[:n_samples]
    
    # Generate count data with realistic patterns
    count_data = []
    
    for i, gene in enumerate(gene_names):
        gene_counts = []
        
        # Base expression level (log-normal distribution)
        base_expression = np.random.lognormal(mean=5, sigma=2)
        
        for sample in sample_names:
            condition = SAMPLE_METADATA[sample]['condition']
            patient = SAMPLE_METADATA[sample]['patient']
            
            # Start with base expression
            expression = base_expression
            
            # Add condition effects for known genes
            if gene in real_genes:
                # Get expected response for this gene and condition
                for category in DEMO_GENES.values():
                    if gene in category:
                        expected = category[gene].get('expected_response', {})
                        if condition in expected:
                            response = expected[condition]
                            if response == 'up':
                                expression *= np.random.uniform(1.5, 3.0)
                            elif response == 'down':
                                expression *= np.random.uniform(0.3, 0.7)
                            elif response == 'variable':
                                expression *= np.random.uniform(0.8, 1.2)
                        break
            
            # Add patient effects (biological variation)
            patient_effect = np.random.normal(1.0, 0.2)
            expression *= patient_effect
            
            # Add technical noise
            technical_noise = np.random.normal(1.0, 0.1)
            expression *= technical_noise
            
            # Convert to counts (Poisson-like)
            count = np.random.poisson(expression)
            count = max(0, count)  # Ensure non-negative
            
            gene_counts.append(count)
        
        count_data.append(gene_counts)
    
    # Create DataFrame
    count_matrix = pd.DataFrame(count_data, 
                               index=gene_names, 
                               columns=sample_names)
    
    return count_matrix

def get_sample_metadata_df():
    """
    Convert sample metadata to DataFrame format.
    
    Returns:
        pd.DataFrame: Sample metadata
    """
    metadata_df = pd.DataFrame.from_dict(SAMPLE_METADATA, orient='index')
    
    # Convert to categorical
    for col in metadata_df.columns:
        metadata_df[col] = metadata_df[col].astype('category')
    
    return metadata_df

def get_pathway_gene_sets():
    """
    Get pathway gene sets in format suitable for analysis.
    
    Returns:
        dict: Pathway name -> gene list mapping
    """
    gene_sets = {}
    for pathway, info in PATHWAY_GENE_SETS.items():
        gene_sets[pathway] = info['genes']
    
    return gene_sets

def get_treatment_contrasts():
    """
    Get standard treatment contrasts for analysis.
    
    Returns:
        list: List of contrast tuples (treatment, control)
    """
    contrasts = [
        ('OHT', 'Control'),
        ('DPN', 'Control'),
        ('OHT', 'DPN')
    ]
    
    return contrasts

def get_demo_gene_info(gene_name):
    """
    Get information about a demo gene.
    
    Args:
        gene_name (str): Gene name
        
    Returns:
        dict: Gene information or None if not found
    """
    for category, genes in DEMO_GENES.items():
        if gene_name in genes:
            info = genes[gene_name].copy()
            info['category'] = category
            return info
    
    return None

def validate_count_matrix(count_matrix):
    """
    Validate count matrix format and content.
    
    Args:
        count_matrix (pd.DataFrame): Count matrix to validate
        
    Returns:
        dict: Validation results
    """
    results = {
        'valid': True,
        'warnings': [],
        'errors': []
    }
    
    # Check for negative values
    if (count_matrix < 0).any().any():
        results['errors'].append("Count matrix contains negative values")
        results['valid'] = False
    
    # Check for non-integer values
    if not count_matrix.dtypes.apply(lambda x: np.issubdtype(x, np.integer)).all():
        results['warnings'].append("Count matrix contains non-integer values")
    
    # Check for all-zero genes
    zero_genes = (count_matrix.sum(axis=1) == 0).sum()
    if zero_genes > 0:
        results['warnings'].append(f"{zero_genes} genes have zero counts across all samples")
    
    # Check library sizes
    library_sizes = count_matrix.sum(axis=0)
    min_lib_size = ANALYSIS_PARAMS['quality_control']['min_library_size']
    low_lib_samples = (library_sizes < min_lib_size).sum()
    if low_lib_samples > 0:
        results['warnings'].append(f"{low_lib_samples} samples have library size < {min_lib_size}")
    
    return results

# Data source URLs and references
DATA_SOURCES = {
    'GEO': 'https://www.ncbi.nlm.nih.gov/geo/',
    'TCGA': 'https://portal.gdc.cancer.gov/',
    'GTEx': 'https://gtexportal.org/',
    'Human_Cell_Atlas': 'https://www.humancellatlas.org/',
    'ENCODE': 'https://www.encodeproject.org/'
}

# Software and tool references
TOOL_REFERENCES = {
    'DESeq2': 'Love, M.I., Huber, W., Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology.',
    'edgeR': 'Robinson, M.D., McCarthy, D.J., Smyth, G.K. (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics.',
    'limma': 'Ritchie, M.E., et al. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research.',
    'GSEA': 'Subramanian, A., et al. (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. PNAS.',
    'STAR': 'Dobin, A., et al. (2013). STAR: ultrafast universal RNA-seq aligner. Bioinformatics.',
    'RSEM': 'Li, B., Dewey, C.N. (2011). RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. BMC Bioinformatics.'
} 