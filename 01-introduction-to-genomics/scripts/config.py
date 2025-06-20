"""
Configuration file for Genomic Data Analysis
============================================

This module contains configuration settings and constants for the
genomic data analysis module.

Author: Emediong Francis
Part of: Stanford Data Ocean - Data Science in Precision Medicine
"""

import os
from typing import Dict, List, Optional

# Database table mappings for AWS Athena
DATABASE_TABLES = {
    # 1000 Genomes Project datasets
    'g1000_csv': 'default.g1000vcf_csv',
    'g1000_csv_int': 'default.g1000vcf_csv_int',
    'g1000_parquet': 'default.g1000vcf_parquet', 
    'g1000_partitioned': 'default.g1000vcf_partitioned',
    
    # Annotation databases
    'cosmic68': '1000_genomes.hg19_cosmic68_int',
    'refgene': '1000_genomes.hg19_ucsc_refgene_int',
}

# Population code mappings for VCF INFO field
POPULATION_CODES = {
    'EAS_AF': 'East Asian',
    'AMR_AF': 'Ad Mixed American',
    'EUR_AF': 'European', 
    'AFR_AF': 'African',
    'SAS_AF': 'South Asian',
    'AF': 'Global',
}

# Well-known SNPs for demonstration
DEMO_SNPS = {
    'rs12913832': {
        'name': 'Eye Color (Blue/Brown)',
        'chromosome': '15',
        'phenotype': 'Eye pigmentation',
        'description': 'Associated with blue vs brown eye color, shows strong population differences'
    },
    'rs53576': {
        'name': 'Empathy and Social Behavior',
        'chromosome': '3', 
        'phenotype': 'Behavioral traits',
        'description': 'Associated with empathy, social behavior, and stress response'
    },
    'rs72921001': {
        'name': 'Cilantro Taste Preference',
        'chromosome': '11',
        'phenotype': 'Taste perception', 
        'description': 'Associated with cilantro taste preference (love it or hate it)'
    },
    'rs28936679': {
        'name': 'Sleep Patterns',
        'chromosome': '19',
        'phenotype': 'Circadian rhythm',
        'description': 'Associated with sleep duration and circadian rhythm disorders'
    },
    'rs1805009': {
        'name': 'Skin Cancer Risk',
        'chromosome': '16',
        'phenotype': 'Cancer susceptibility',
        'description': 'Associated with melanoma and skin cancer risk'
    }
}

# Gene symbols for interval-based analysis
TUMOR_SUPPRESSOR_GENES = {
    'TP53': {
        'chromosome': '17',
        'name': 'Tumor Protein p53',
        'description': 'Guardian of the Genome - controls cell cycle and apoptosis',
        'clinical_significance': 'Mutations found in >50% of human cancers'
    },
    'BRCA1': {
        'chromosome': '17', 
        'name': 'Breast Cancer gene 1',
        'description': 'DNA repair gene associated with hereditary breast and ovarian cancer',
        'clinical_significance': 'High-penetrance cancer predisposition gene'
    },
    'BRCA2': {
        'chromosome': '13',
        'name': 'Breast Cancer gene 2', 
        'description': 'DNA repair gene associated with hereditary breast and ovarian cancer',
        'clinical_significance': 'High-penetrance cancer predisposition gene'
    },
    'APC': {
        'chromosome': '5',
        'name': 'Adenomatous Polyposis Coli',
        'description': 'Tumor suppressor gene involved in colorectal cancer',
        'clinical_significance': 'Associated with familial adenomatous polyposis'
    }
}

# Analysis parameters
ANALYSIS_PARAMS = {
    'default_sample_size': 10,
    'max_query_limit': 10000,
    'quality_threshold': 20.0,
    'frequency_threshold': 0.01,  # 1%
    'visualization_dpi': 300,
    'figure_width': 12,
    'figure_height': 8,
}

# Color schemes for visualizations
COLOR_SCHEMES = {
    'populations': {
        'East Asian': '#FF6B6B',
        'Ad Mixed American': '#4ECDC4', 
        'European': '#45B7D1',
        'African': '#96CEB4',
        'South Asian': '#FFEAA7',
        'Global': '#DDA0DD'
    },
    'chromosomes': [
        '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
        '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe',
        '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000',
        '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080'
    ]
}

def get_database_table(dataset_key: str) -> str:
    """
    Get fully qualified table name for dataset.
    
    Args:
        dataset_key: Key for the dataset (e.g., 'g1000_csv_int')
        
    Returns:
        Fully qualified table name
        
    Raises:
        ValueError: If dataset key is not found
    """
    if dataset_key not in DATABASE_TABLES:
        available_keys = list(DATABASE_TABLES.keys())
        raise ValueError(f"Unknown dataset '{dataset_key}'. Available: {available_keys}")
    
    return DATABASE_TABLES[dataset_key]

def get_demo_snp_info(rsid: str) -> Dict:
    """
    Get information about a demonstration SNP.
    
    Args:
        rsid: SNP identifier (e.g., 'rs12913832')
        
    Returns:
        Dictionary with SNP information
    """
    return DEMO_SNPS.get(rsid, {
        'name': 'Unknown SNP',
        'chromosome': 'Unknown',
        'phenotype': 'Unknown',
        'description': 'No description available'
    })

def validate_chromosome(chromosome: str) -> bool:
    """
    Validate chromosome identifier.
    
    Args:
        chromosome: Chromosome identifier
        
    Returns:
        True if valid, False otherwise
    """
    valid_chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']
    return chromosome in valid_chromosomes

def validate_rsid(rsid: str) -> bool:
    """
    Validate SNP identifier format.
    
    Args:
        rsid: SNP identifier
        
    Returns:
        True if valid format, False otherwise
    """
    return rsid.startswith('rs') and rsid[2:].isdigit()

# AWS configuration
def get_aws_config() -> Dict[str, str]:
    """
    Get AWS configuration from environment variables.
    
    Returns:
        Dictionary with AWS configuration
    """
    return {
        'region': os.getenv('AWS_DEFAULT_REGION', 'us-east-1'),
        's3_staging_dir': os.getenv('AWS_S3_STAGING_DIR'),
        'profile_name': os.getenv('AWS_PROFILE')
    } 