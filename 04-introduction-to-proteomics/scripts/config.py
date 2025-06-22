"""
Configuration file for proteomic data analysis.

Contains constants, sample data, and configuration parameters
for the Introduction to Proteomics module.
"""

import pandas as pd
import numpy as np

# Sample metadata for proteomics study (iPOP-inspired design)
SAMPLE_METADATA = {
    'Patient_Z_Baseline': {'condition': 'Healthy', 'timepoint': 'Baseline', 'patient': 'Patient_Z', 'infection_status': 'None'},
    'Patient_Z_Infection_Early': {'condition': 'Infection', 'timepoint': 'Early', 'patient': 'Patient_Z', 'infection_status': 'Acute'},
    'Patient_Z_Infection_Peak': {'condition': 'Infection', 'timepoint': 'Peak', 'patient': 'Patient_Z', 'infection_status': 'Acute'},
    'Patient_Z_Recovery_1': {'condition': 'Recovery', 'timepoint': 'Recovery_1', 'patient': 'Patient_Z', 'infection_status': 'Recovering'},
    'Patient_Z_Recovery_2': {'condition': 'Recovery', 'timepoint': 'Recovery_2', 'patient': 'Patient_Z', 'infection_status': 'Recovering'},
    'Patient_Z_Recovered': {'condition': 'Healthy', 'timepoint': 'Recovered', 'patient': 'Patient_Z', 'infection_status': 'None'},
    'Control_1': {'condition': 'Control', 'timepoint': 'Single', 'patient': 'Control_1', 'infection_status': 'None'},
    'Control_2': {'condition': 'Control', 'timepoint': 'Single', 'patient': 'Control_2', 'infection_status': 'None'},
    'Control_3': {'condition': 'Control', 'timepoint': 'Single', 'patient': 'Control_3', 'infection_status': 'None'},
    'Control_4': {'condition': 'Control', 'timepoint': 'Single', 'patient': 'Control_4', 'infection_status': 'None'},
    'Patient_A_Healthy': {'condition': 'Healthy', 'timepoint': 'Single', 'patient': 'Patient_A', 'infection_status': 'None'},
    'Patient_B_Healthy': {'condition': 'Healthy', 'timepoint': 'Single', 'patient': 'Patient_B', 'infection_status': 'None'}
}

# Protein pathway definitions for enrichment analysis
PROTEIN_PATHWAYS = {
    'Acute Phase Response': {
        'description': 'Proteins involved in acute phase inflammatory response',
        'proteins': ['CRP', 'SAA1', 'SAA2', 'HP', 'HPX', 'FGA', 'FGB', 'FGG', 
                    'APOA1', 'APOA2', 'TTR', 'ALB', 'TF', 'CP', 'ORM1', 'ORM2'],
        'source': 'KEGG'
    },
    'Complement System': {
        'description': 'Complement cascade proteins and regulators',
        'proteins': ['C1QA', 'C1QB', 'C1QC', 'C1R', 'C1S', 'C2', 'C3', 'C4A', 'C4B', 
                    'C5', 'C6', 'C7', 'C8A', 'C8B', 'C9', 'CFB', 'CFD', 'CFH', 'CFI'],
        'source': 'KEGG'
    },
    'Blood Coagulation': {
        'description': 'Coagulation cascade and fibrinolysis proteins',
        'proteins': ['F2', 'F5', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12', 'F13A1', 'F13B',
                    'FGA', 'FGB', 'FGG', 'PLG', 'SERPINC1', 'SERPIND1', 'SERPINE1', 'PROC'],
        'source': 'KEGG'
    },
    'Immunoglobulins': {
        'description': 'Antibody heavy and light chains',
        'proteins': ['IGHG1', 'IGHG2', 'IGHG3', 'IGHG4', 'IGHA1', 'IGHA2', 'IGHM', 'IGHD',
                    'IGKC', 'IGKV1', 'IGKV2', 'IGKV3', 'IGLC1', 'IGLC2', 'IGLV1', 'IGLV2'],
        'source': 'GO'
    },
    'Lipid Transport': {
        'description': 'Lipoproteins and lipid transport proteins',
        'proteins': ['APOA1', 'APOA2', 'APOA4', 'APOB', 'APOC1', 'APOC2', 'APOC3', 'APOE',
                    'LCAT', 'CETP', 'PLTP', 'LPL', 'LDLR', 'ABCA1', 'ABCG1', 'SCARB1'],
        'source': 'KEGG'
    },
    'Protein Folding': {
        'description': 'Molecular chaperones and protein folding machinery',
        'proteins': ['HSP90AA1', 'HSP90AB1', 'HSP90B1', 'HSPA1A', 'HSPA1B', 'HSPA5', 'HSPA8',
                    'HSPB1', 'HSPE1', 'HSPD1', 'CALR', 'CANX', 'PDIA3', 'PDIA4', 'PDIA6', 'ERP29'],
        'source': 'GO'
    },
    'Proteolysis': {
        'description': 'Proteases and protein degradation',
        'proteins': ['CTSD', 'CTSB', 'CTSL', 'CTSK', 'CTSZ', 'CTSS', 'CTSH', 'CTSF',
                    'MMP1', 'MMP2', 'MMP3', 'MMP9', 'TIMP1', 'TIMP2', 'A2M', 'SERPINA1'],
        'source': 'GO'
    },
    'Oxidative Stress': {
        'description': 'Antioxidant enzymes and oxidative stress response',
        'proteins': ['SOD1', 'SOD2', 'SOD3', 'CAT', 'GPX1', 'GPX3', 'PRDX1', 'PRDX2', 'PRDX6',
                    'TXN', 'TXNRD1', 'GSR', 'GCLC', 'HMOX1', 'NQO1', 'GSTP1'],
        'source': 'GO'
    },
    'Cell Adhesion': {
        'description': 'Cell adhesion molecules and extracellular matrix proteins',
        'proteins': ['ITGA1', 'ITGA2', 'ITGB1', 'ITGB3', 'CDH1', 'CDH2', 'VCAM1', 'ICAM1',
                    'PECAM1', 'NCAM1', 'FN1', 'VTN', 'THBS1', 'COMP', 'TNC', 'LAMB1'],
        'source': 'GO'
    },
    'Metabolic Enzymes': {
        'description': 'Key metabolic pathway enzymes',
        'proteins': ['GAPDH', 'PKM', 'LDHA', 'LDHB', 'ENO1', 'ENO2', 'PGK1', 'ALDOA',
                    'TPI1', 'PFKL', 'PFKM', 'HK1', 'HK2', 'CS', 'IDH1', 'MDH2'],
        'source': 'KEGG'
    }
}

# Demo proteins with detailed information
DEMO_PROTEINS = {
    'acute_phase': {
        'CRP': {
            'name': 'C-Reactive Protein',
            'description': 'Major acute phase protein, inflammatory marker',
            'molecular_weight': 25.0,
            'expected_response': {'Infection': 'up', 'Recovery': 'down'}
        },
        'SAA1': {
            'name': 'Serum Amyloid A1',
            'description': 'Acute phase protein, dramatically increases during inflammation',
            'molecular_weight': 13.5,
            'expected_response': {'Infection': 'up', 'Recovery': 'down'}
        },
        'HP': {
            'name': 'Haptoglobin',
            'description': 'Acute phase protein, binds free hemoglobin',
            'molecular_weight': 45.2,
            'expected_response': {'Infection': 'up', 'Recovery': 'down'}
        },
        'ALB': {
            'name': 'Albumin',
            'description': 'Most abundant plasma protein, negative acute phase',
            'molecular_weight': 69.4,
            'expected_response': {'Infection': 'down', 'Recovery': 'up'}
        }
    },
    'complement': {
        'C3': {
            'name': 'Complement Component 3',
            'description': 'Central component of complement system',
            'molecular_weight': 187.0,
            'expected_response': {'Infection': 'variable', 'Recovery': 'stable'}
        },
        'CFD': {
            'name': 'Complement Factor D',
            'description': 'Alternative pathway complement factor',
            'molecular_weight': 27.0,
            'expected_response': {'Infection': 'down', 'Recovery': 'up'}
        },
        'C4A': {
            'name': 'Complement Component 4A',
            'description': 'Classical pathway complement component',
            'molecular_weight': 192.8,
            'expected_response': {'Infection': 'variable', 'Recovery': 'stable'}
        }
    },
    'coagulation': {
        'FGA': {
            'name': 'Fibrinogen Alpha Chain',
            'description': 'Fibrinogen component, acute phase protein',
            'molecular_weight': 94.9,
            'expected_response': {'Infection': 'up', 'Recovery': 'down'}
        },
        'FGB': {
            'name': 'Fibrinogen Beta Chain',
            'description': 'Fibrinogen component, coagulation factor',
            'molecular_weight': 55.9,
            'expected_response': {'Infection': 'up', 'Recovery': 'down'}
        },
        'PLG': {
            'name': 'Plasminogen',
            'description': 'Fibrinolysis protein, converts to plasmin',
            'molecular_weight': 90.6,
            'expected_response': {'Infection': 'variable', 'Recovery': 'stable'}
        }
    },
    'lipoproteins': {
        'APOA1': {
            'name': 'Apolipoprotein A1',
            'description': 'Major HDL apolipoprotein, negative acute phase',
            'molecular_weight': 30.8,
            'expected_response': {'Infection': 'down', 'Recovery': 'up'}
        },
        'APOB': {
            'name': 'Apolipoprotein B',
            'description': 'Major LDL apolipoprotein',
            'molecular_weight': 515.6,
            'expected_response': {'Infection': 'down', 'Recovery': 'stable'}
        },
        'APOE': {
            'name': 'Apolipoprotein E',
            'description': 'Lipoprotein component, lipid transport',
            'molecular_weight': 36.2,
            'expected_response': {'Infection': 'variable', 'Recovery': 'stable'}
        }
    },
    'immunoglobulins': {
        'IGHG1': {
            'name': 'Immunoglobulin Heavy Chain Gamma 1',
            'description': 'IgG1 heavy chain, most abundant antibody',
            'molecular_weight': 51.3,
            'expected_response': {'Infection': 'up', 'Recovery': 'stable'}
        },
        'IGKC': {
            'name': 'Immunoglobulin Kappa Constant',
            'description': 'Kappa light chain constant region',
            'molecular_weight': 11.8,
            'expected_response': {'Infection': 'up', 'Recovery': 'stable'}
        }
    }
}

# iPOP study information
IPOP_STUDY_INFO = {
    'name': 'Integrated Personal Omics Profiling',
    'institution': 'Stanford University',
    'description': 'Longitudinal multi-omics study of health and disease',
    'patient_count': 106,
    'focus': 'Personal precision medicine through omics profiling',
    'timepoints': 'Regular intervals over several years',
    'conditions': ['Healthy', 'Pre-diabetic', 'Diabetic', 'Infection', 'Recovery']
}

# Mass spectrometry parameters
MS_PARAMETERS = {
    'instrument_type': 'Orbitrap',
    'acquisition_mode': 'DDA',  # Data Dependent Acquisition
    'mass_range': '350-1600 m/z',
    'resolution': 70000,
    'dynamic_exclusion': '30 seconds',
    'normalization_method': 'Total Ion Current',
    'quantification_method': 'Label-free'
}

# Analysis parameters
ANALYSIS_PARAMS = {
    'differential_expression': {
        'p_value_threshold': 0.05,
        'fold_change_threshold': 1.5,
        'effect_size_threshold': 0.5,
        'min_detection_frequency': 0.5,
        'multiple_testing_method': 'fdr_bh'
    },
    'pathway_enrichment': {
        'p_value_threshold': 0.05,
        'min_pathway_size': 3,
        'max_pathway_size': 500,
        'background_correction': True
    },
    'quality_control': {
        'max_missing_percentage': 50,
        'max_cv': 1.0,
        'min_abundance': 0.1,
        'outlier_threshold': 3
    },
    'normalization': {
        'default_method': 'median_centering',
        'missing_value_method': 'knn',
        'log_transform': True
    }
}

# File format specifications
FILE_FORMATS = {
    'abundance_matrix': {
        'format': 'TSV',
        'separator': '\t',
        'protein_id_column': 0,
        'sample_start_column': 1,
        'expected_columns': ['Protein.IDs', 'Gene.names', 'sample_columns...']
    },
    'sample_metadata': {
        'format': 'TSV',
        'separator': '\t',
        'sample_id_column': 0,
        'required_columns': ['condition', 'timepoint'],
        'optional_columns': ['patient', 'infection_status', 'batch']
    }
}

def generate_sample_abundance_matrix(n_proteins=5000, n_samples=12, seed=42):
    """
    Generate sample protein abundance matrix for demonstration.
    
    Args:
        n_proteins (int): Number of proteins to simulate
        n_samples (int): Number of samples to simulate
        seed (int): Random seed for reproducibility
        
    Returns:
        pd.DataFrame: Simulated abundance matrix
    """
    np.random.seed(seed)
    
    # Generate protein names
    protein_names = [f"PROTEIN_{i+1:05d}" for i in range(n_proteins)]
    
    # Add real protein names for demonstration
    real_proteins = []
    for category in DEMO_PROTEINS.values():
        real_proteins.extend(list(category.keys()))
    
    # Replace first proteins with real protein names
    for i, protein in enumerate(real_proteins[:min(len(real_proteins), n_proteins)]):
        protein_names[i] = protein
    
    # Generate sample names
    sample_names = list(SAMPLE_METADATA.keys())[:n_samples]
    
    # Generate abundance data with realistic patterns
    abundance_data = []
    
    for i, protein in enumerate(protein_names):
        protein_abundances = []
        
        # Base abundance level (log-normal distribution, pre-log transformed)
        base_abundance = np.random.normal(mean=20, sigma=3)
        
        for sample in sample_names:
            condition = SAMPLE_METADATA[sample]['condition']
            timepoint = SAMPLE_METADATA[sample]['timepoint']
            
            # Start with base abundance
            abundance = base_abundance
            
            # Add condition effects for known proteins
            if protein in real_proteins:
                # Get expected response for this protein and condition
                for category in DEMO_PROTEINS.values():
                    if protein in category:
                        expected = category[protein].get('expected_response', {})
                        if condition in expected:
                            response = expected[condition]
                            if response == 'up':
                                abundance += np.random.uniform(2, 5)
                            elif response == 'down':
                                abundance -= np.random.uniform(2, 5)
                            elif response == 'variable':
                                abundance += np.random.uniform(-1, 1)
                        break
            
            # Add infection timeline effects for Patient Z
            if 'Patient_Z' in sample:
                if 'Infection' in timepoint:
                    # Acute phase proteins go up, others may go down
                    if protein in ['CRP', 'SAA1', 'HP', 'FGA', 'FGB']:
                        abundance += np.random.uniform(3, 6)
                    elif protein in ['ALB', 'APOA1', 'APOB']:
                        abundance -= np.random.uniform(1, 3)
                elif 'Recovery' in timepoint:
                    # Gradual return to baseline
                    if protein in ['CRP', 'SAA1', 'HP']:
                        abundance += np.random.uniform(0, 2)  # Still elevated but decreasing
                    elif 'CFD' in protein:
                        abundance += np.random.uniform(1, 3)  # Recovery marker
            
            # Add biological variation
            bio_variation = np.random.normal(0, 1)
            abundance += bio_variation
            
            # Add technical noise
            tech_noise = np.random.normal(0, 0.5)
            abundance += tech_noise
            
            # Add some missing values (common in proteomics)
            if np.random.random() < 0.05:  # 5% missing values
                abundance = np.nan
            
            protein_abundances.append(abundance)
        
        abundance_data.append(protein_abundances)
    
    # Create DataFrame
    abundance_matrix = pd.DataFrame(abundance_data, 
                                  index=protein_names, 
                                  columns=sample_names)
    
    return abundance_matrix

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

def get_protein_pathways():
    """
    Get protein pathways in format suitable for analysis.
    
    Returns:
        dict: Pathway name -> protein list mapping
    """
    pathways = {}
    for pathway, info in PROTEIN_PATHWAYS.items():
        pathways[pathway] = info['proteins']
    
    return pathways

def get_analysis_contrasts():
    """
    Get standard analysis contrasts for proteomics study.
    
    Returns:
        list: List of contrast tuples (condition1, condition2)
    """
    contrasts = [
        ('Infection', 'Healthy'),
        ('Recovery', 'Healthy'),
        ('Recovery', 'Infection')
    ]
    
    return contrasts

def get_demo_protein_info(protein_name):
    """
    Get information about a demo protein.
    
    Args:
        protein_name (str): Protein name
        
    Returns:
        dict: Protein information or None if not found
    """
    for category, proteins in DEMO_PROTEINS.items():
        if protein_name in proteins:
            info = proteins[protein_name].copy()
            info['category'] = category
            return info
    
    return None

def validate_abundance_matrix(abundance_matrix):
    """
    Validate abundance matrix format and content.
    
    Args:
        abundance_matrix (pd.DataFrame): Abundance matrix to validate
        
    Returns:
        dict: Validation results
    """
    results = {
        'valid': True,
        'warnings': [],
        'errors': []
    }
    
    # Check for infinite values
    if np.isinf(abundance_matrix.select_dtypes(include=[np.number])).any().any():
        results['errors'].append("Abundance matrix contains infinite values")
        results['valid'] = False
    
    # Check missing value percentage
    missing_pct = abundance_matrix.isnull().sum().sum() / abundance_matrix.size
    if missing_pct > 0.5:
        results['warnings'].append(f"High missing value percentage: {missing_pct*100:.1f}%")
    
    # Check for all-missing proteins
    all_missing_proteins = abundance_matrix.isnull().all(axis=1).sum()
    if all_missing_proteins > 0:
        results['errors'].append(f"{all_missing_proteins} proteins have no quantified values")
        results['valid'] = False
    
    # Check abundance range (should be reasonable for log-transformed data)
    numeric_data = abundance_matrix.select_dtypes(include=[np.number])
    if not numeric_data.empty:
        min_val = numeric_data.min().min()
        max_val = numeric_data.max().max()
        
        if min_val < 0 or max_val > 50:
            results['warnings'].append(f"Unusual abundance range: {min_val:.2f} to {max_val:.2f}")
    
    return results

def create_design_matrix(metadata, formula="~ condition"):
    """
    Create design matrix for statistical analysis.
    
    Args:
        metadata (pd.DataFrame): Sample metadata
        formula (str): Design formula
        
    Returns:
        pd.DataFrame: Design matrix
    """
    # Simple design matrix creation
    if formula == "~ condition":
        # One-hot encode condition
        design = pd.get_dummies(metadata['condition'], prefix='condition')
        design.index = metadata.index
        return design
    
    elif formula == "~ condition + patient":
        # Include patient as blocking factor
        condition_design = pd.get_dummies(metadata['condition'], prefix='condition')
        patient_design = pd.get_dummies(metadata['patient'], prefix='patient')
        
        design = pd.concat([condition_design, patient_design], axis=1)
        design.index = metadata.index
        return design
    
    else:
        raise ValueError(f"Unsupported formula: {formula}")

# Data source URLs and references
DATA_SOURCES = {
    'iPOP': 'https://med.stanford.edu/ipop.html',
    'PRIDE': 'https://www.ebi.ac.uk/pride/',
    'ProteomeXchange': 'http://www.proteomexchange.org/',
    'Human_Proteome_Map': 'http://www.humanproteomemap.org/',
    'MaxQuant': 'https://www.maxquant.org/'
}

# Software and tool references
TOOL_REFERENCES = {
    'MaxQuant': 'Cox, J., Mann, M. (2008). MaxQuant enables high peptide identification rates. Nature Biotechnology.',
    'Perseus': 'Tyanova, S., et al. (2016). The Perseus computational platform for comprehensive analysis of (prote)omics data. Nature Methods.',
    'MSstats': 'Choi, M., et al. (2014). MSstats: an R package for statistical analysis of quantitative mass spectrometry-based proteomic experiments. Bioinformatics.',
    'limma': 'Ritchie, M.E., et al. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research.',
    'STRING': 'Szklarczyk, D., et al. (2019). STRING v11: protein-protein association networks with increased coverage. Nucleic Acids Research.',
    'UniProt': 'The UniProt Consortium (2019). UniProt: a worldwide hub of protein knowledge. Nucleic Acids Research.'
} 