"""
Configuration file for Clinical Data Analysis Module
Stanford Data Science in Precision Medicine - Module 5

Contains demo data configurations, reference ranges, and analysis parameters
based on the iPOP study and clinical laboratory standards.
"""

import numpy as np

# Demo configuration based on iPOP study Patient Z and similar cases
DEMO_CONFIG = {
    "study_info": {
        "name": "iPOP Study - Clinical Laboratory Analysis",
        "description": "Longitudinal clinical biomarker analysis from Stanford iPOP study",
        "reference": "Schüssler-Fiorenza Rose et al., Nature Medicine 2019",
        "focus": "Diabetes progression through different pathways",
        "timeframe": "Up to 8 years longitudinal follow-up"
    },
    
    "patients": [
        {
            "id": "ZNDMXI3",
            "pattern": "ogtt_progression",
            "description": "Initial abnormality of OGTT2HR - Oral Glucose Tolerance Test",
            "clinical_notes": "First developed OGTT in diabetic range (>200 mg/dl)"
        },
        {
            "id": "ZNED4XZ", 
            "pattern": "fpg_progression",
            "description": "Initial abnormality of FPG - Fasting Plasma Glucose",
            "clinical_notes": "First had abnormality in fasting glucose homeostasis (>126 mg/dl)"
        },
        {
            "id": "ZOZOWIT",
            "pattern": "viral_infection",
            "description": "Initial improvement followed by progression after viral infections",
            "clinical_notes": "Diabetic range after viral infection, improved with lifestyle, worsened after second infection"
        },
        {
            "id": "HEALTHY01",
            "pattern": "normal_control",
            "description": "Healthy control participant",
            "clinical_notes": "Maintained normal glucose homeostasis throughout study"
        },
        {
            "id": "PREDIAB01",
            "pattern": "prediabetic_stable",
            "description": "Stable prediabetic participant",
            "clinical_notes": "Consistent prediabetic markers without progression to diabetes"
        }
    ],
    
    "timepoints": [0, 90, 180, 270, 360, 450, 540, 630, 720, 810],  # Days in study
    
    "clinical_markers": {
        "FPG": {
            "name": "Fasting Plasma Glucose",
            "unit": "mg/dl",
            "description": "Measures glucose homeostasis in fasting state",
            "reference_ranges": {
                "normal": "70-99 mg/dl",
                "prediabetic": "100-125 mg/dl", 
                "diabetic": "≥126 mg/dl"
            },
            "clinical_significance": "Primary marker for diabetes diagnosis"
        },
        "OGTT2HR": {
            "name": "Oral Glucose Tolerance Test (2-hour)",
            "unit": "mg/dl", 
            "description": "Measures glucose processing after 75g glucose load",
            "reference_ranges": {
                "normal": "<140 mg/dl",
                "prediabetic": "140-199 mg/dl",
                "diabetic": "≥200 mg/dl"
            },
            "clinical_significance": "Assesses glucose tolerance and insulin response"
        },
        "SSPG": {
            "name": "Steady-State Plasma Glucose",
            "unit": "mg/dl",
            "description": "Measures insulin resistance",
            "reference_ranges": {
                "normal": "<150 mg/dl",
                "insulin_resistant": "≥150 mg/dl"
            },
            "clinical_significance": "Direct measure of insulin sensitivity"
        },
        "HbA1C": {
            "name": "Hemoglobin A1C (Glycated Hemoglobin)",
            "unit": "%",
            "description": "Average blood glucose over 2-3 months",
            "reference_ranges": {
                "normal": "4.0-5.6%",
                "prediabetic": "5.7-6.4%",
                "diabetic": "≥6.5%"
            },
            "clinical_significance": "Long-term glycemic control marker"
        }
    },
    
    "diabetes_pathways": {
        "glucose_homeostasis": {
            "description": "Fasting glucose regulation pathway",
            "primary_marker": "FPG",
            "associated_processes": ["hepatic glucose production", "insulin sensitivity", "pancreatic beta-cell function"]
        },
        "glucose_tolerance": {
            "description": "Post-prandial glucose processing",
            "primary_marker": "OGTT2HR", 
            "associated_processes": ["insulin secretion", "peripheral glucose uptake", "incretin response"]
        },
        "insulin_resistance": {
            "description": "Cellular insulin response pathway",
            "primary_marker": "SSPG",
            "associated_processes": ["muscle glucose uptake", "adipose tissue metabolism", "liver insulin sensitivity"]
        },
        "chronic_glycemia": {
            "description": "Long-term glucose exposure",
            "primary_marker": "HbA1C",
            "associated_processes": ["protein glycation", "oxidative stress", "vascular complications"]
        }
    },
    
    "clinical_scenarios": {
        "scenario_1": {
            "title": "OGTT-First Diabetes Progression",
            "patient": "ZNDMXI3",
            "description": "Patient develops diabetes through glucose tolerance pathway first",
            "learning_objectives": [
                "Understand different pathways to diabetes",
                "Recognize OGTT as early diabetes indicator",
                "Analyze progression patterns"
            ]
        },
        "scenario_2": {
            "title": "FPG-First Diabetes Progression", 
            "patient": "ZNED4XZ",
            "description": "Patient develops diabetes through fasting glucose pathway first",
            "learning_objectives": [
                "Identify fasting glucose abnormalities",
                "Understand glucose homeostasis disruption",
                "Correlate with HbA1C changes"
            ]
        },
        "scenario_3": {
            "title": "Infection-Triggered Diabetes",
            "patient": "ZOZOWIT", 
            "description": "Viral infections trigger diabetes with recovery and relapse pattern",
            "learning_objectives": [
                "Understand environmental triggers",
                "Analyze recovery and relapse patterns",
                "Evaluate lifestyle intervention effects"
            ]
        }
    },
    
    "analysis_parameters": {
        "missing_value_threshold": 0.3,  # Maximum proportion of missing values allowed
        "outlier_detection_method": "IQR",  # Interquartile range method
        "correlation_threshold": 0.5,  # Minimum correlation for significance
        "diabetes_prediction_features": ["FPG", "OGTT2HR", "SSPG", "HbA1C"],
        "longitudinal_trend_window": 180,  # Days for trend analysis
        "risk_stratification_weights": {
            "FPG_diabetic": 3,
            "FPG_prediabetic": 2,
            "OGTT_diabetic": 3,
            "OGTT_prediabetic": 2,
            "HbA1C_diabetic": 3,
            "HbA1C_prediabetic": 2,
            "SSPG_resistant": 2
        }
    },
    
    "visualization_settings": {
        "color_palette": {
            "normal": "#2ECC71",
            "prediabetic": "#F39C12", 
            "diabetic": "#E74C3C",
            "insulin_resistant": "#8E44AD"
        },
        "plot_style": "seaborn-v0_8",
        "figure_size": (12, 8),
        "font_sizes": {
            "title": 16,
            "labels": 14,
            "annotations": 12,
            "diabetic_annotations": 16
        },
        "marker_styles": {
            "FPG": {"color": "#F1C40F", "marker": "o", "size": 8},
            "OGTT2HR": {"color": "#3498DB", "marker": "o", "size": 10},
            "SSPG": {"color": "#E74C3C", "marker": "s", "size": 10},
            "HbA1C": {"color": "lightgray", "marker": "o", "size": 8}
        }
    }
}

# Sample clinical data generation parameters
CLINICAL_DATA_PARAMS = {
    "noise_levels": {
        "FPG": {"mean": 0, "std": 8},
        "OGTT2HR": {"mean": 0, "std": 15},
        "SSPG": {"mean": 0, "std": 12},
        "HbA1C": {"mean": 0, "std": 0.2}
    },
    
    "progression_rates": {
        "ogtt_progression": {
            "FPG_increase_per_timepoint": 2,
            "OGTT_increase_per_timepoint": 15,
            "SSPG_increase_per_timepoint": 5,
            "HbA1C_increase_per_timepoint": 0.2
        },
        "fpg_progression": {
            "FPG_increase_per_timepoint": 8,
            "OGTT_increase_per_timepoint": 5,
            "SSPG_increase_per_timepoint": 3,
            "HbA1C_increase_per_timepoint": 0.15
        }
    },
    
    "baseline_values": {
        "normal": {
            "FPG": {"mean": 85, "std": 8},
            "OGTT2HR": {"mean": 120, "std": 15},
            "SSPG": {"mean": 100, "std": 15},
            "HbA1C": {"mean": 5.2, "std": 0.3}
        },
        "prediabetic": {
            "FPG": {"mean": 110, "std": 10},
            "OGTT2HR": {"mean": 160, "std": 20},
            "SSPG": {"mean": 130, "std": 20},
            "HbA1C": {"mean": 6.0, "std": 0.3}
        }
    }
}

# Quality control thresholds
QC_THRESHOLDS = {
    "FPG": {"min": 50, "max": 400},
    "OGTT2HR": {"min": 50, "max": 500},
    "SSPG": {"min": 50, "max": 300},
    "HbA1C": {"min": 3.0, "max": 15.0}
}

# Machine learning parameters
ML_CONFIG = {
    "features": ["FPG", "OGTT2HR", "SSPG", "HbA1C"],
    "target": "Diabetes",
    "test_size": 0.3,
    "random_state": 42,
    "model_parameters": {
        "n_estimators": 100,
        "max_depth": 10,
        "min_samples_split": 5,
        "min_samples_leaf": 2
    },
    "cross_validation_folds": 5
}

def validate_clinical_data(data):
    """Validate clinical data against quality control thresholds."""
    validation_results = {"valid": True, "warnings": [], "errors": []}
    
    for marker, thresholds in QC_THRESHOLDS.items():
        marker_col = f"{marker} (mg/dl)" if marker != "HbA1C" else f"{marker} (%)"
        
        if marker_col in data.columns:
            # Check for values outside acceptable range
            out_of_range = data[
                (data[marker_col] < thresholds["min"]) | 
                (data[marker_col] > thresholds["max"])
            ]
            
            if len(out_of_range) > 0:
                validation_results["warnings"].append(
                    f"{marker}: {len(out_of_range)} values outside normal range ({thresholds['min']}-{thresholds['max']})"
                )
    
    return validation_results

def get_reference_ranges():
    """Return clinical reference ranges for all markers."""
    return {
        marker: info["reference_ranges"]
        for marker, info in DEMO_CONFIG["clinical_markers"].items()
    }

def get_diabetes_thresholds():
    """Return diabetes diagnostic thresholds."""
    return {
        "FPG": 126,
        "OGTT2HR": 200,
        "HbA1C": 6.5,
        "SSPG": 150  # Insulin resistance threshold
    }
