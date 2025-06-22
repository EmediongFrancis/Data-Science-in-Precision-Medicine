#!/usr/bin/env python3
"""
Clinical Data Analysis Demo Script
Stanford Data Science in Precision Medicine - Module 5

This comprehensive demo showcases clinical laboratory data analysis using
the iPOP study framework, focusing on diabetes progression patterns and
precision medicine applications.

Author: Stanford Data Science in Precision Medicine Team
"""

import sys
import os
sys.path.append('./scripts')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# Import our clinical analysis utilities
from clinical_utils import (
    ClinicalDataProcessor,
    ClinicalBiomarkerAnalyzer,
    ClinicalDataVisualizer,
    ClinicalPredictiveModeling
)
from config import DEMO_CONFIG, get_diabetes_thresholds

def main():
    """
    Main demo function showcasing comprehensive clinical data analysis.
    """
    print("=" * 80)
    print("STANFORD DATA SCIENCE IN PRECISION MEDICINE")
    print("Module 5: Introduction to Clinical Data Analysis")
    print("=" * 80)
    print()
    
    print("🏥 iPOP Study: Clinical Laboratory Analysis Demo")
    print("📊 Analyzing longitudinal biomarker data for diabetes progression")
    print("🎯 Focus: Precision medicine through clinical data interpretation")
    print()
    
    # Initialize analysis components
    print("🔧 Initializing Clinical Analysis Components...")
    processor = ClinicalDataProcessor()
    analyzer = ClinicalBiomarkerAnalyzer()
    visualizer = ClinicalDataVisualizer()
    ml_model = ClinicalPredictiveModeling()
    print("✓ Components initialized successfully")
    print()
    
    # Generate comprehensive demo data
    print("📋 Generating iPOP Study Demo Data...")
    print("   Based on real patient cases from Stanford iPOP study:")
    for patient in DEMO_CONFIG['patients']:
        print(f"   • {patient['id']}: {patient['description']}")
    print()
    
    clinical_data = processor.generate_demo_data(DEMO_CONFIG)
    print(f"✓ Generated clinical data: {clinical_data.shape[0]} samples, {clinical_data.shape[1]} features")
    print()
    
    # Display data overview
    print("📊 Clinical Data Overview:")
    print("=" * 50)
    print(clinical_data.head(10))
    print()
    print("📈 Data Summary Statistics:")
    print(clinical_data.describe())
    print()
    
    # Quality control analysis
    print("🔍 Performing Quality Control Analysis...")
    qc_results = processor.quality_control(clinical_data)
    
    print("Quality Control Results:")
    print(f"   • Total samples: {qc_results['total_samples']}")
    print(f"   • Complete cases: {qc_results['complete_cases']}")
    print(f"   • Completeness rate: {qc_results['completeness_rate']:.1f}%")
    
    if qc_results['missing_values']:
        print("   • Missing values detected:")
        for col, percent in qc_results['missing_values'].items():
            print(f"     - {col}: {percent:.1f}%")
    
    if qc_results['outliers']:
        print("   • Outliers detected:")
        for col, count in qc_results['outliers'].items():
            print(f"     - {col}: {count} outliers")
    print()
    
    # Diabetes progression analysis
    print("🩺 Analyzing Diabetes Progression Patterns...")
    progression_results = analyzer.diabetes_progression_analysis(clinical_data)
    
    print("Diabetes Progression Analysis:")
    for patient, progression in progression_results.items():
        print(f"   Patient {patient}:")
        if progression['first_diabetic_marker']:
            print(f"     • First diabetic marker: {progression['first_diabetic_marker']} (Day {progression['first_diabetic_day']})")
            print(f"     • Total diabetic markers: {progression['total_diabetic_markers']}")
            print(f"     • Progression events: {len(progression['progression_pattern'])}")
        else:
            print(f"     • No diabetic markers detected")
    print()
    
    # Biomarker correlation analysis
    print("🔗 Analyzing Biomarker Correlations...")
    correlation_results = analyzer.biomarker_correlation_analysis(clinical_data)
    
    print("Significant Biomarker Correlations:")
    for corr in correlation_results['significant_correlations']:
        print(f"   • {corr['marker1']} ↔ {corr['marker2']}: r={corr['correlation']:.3f} ({corr['strength']})")
    print()
    
    # Risk stratification
    print("⚠️  Performing Risk Stratification...")
    risk_results = analyzer.risk_stratification(clinical_data)
    
    print("Patient Risk Stratification:")
    for patient, risk_info in risk_results.items():
        print(f"   Patient {patient}: {risk_info['risk_category']} (Score: {risk_info['risk_score']})")
        if risk_info['risk_factors']:
            print(f"     Risk factors: {', '.join(risk_info['risk_factors'])}")
    print()
    
    # Machine learning analysis
    print("🤖 Training Diabetes Prediction Model...")
    feature_data = ml_model.prepare_features(clinical_data)
    ml_results = ml_model.train_diabetes_classifier(feature_data)
    
    print("Machine Learning Results:")
    print(f"   • Model accuracy: {ml_results['test_accuracy']:.3f}")
    print(f"   • ROC AUC: {ml_results['roc_auc']:.3f}")
    print("   • Feature importance:")
    for feature, importance in ml_results['feature_importance'].items():
        print(f"     - {feature}: {importance:.3f}")
    print()
    
    # Comprehensive visualizations
    print("📊 Creating Comprehensive Visualizations...")
    
    # Create figure directory if it doesn't exist
    os.makedirs('figures', exist_ok=True)
    
    # 1. Individual patient longitudinal plots
    print("   Creating longitudinal biomarker plots...")
    for i, patient_info in enumerate(DEMO_CONFIG['patients'][:3]):  # First 3 patients
        patient_id = patient_info['id']
        fig, ax = visualizer.plot_longitudinal_biomarkers(
            clinical_data, 
            patient_id,
            title=f"Patient {patient_id}: {patient_info['description']}"
        )
        plt.savefig(f'figures/patient_{patient_id}_longitudinal.png', dpi=300, bbox_inches='tight')
        print(f"     ✓ Saved: patient_{patient_id}_longitudinal.png")
        plt.show()
    
    # 2. Diabetes progression heatmap
    print("   Creating diabetes progression heatmap...")
    fig, ax = visualizer.plot_diabetes_progression_heatmap(progression_results)
    plt.savefig('figures/diabetes_progression_heatmap.png', dpi=300, bbox_inches='tight')
    print("     ✓ Saved: diabetes_progression_heatmap.png")
    plt.show()
    
    # 3. Biomarker correlation heatmap
    print("   Creating biomarker correlation heatmap...")
    fig, ax = visualizer.plot_biomarker_correlations(correlation_results['correlation_matrix'])
    plt.savefig('figures/biomarker_correlations.png', dpi=300, bbox_inches='tight')
    print("     ✓ Saved: biomarker_correlations.png")
    plt.show()
    
    # 4. Risk stratification plots
    print("   Creating risk stratification visualization...")
    fig, axes = visualizer.plot_risk_stratification(risk_results)
    plt.savefig('figures/risk_stratification.png', dpi=300, bbox_inches='tight')
    print("     ✓ Saved: risk_stratification.png")
    plt.show()
    
    print()
    
    # Clinical insights and interpretation
    print("🏥 Clinical Insights and Interpretation:")
    print("=" * 50)
    
    print("\n1. DIABETES PROGRESSION PATHWAYS:")
    print("   The analysis reveals three distinct pathways to diabetes:")
    print("   • OGTT-First Progression: Glucose tolerance impairment precedes fasting abnormalities")
    print("   • FPG-First Progression: Fasting glucose homeostasis disruption occurs first")
    print("   • Infection-Triggered: Environmental factors (viral infections) trigger metabolic dysfunction")
    
    print("\n2. BIOMARKER RELATIONSHIPS:")
    print("   Strong correlations between clinical markers indicate:")
    print("   • FPG and HbA1C: Reflect chronic glycemic control")
    print("   • OGTT and SSPG: Both assess insulin resistance pathways")
    print("   • Temporal progression: Different markers become abnormal at different stages")
    
    print("\n3. PRECISION MEDICINE APPLICATIONS:")
    print("   Clinical data analysis enables:")
    print("   • Early Detection: Identify pre-diabetic states before clinical diagnosis")
    print("   • Risk Stratification: Personalized risk assessment based on biomarker profiles")
    print("   • Treatment Monitoring: Track intervention effectiveness through biomarker changes")
    print("   • Pathway-Specific Interventions: Target specific metabolic pathways based on progression pattern")
    
    print("\n4. MACHINE LEARNING INSIGHTS:")
    print(f"   • Model successfully predicts diabetes with {ml_results['test_accuracy']:.1%} accuracy")
    print("   • Feature importance reveals most predictive biomarkers:")
    top_features = sorted(ml_results['feature_importance'].items(), key=lambda x: x[1], reverse=True)
    for feature, importance in top_features[:2]:
        print(f"     - {feature}: {importance:.1%} contribution")
    
    print("\n5. LONGITUDINAL PATTERNS:")
    print("   Time-series analysis shows:")
    print("   • Individual Variability: Each patient follows unique progression trajectory")
    print("   • Critical Windows: Certain timepoints show rapid biomarker changes")
    print("   • Recovery Potential: Some patients show improvement with lifestyle interventions")
    print("   • Relapse Risk: Environmental triggers can cause disease progression")
    
    print("\n" + "=" * 80)
    print("SUMMARY: CLINICAL DATA IN PRECISION MEDICINE")
    print("=" * 80)
    
    print("\nThis analysis demonstrates how clinical laboratory data serves as the foundation")
    print("for precision medicine approaches:")
    print()
    print("🎯 PERSONALIZED RISK ASSESSMENT:")
    print("   Individual biomarker profiles enable tailored risk stratification")
    print()
    print("📊 PATHWAY-SPECIFIC INTERVENTIONS:")
    print("   Understanding progression patterns guides targeted therapeutic approaches")
    print()
    print("🔬 PREDICTIVE MODELING:")
    print("   Machine learning models can predict disease outcomes and treatment responses")
    print()
    print("📈 LONGITUDINAL MONITORING:")
    print("   Continuous biomarker tracking enables real-time health optimization")
    print()
    print("🏥 CLINICAL DECISION SUPPORT:")
    print("   Data-driven insights enhance clinical decision-making and patient care")
    
    print("\n" + "=" * 80)
    print("Demo completed successfully! 🎉")
    print("All visualizations saved in the 'figures/' directory")
    print("=" * 80)

if __name__ == "__main__":
    # Create directories if they don't exist
    os.makedirs('figures', exist_ok=True)
    os.makedirs('data', exist_ok=True)
    
    # Run the comprehensive demo
    main()
