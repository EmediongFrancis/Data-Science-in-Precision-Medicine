"""
Clinical Data Analysis Utilities for Stanford Data Science in Precision Medicine
Module 5: Introduction to Clinical Data

This module provides comprehensive tools for analyzing clinical laboratory data,
focusing on longitudinal biomarker analysis, diabetes progression tracking,
and precision medicine applications based on the Stanford iPOP study.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, roc_auc_score, roc_curve
import warnings
warnings.filterwarnings('ignore')

class ClinicalDataProcessor:
    """
    Comprehensive clinical data processing and quality control.
    
    Handles longitudinal clinical laboratory data with focus on:
    - Data validation and quality control
    - Missing value analysis and imputation
    - Temporal trend analysis
    - Reference range validation
    - Multi-timepoint data alignment
    """
    
    def __init__(self):
        self.reference_ranges = {
            'FPG': {'normal': (70, 99), 'prediabetic': (100, 125), 'diabetic': 126},
            'OGTT2HR': {'normal': (0, 139), 'prediabetic': (140, 199), 'diabetic': 200},
            'SSPG': {'normal': (0, 149), 'insulin_resistant': 150},
            'HbA1C': {'normal': (4.0, 5.6), 'prediabetic': (5.7, 6.4), 'diabetic': 6.5}
        }
        
    def load_clinical_data(self, file_path):
        """Load and validate clinical data from CSV file."""
        try:
            data = pd.read_csv(file_path)
            print(f"✓ Data loaded successfully: {data.shape}")
            return self.validate_data_structure(data)
        except Exception as e:
            print(f"✗ Error loading data: {e}")
            return None
    
    def validate_data_structure(self, data):
        """Validate clinical data structure and required columns."""
        required_columns = ['Days', 'PatientID']
        clinical_markers = ['FPG', 'OGTT2HR', 'SSPG', 'HbA1C']
        
        missing_cols = [col for col in required_columns if col not in data.columns]
        if missing_cols:
            print(f"⚠️  Missing required columns: {missing_cols}")
        
        available_markers = [col for col in clinical_markers 
                           if any(marker in col for marker in clinical_markers)]
        print(f"✓ Available clinical markers: {len(available_markers)}")
        
        return data
    
    def quality_control(self, data):
        """Comprehensive quality control analysis."""
        qc_results = {}
        
        # Missing value analysis
        missing_percent = (data.isnull().sum() / len(data)) * 100
        qc_results['missing_values'] = missing_percent[missing_percent > 0].to_dict()
        
        # Outlier detection using IQR method
        numeric_cols = data.select_dtypes(include=[np.number]).columns
        outliers = {}
        
        for col in numeric_cols:
            if col != 'Days':
                Q1 = data[col].quantile(0.25)
                Q3 = data[col].quantile(0.75)
                IQR = Q3 - Q1
                lower_bound = Q1 - 1.5 * IQR
                upper_bound = Q3 + 1.5 * IQR
                
                outlier_count = len(data[(data[col] < lower_bound) | (data[col] > upper_bound)])
                if outlier_count > 0:
                    outliers[col] = outlier_count
        
        qc_results['outliers'] = outliers
        
        # Data completeness
        qc_results['total_samples'] = len(data)
        qc_results['complete_cases'] = len(data.dropna())
        qc_results['completeness_rate'] = (qc_results['complete_cases'] / qc_results['total_samples']) * 100
        
        return qc_results
    
    def categorize_values(self, data, marker):
        """Categorize clinical values based on reference ranges."""
        if marker not in self.reference_ranges:
            return data
        
        ranges = self.reference_ranges[marker]
        categories = []
        
        for value in data[f'{marker} (mg/dl)'] if f'{marker} (mg/dl)' in data.columns else data[f'{marker} (%)']:
            if pd.isna(value):
                categories.append('Unknown')
            elif marker in ['FPG', 'OGTT2HR', 'HbA1C']:
                if marker == 'HbA1C':
                    if value <= ranges['normal'][1]:
                        categories.append('Normal')
                    elif value <= ranges['prediabetic'][1]:
                        categories.append('Prediabetic')
                    else:
                        categories.append('Diabetic')
                else:
                    if value <= ranges['normal'][1]:
                        categories.append('Normal')
                    elif value <= ranges['prediabetic'][1]:
                        categories.append('Prediabetic')
                    else:
                        categories.append('Diabetic')
            elif marker == 'SSPG':
                if value < ranges['insulin_resistant']:
                    categories.append('Normal')
                else:
                    categories.append('Insulin Resistant')
        
        return categories
    
    def temporal_analysis(self, data, patient_id=None):
        """Analyze temporal trends in clinical markers."""
        if patient_id:
            patient_data = data[data['PatientID'] == patient_id].copy()
        else:
            patient_data = data.copy()
        
        patient_data = patient_data.sort_values('Days')
        
        trends = {}
        numeric_cols = patient_data.select_dtypes(include=[np.number]).columns
        
        for col in numeric_cols:
            if col != 'Days' and not patient_data[col].isna().all():
                # Calculate trend using linear regression
                valid_data = patient_data[[col, 'Days']].dropna()
                if len(valid_data) > 1:
                    slope, intercept, r_value, p_value, std_err = stats.linregress(
                        valid_data['Days'], valid_data[col]
                    )
                    trends[col] = {
                        'slope': slope,
                        'r_squared': r_value**2,
                        'p_value': p_value,
                        'trend': 'increasing' if slope > 0 else 'decreasing' if slope < 0 else 'stable'
                    }
        
        return trends
    
    def generate_demo_data(self, config):
        """Generate comprehensive demo clinical data based on iPOP study patterns."""
        np.random.seed(42)
        
        patients = config['patients']
        timepoints = config['timepoints']
        
        all_data = []
        
        for patient_info in patients:
            patient_id = patient_info['id']
            pattern = patient_info['pattern']
            
            for i, day in enumerate(timepoints):
                row = {'PatientID': patient_id, 'Days': day}
                
                # Generate realistic clinical values based on patient pattern
                if pattern == 'ogtt_progression':
                    # Patient with OGTT2HR progression (like ZNDMXI3)
                    row['FPG (mg/dl)'] = np.random.normal(95, 10) + (i * 2)
                    row['OGTT2HR (mg/dl)'] = np.random.normal(150, 20) + (i * 15)
                    if i > 3:  # Diabetic range after day 3
                        row['OGTT2HR (mg/dl)'] = max(row['OGTT2HR (mg/dl)'], 200)
                    row['SSPG (mg/dl)'] = np.random.normal(120, 15) + (i * 5)
                    row['HbA1C (%)'] = np.random.normal(5.5, 0.3) + (i * 0.2)
                    
                elif pattern == 'fpg_progression':
                    # Patient with FPG progression (like ZNED4XZ)
                    row['FPG (mg/dl)'] = np.random.normal(110, 15) + (i * 8)
                    if i > 2:  # Diabetic range after day 2
                        row['FPG (mg/dl)'] = max(row['FPG (mg/dl)'], 126)
                    row['OGTT2HR (mg/dl)'] = np.random.normal(160, 25)
                    row['SSPG (mg/dl)'] = np.random.normal(130, 20)
                    row['HbA1C (%)'] = np.random.normal(6.0, 0.4) + (i * 0.15)
                    
                elif pattern == 'viral_infection':
                    # Patient with viral infection pattern (like ZOZOWIT)
                    if i < 2:  # Initial normal
                        row['FPG (mg/dl)'] = np.random.normal(90, 8)
                        row['HbA1C (%)'] = np.random.normal(5.2, 0.2)
                    elif i < 4:  # Post-viral spike
                        row['FPG (mg/dl)'] = np.random.normal(140, 15)
                        row['HbA1C (%)'] = np.random.normal(6.8, 0.3)
                    elif i < 6:  # Recovery period
                        row['FPG (mg/dl)'] = np.random.normal(100, 10)
                        row['HbA1C (%)'] = np.random.normal(5.8, 0.2)
                    else:  # Second infection
                        row['FPG (mg/dl)'] = np.random.normal(150, 20)
                        row['HbA1C (%)'] = np.random.normal(7.2, 0.4)
                        row['OGTT2HR (mg/dl)'] = np.random.normal(220, 30)
                    
                    row['SSPG (mg/dl)'] = np.random.normal(140, 25)
                    if 'OGTT2HR (mg/dl)' not in row:
                        row['OGTT2HR (mg/dl)'] = np.random.normal(170, 30)
                
                # Add some missing values randomly
                for marker in ['FPG (mg/dl)', 'OGTT2HR (mg/dl)', 'SSPG (mg/dl)', 'HbA1C (%)']:
                    if np.random.random() < 0.1:  # 10% missing rate
                        row[marker] = np.nan
                
                all_data.append(row)
        
        return pd.DataFrame(all_data)


class ClinicalBiomarkerAnalyzer:
    """
    Advanced biomarker analysis for clinical data.
    
    Provides statistical analysis tools for:
    - Diabetes progression modeling
    - Biomarker correlation analysis
    - Risk stratification
    - Longitudinal trend analysis
    """
    
    def __init__(self):
        self.diabetes_thresholds = {
            'FPG': 126,
            'OGTT2HR': 200,
            'HbA1C': 6.5,
            'SSPG': 150
        }
    
    def diabetes_progression_analysis(self, data):
        """Analyze diabetes progression patterns across patients."""
        progression_results = {}
        
        patients = data['PatientID'].unique()
        
        for patient in patients:
            patient_data = data[data['PatientID'] == patient].sort_values('Days')
            
            # Check for diabetes progression
            progression = {
                'first_diabetic_marker': None,
                'first_diabetic_day': None,
                'progression_pattern': [],
                'total_diabetic_markers': 0
            }
            
            for _, row in patient_data.iterrows():
                day = row['Days']
                diabetic_markers = []
                
                if not pd.isna(row.get('FPG (mg/dl)', np.nan)) and row['FPG (mg/dl)'] >= self.diabetes_thresholds['FPG']:
                    diabetic_markers.append('FPG')
                if not pd.isna(row.get('OGTT2HR (mg/dl)', np.nan)) and row['OGTT2HR (mg/dl)'] >= self.diabetes_thresholds['OGTT2HR']:
                    diabetic_markers.append('OGTT2HR')
                if not pd.isna(row.get('HbA1C (%)', np.nan)) and row['HbA1C (%)'] >= self.diabetes_thresholds['HbA1C']:
                    diabetic_markers.append('HbA1C')
                if not pd.isna(row.get('SSPG (mg/dl)', np.nan)) and row['SSPG (mg/dl)'] >= self.diabetes_thresholds['SSPG']:
                    diabetic_markers.append('SSPG')
                
                if diabetic_markers and progression['first_diabetic_marker'] is None:
                    progression['first_diabetic_marker'] = diabetic_markers[0]
                    progression['first_diabetic_day'] = day
                
                if diabetic_markers:
                    progression['progression_pattern'].append({
                        'day': day,
                        'markers': diabetic_markers
                    })
            
            progression['total_diabetic_markers'] = len(set([
                marker for pattern in progression['progression_pattern'] 
                for marker in pattern['markers']
            ]))
            
            progression_results[patient] = progression
        
        return progression_results
    
    def biomarker_correlation_analysis(self, data):
        """Analyze correlations between clinical biomarkers."""
        numeric_cols = [col for col in data.columns if 'mg/dl' in col or '%' in col]
        correlation_matrix = data[numeric_cols].corr()
        
        # Find significant correlations
        significant_correlations = []
        for i in range(len(correlation_matrix.columns)):
            for j in range(i+1, len(correlation_matrix.columns)):
                corr_value = correlation_matrix.iloc[i, j]
                if abs(corr_value) > 0.5:  # Threshold for significant correlation
                    significant_correlations.append({
                        'marker1': correlation_matrix.columns[i],
                        'marker2': correlation_matrix.columns[j],
                        'correlation': corr_value,
                        'strength': 'strong' if abs(corr_value) > 0.7 else 'moderate'
                    })
        
        return {
            'correlation_matrix': correlation_matrix,
            'significant_correlations': significant_correlations
        }
    
    def risk_stratification(self, data):
        """Stratify patients based on diabetes risk factors."""
        risk_scores = {}
        
        patients = data['PatientID'].unique()
        
        for patient in patients:
            patient_data = data[data['PatientID'] == patient]
            
            risk_score = 0
            risk_factors = []
            
            # Latest measurements
            latest_data = patient_data.iloc[-1]
            
            # FPG risk assessment
            fpg = latest_data.get('FPG (mg/dl)', np.nan)
            if not pd.isna(fpg):
                if fpg >= 126:
                    risk_score += 3
                    risk_factors.append('Diabetic FPG')
                elif fpg >= 100:
                    risk_score += 2
                    risk_factors.append('Prediabetic FPG')
            
            # HbA1C risk assessment
            hba1c = latest_data.get('HbA1C (%)', np.nan)
            if not pd.isna(hba1c):
                if hba1c >= 6.5:
                    risk_score += 3
                    risk_factors.append('Diabetic HbA1C')
                elif hba1c >= 5.7:
                    risk_score += 2
                    risk_factors.append('Prediabetic HbA1C')
            
            # OGTT risk assessment
            ogtt = latest_data.get('OGTT2HR (mg/dl)', np.nan)
            if not pd.isna(ogtt):
                if ogtt >= 200:
                    risk_score += 3
                    risk_factors.append('Diabetic OGTT')
                elif ogtt >= 140:
                    risk_score += 2
                    risk_factors.append('Prediabetic OGTT')
            
            # SSPG insulin resistance
            sspg = latest_data.get('SSPG (mg/dl)', np.nan)
            if not pd.isna(sspg) and sspg >= 150:
                risk_score += 2
                risk_factors.append('Insulin Resistant')
            
            # Risk categorization
            if risk_score >= 6:
                risk_category = 'High Risk'
            elif risk_score >= 3:
                risk_category = 'Moderate Risk'
            else:
                risk_category = 'Low Risk'
            
            risk_scores[patient] = {
                'risk_score': risk_score,
                'risk_category': risk_category,
                'risk_factors': risk_factors
            }
        
        return risk_scores


class ClinicalDataVisualizer:
    """
    Comprehensive visualization tools for clinical data.
    
    Creates publication-ready plots including:
    - Longitudinal biomarker tracking
    - Diabetes progression visualization
    - Risk stratification plots
    - Correlation heatmaps
    """
    
    def __init__(self):
        plt.style.use('seaborn-v0_8')
        self.colors = {
            'normal': '#2ECC71',
            'prediabetic': '#F39C12',
            'diabetic': '#E74C3C',
            'insulin_resistant': '#8E44AD'
        }
    
    def plot_longitudinal_biomarkers(self, data, patient_id, title=None):
        """Create longitudinal biomarker plot for a specific patient."""
        patient_data = data[data['PatientID'] == patient_id].sort_values('Days')
        
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Plot FPG
        fpg_data = patient_data.dropna(subset=['FPG (mg/dl)'])
        if not fpg_data.empty:
            fpg_line, = ax.plot(fpg_data['Days'], fpg_data['FPG (mg/dl)'], 
                               color='#F1C40F', marker='o', linewidth=2, markersize=8, label='FPG')
            
            # Annotate FPG values
            for _, row in fpg_data.iterrows():
                color = '#E74C3C' if row['FPG (mg/dl)'] >= 126 else 'black'
                fontsize = 14 if row['FPG (mg/dl)'] >= 126 else 12
                ax.annotate(f"{row['FPG (mg/dl)']:.0f}", 
                           (row['Days'], row['FPG (mg/dl)']),
                           textcoords="offset points", xytext=(0,-25), 
                           ha='center', color=color, fontsize=fontsize, weight='bold' if color == '#E74C3C' else 'normal')
        
        # Plot OGTT2HR
        ogtt_data = patient_data.dropna(subset=['OGTT2HR (mg/dl)'])
        if not ogtt_data.empty:
            ogtt_line, = ax.plot(ogtt_data['Days'], ogtt_data['OGTT2HR (mg/dl)'], 
                                color='#3498DB', marker='o', linewidth=2, markersize=10, 
                                linestyle='', label='OGTT2HR')
            
            # Annotate OGTT values
            for _, row in ogtt_data.iterrows():
                color = '#E74C3C' if row['OGTT2HR (mg/dl)'] >= 200 else 'black'
                fontsize = 16 if row['OGTT2HR (mg/dl)'] >= 200 else 12
                ax.annotate(f"{row['OGTT2HR (mg/dl)']:.0f}", 
                           (row['Days'], row['OGTT2HR (mg/dl)']),
                           textcoords="offset points", xytext=(0,15), 
                           ha='center', color=color, fontsize=fontsize, weight='bold' if color == '#E74C3C' else 'normal')
        
        # Plot SSPG
        sspg_data = patient_data.dropna(subset=['SSPG (mg/dl)'])
        if not sspg_data.empty:
            sspg_line, = ax.plot(sspg_data['Days'], sspg_data['SSPG (mg/dl)'], 
                                color='#E74C3C', marker='s', linewidth=2, markersize=10, 
                                linestyle='', label='SSPG')
            
            # Annotate SSPG values
            for _, row in sspg_data.iterrows():
                color = '#E74C3C' if row['SSPG (mg/dl)'] >= 150 else 'black'
                fontsize = 16 if row['SSPG (mg/dl)'] >= 150 else 12
                ax.annotate(f"{row['SSPG (mg/dl)']:.0f}", 
                           (row['Days'], row['SSPG (mg/dl)']),
                           textcoords="offset points", xytext=(0,15), 
                           ha='center', color=color, fontsize=fontsize, weight='bold' if color == '#E74C3C' else 'normal')
        
        # Secondary y-axis for HbA1C
        ax2 = ax.twinx()
        hba1c_data = patient_data.dropna(subset=['HbA1C (%)'])
        if not hba1c_data.empty:
            hba1c_line, = ax2.plot(hba1c_data['Days'], hba1c_data['HbA1C (%)'], 
                                  color='lightgray', marker='o', linewidth=2, markersize=8, label='HbA1C')
            
            # Annotate HbA1C values
            for _, row in hba1c_data.iterrows():
                color = '#E74C3C' if row['HbA1C (%)'] >= 6.5 else 'black'
                fontsize = 16 if row['HbA1C (%)'] >= 6.5 else 12
                ax2.annotate(f"{row['HbA1C (%)']:.1f}", 
                            (row['Days'], row['HbA1C (%)']),
                            textcoords="offset points", xytext=(0,10), 
                            ha='center', color=color, fontsize=fontsize, weight='bold' if color == '#E74C3C' else 'normal')
        
        # Formatting
        ax.set_xlabel('Days in Study', fontsize=14)
        ax.set_ylabel('Glucose Levels (mg/dl)', fontsize=14)
        ax2.set_ylabel('HbA1C (%)', fontsize=14)
        
        # Add reference lines
        ax.axhline(y=126, color='red', linestyle='--', alpha=0.5, label='FPG Diabetic Threshold')
        ax.axhline(y=200, color='blue', linestyle='--', alpha=0.5, label='OGTT Diabetic Threshold')
        ax.axhline(y=150, color='purple', linestyle='--', alpha=0.5, label='SSPG Insulin Resistance')
        ax2.axhline(y=6.5, color='gray', linestyle='--', alpha=0.5, label='HbA1C Diabetic Threshold')
        
        # Legend
        lines1, labels1 = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax.legend(lines1 + lines2, labels1 + labels2, 
                 bbox_to_anchor=(0., 1.02, 1., .102), ncol=4,
                 mode="expand", borderaxespad=0., fontsize=12)
        
        if title:
            plt.title(title, fontsize=16, pad=20)
        else:
            plt.title(f'Longitudinal Clinical Biomarkers - Patient {patient_id}', fontsize=16, pad=20)
        
        plt.tight_layout()
        return fig, ax
    
    def plot_diabetes_progression_heatmap(self, progression_data):
        """Create heatmap showing diabetes progression patterns."""
        patients = list(progression_data.keys())
        markers = ['FPG', 'OGTT2HR', 'HbA1C', 'SSPG']
        
        # Create matrix for heatmap
        progression_matrix = np.zeros((len(patients), len(markers)))
        
        for i, patient in enumerate(patients):
            patient_prog = progression_data[patient]
            for pattern in patient_prog['progression_pattern']:
                for marker in pattern['markers']:
                    if marker in markers:
                        j = markers.index(marker)
                        progression_matrix[i, j] = 1
        
        fig, ax = plt.subplots(figsize=(8, 6))
        sns.heatmap(progression_matrix, 
                   xticklabels=markers, 
                   yticklabels=patients,
                   cmap='Reds', 
                   cbar_kws={'label': 'Diabetic Range Detected'},
                   ax=ax)
        
        plt.title('Diabetes Progression Patterns by Patient', fontsize=14)
        plt.xlabel('Clinical Markers', fontsize=12)
        plt.ylabel('Patient ID', fontsize=12)
        plt.tight_layout()
        return fig, ax
    
    def plot_biomarker_correlations(self, correlation_matrix):
        """Create correlation heatmap for clinical biomarkers."""
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Create mask for upper triangle
        mask = np.triu(np.ones_like(correlation_matrix, dtype=bool))
        
        sns.heatmap(correlation_matrix, 
                   mask=mask,
                   annot=True, 
                   cmap='RdBu_r', 
                   center=0,
                   square=True,
                   fmt='.2f',
                   cbar_kws={'label': 'Correlation Coefficient'},
                   ax=ax)
        
        plt.title('Clinical Biomarker Correlations', fontsize=16)
        plt.tight_layout()
        return fig, ax
    
    def plot_risk_stratification(self, risk_data):
        """Create risk stratification visualization."""
        patients = list(risk_data.keys())
        risk_scores = [risk_data[p]['risk_score'] for p in patients]
        risk_categories = [risk_data[p]['risk_category'] for p in patients]
        
        # Color mapping
        color_map = {'Low Risk': '#2ECC71', 'Moderate Risk': '#F39C12', 'High Risk': '#E74C3C'}
        colors = [color_map[cat] for cat in risk_categories]
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Risk score distribution
        ax1.bar(patients, risk_scores, color=colors)
        ax1.set_xlabel('Patient ID', fontsize=12)
        ax1.set_ylabel('Risk Score', fontsize=12)
        ax1.set_title('Diabetes Risk Scores by Patient', fontsize=14)
        ax1.tick_params(axis='x', rotation=45)
        
        # Risk category pie chart
        category_counts = pd.Series(risk_categories).value_counts()
        ax2.pie(category_counts.values, 
               labels=category_counts.index, 
               colors=[color_map[cat] for cat in category_counts.index],
               autopct='%1.1f%%', 
               startangle=90)
        ax2.set_title('Risk Category Distribution', fontsize=14)
        
        plt.tight_layout()
        return fig, (ax1, ax2)


class ClinicalPredictiveModeling:
    """
    Machine learning models for clinical prediction.
    
    Provides tools for:
    - Diabetes risk prediction
    - Biomarker-based classification
    - Longitudinal outcome prediction
    - Model interpretation and validation
    """
    
    def __init__(self):
        self.models = {}
        self.scalers = {}
    
    def prepare_features(self, data):
        """Prepare features for machine learning models."""
        # Get latest measurements for each patient
        feature_data = []
        
        patients = data['PatientID'].unique()
        for patient in patients:
            patient_data = data[data['PatientID'] == patient].iloc[-1]  # Latest measurement
            
            features = {
                'PatientID': patient,
                'FPG': patient_data.get('FPG (mg/dl)', np.nan),
                'OGTT2HR': patient_data.get('OGTT2HR (mg/dl)', np.nan),
                'SSPG': patient_data.get('SSPG (mg/dl)', np.nan),
                'HbA1C': patient_data.get('HbA1C (%)', np.nan)
            }
            
            # Create diabetes label based on any marker in diabetic range
            diabetes_indicators = [
                features['FPG'] >= 126 if not pd.isna(features['FPG']) else False,
                features['OGTT2HR'] >= 200 if not pd.isna(features['OGTT2HR']) else False,
                features['HbA1C'] >= 6.5 if not pd.isna(features['HbA1C']) else False
            ]
            features['Diabetes'] = int(any(diabetes_indicators))
            
            feature_data.append(features)
        
        return pd.DataFrame(feature_data)
    
    def train_diabetes_classifier(self, feature_data):
        """Train random forest classifier for diabetes prediction."""
        # Prepare data
        feature_cols = ['FPG', 'OGTT2HR', 'SSPG', 'HbA1C']
        X = feature_data[feature_cols].fillna(feature_data[feature_cols].mean())
        y = feature_data['Diabetes']
        
        # Split data
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
        
        # Scale features
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)
        
        # Train model
        model = RandomForestClassifier(n_estimators=100, random_state=42)
        model.fit(X_train_scaled, y_train)
        
        # Evaluate
        y_pred = model.predict(X_test_scaled)
        y_pred_proba = model.predict_proba(X_test_scaled)[:, 1]
        
        # Store model and scaler
        self.models['diabetes_classifier'] = model
        self.scalers['diabetes_classifier'] = scaler
        
        # Return evaluation metrics
        results = {
            'classification_report': classification_report(y_test, y_pred),
            'roc_auc': roc_auc_score(y_test, y_pred_proba),
            'feature_importance': dict(zip(feature_cols, model.feature_importances_)),
            'test_accuracy': model.score(X_test_scaled, y_test)
        }
        
        return results
    
    def predict_diabetes_risk(self, patient_features):
        """Predict diabetes risk for new patient data."""
        if 'diabetes_classifier' not in self.models:
            raise ValueError("Diabetes classifier not trained. Call train_diabetes_classifier first.")
        
        model = self.models['diabetes_classifier']
        scaler = self.scalers['diabetes_classifier']
        
        # Prepare features
        feature_cols = ['FPG', 'OGTT2HR', 'SSPG', 'HbA1C']
        X = patient_features[feature_cols].fillna(patient_features[feature_cols].mean())
        X_scaled = scaler.transform(X)
        
        # Predict
        predictions = model.predict_proba(X_scaled)[:, 1]
        
        return predictions
