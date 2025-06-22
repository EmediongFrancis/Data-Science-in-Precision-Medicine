# Getting Started with Clinical Data Analysis

**Stanford Data Science in Precision Medicine - Module 5**

This guide will help you quickly get started with clinical data analysis using our comprehensive toolkit.

## 🚀 Quick Setup

### 1. Installation
```bash
# Install dependencies
pip install -r requirements.txt

# Install the module
pip install -e .
```

### 2. Verify Installation
```bash
# Run the demo script
python demo.py
```

## 📊 First Analysis

### Basic Clinical Data Loading
```python
from scripts.clinical_utils import ClinicalDataProcessor
from scripts.config import DEMO_CONFIG

# Initialize processor
processor = ClinicalDataProcessor()

# Generate demo data (or load your own)
clinical_data = processor.generate_demo_data(DEMO_CONFIG)

# Basic data exploration
print(f"Data shape: {clinical_data.shape}")
print(clinical_data.head())
```

### Quality Control Check
```python
# Perform quality control
qc_results = processor.quality_control(clinical_data)

# View results
print(f"Completeness: {qc_results['completeness_rate']:.1f}%")
print(f"Missing values: {qc_results['missing_values']}")
```

### Diabetes Progression Analysis
```python
from scripts.clinical_utils import ClinicalBiomarkerAnalyzer

# Initialize analyzer
analyzer = ClinicalBiomarkerAnalyzer()

# Analyze progression patterns
progression = analyzer.diabetes_progression_analysis(clinical_data)

# View results for each patient
for patient, results in progression.items():
    print(f"Patient {patient}: {results['first_diabetic_marker']}")
```

## 📈 Common Use Cases

### 1. Individual Patient Analysis
```python
from scripts.clinical_utils import ClinicalDataVisualizer

# Initialize visualizer
visualizer = ClinicalDataVisualizer()

# Plot patient trajectory
patient_id = "ZNDMXI3"
fig, ax = visualizer.plot_longitudinal_biomarkers(clinical_data, patient_id)
plt.show()
```

### 2. Risk Stratification
```python
# Calculate risk scores
risk_results = analyzer.risk_stratification(clinical_data)

# View risk categories
for patient, risk in risk_results.items():
    print(f"{patient}: {risk['risk_category']} ({risk['risk_score']} points)")
```

### 3. Predictive Modeling
```python
from scripts.clinical_utils import ClinicalPredictiveModeling

# Initialize ML model
ml_model = ClinicalPredictiveModeling()

# Prepare features and train model
features = ml_model.prepare_features(clinical_data)
results = ml_model.train_diabetes_classifier(features)

print(f"Model accuracy: {results['test_accuracy']:.3f}")
```

## �� Clinical Scenarios

### Scenario 1: OGTT-First Progression
```python
# Focus on patient with OGTT-first progression
patient_data = clinical_data[clinical_data['PatientID'] == 'ZNDMXI3']

# Analyze temporal trends
trends = processor.temporal_analysis(patient_data)
print("OGTT progression pattern:", trends['OGTT2HR (mg/dl)'])
```

### Scenario 2: Infection-Triggered Diabetes
```python
# Analyze viral infection impact
patient_data = clinical_data[clinical_data['PatientID'] == 'ZOZOWIT']

# Plot recovery and relapse pattern
fig, ax = visualizer.plot_longitudinal_biomarkers(
    clinical_data, 'ZOZOWIT',
    title="Viral Infection Impact on Glucose Metabolism"
)
plt.show()
```

## 🔧 Troubleshooting

### Common Issues

**1. Missing Dependencies**
```bash
# If you get import errors
pip install --upgrade -r requirements.txt
```

**2. Data Loading Problems**
```python
# Check data structure
print(clinical_data.columns.tolist())
print(clinical_data.dtypes)
```

**3. Visualization Issues**
```python
# Ensure matplotlib backend is set
import matplotlib
matplotlib.use('TkAgg')  # or 'Qt5Agg'
```

### Performance Optimization

**1. Large Datasets**
```python
# Use chunked processing for large files
chunk_size = 1000
for chunk in pd.read_csv('large_file.csv', chunksize=chunk_size):
    # Process chunk
    pass
```

**2. Memory Management**
```python
# Optimize data types
clinical_data = clinical_data.astype({
    'PatientID': 'category',
    'Days': 'int16'
})
```

## 📚 Learning Path

### Beginner
1. Start with the demo script (`python demo.py`)
2. Explore the Jupyter notebook
3. Try basic data loading and visualization

### Intermediate
1. Implement custom quality control checks
2. Create your own visualization functions
3. Experiment with different ML models

### Advanced
1. Integrate with clinical databases
2. Develop custom biomarker algorithms
3. Build clinical decision support tools

## 🎯 Key Concepts

### Clinical Markers
- **FPG**: Fasting Plasma Glucose (homeostasis)
- **OGTT2HR**: Oral Glucose Tolerance Test (processing)
- **SSPG**: Steady-State Plasma Glucose (resistance)
- **HbA1C**: Hemoglobin A1C (long-term control)

### Reference Ranges
- **Normal FPG**: 70-99 mg/dl
- **Normal OGTT**: <140 mg/dl
- **Normal HbA1C**: 4.0-5.6%
- **Insulin Sensitive**: SSPG <150 mg/dl

### Diabetes Thresholds
- **FPG**: ≥126 mg/dl
- **OGTT2HR**: ≥200 mg/dl
- **HbA1C**: ≥6.5%

## 🤝 Getting Help

### Documentation
- Check the main README.md for comprehensive documentation
- Review the notebook for step-by-step examples
- Examine the demo script for complete workflows

### Common Questions
1. **How do I load my own data?** Use `processor.load_clinical_data('your_file.csv')`
2. **How do I interpret results?** Check the clinical insights section in the demo
3. **How do I customize visualizations?** Modify the `ClinicalDataVisualizer` class methods

### Support
For technical support or questions about the module, please refer to the main documentation or contact the Stanford Data Science in Precision Medicine team.

---

**Next Steps**: Once you're comfortable with the basics, dive into the comprehensive Jupyter notebook for detailed analysis and clinical interpretation!
