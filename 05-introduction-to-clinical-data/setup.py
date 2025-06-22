from setuptools import setup, find_packages

setup(
    name="clinical-data-analysis",
    version="1.0.0",
    description="Clinical Data Analysis for Stanford Data Science in Precision Medicine",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Stanford Data Science in Precision Medicine Team",
    author_email="precision-medicine@stanford.edu",
    url="https://github.com/stanford-medicine/precision-medicine-modules",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Education",
        "Intended Audience :: Healthcare Industry",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Education",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    keywords=[
        "clinical-data",
        "precision-medicine", 
        "biomarkers",
        "diabetes",
        "longitudinal-analysis",
        "machine-learning",
        "healthcare-analytics",
        "ipop-study",
        "stanford-medicine"
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.21.0",
        "pandas>=1.3.0",
        "scipy>=1.7.0",
        "scikit-learn>=1.0.0",
        "matplotlib>=3.4.0",
        "seaborn>=0.11.0",
        "statsmodels>=0.12.0",
        "tqdm>=4.62.0",
        "jupyter>=1.0.0"
    ],
    extras_require={
        "dev": [
            "pytest>=6.0.0",
            "black>=21.0.0",
            "flake8>=3.9.0",
            "mypy>=0.910"
        ],
        "advanced": [
            "xgboost>=1.5.0",
            "plotly>=5.0.0",
            "lifelines>=0.27.0",
            "pingouin>=0.5.0"
        ]
    },
    entry_points={
        "console_scripts": [
            "clinical-demo=demo:main",
        ],
    },
    include_package_data=True,
    package_data={
        "": ["*.md", "*.txt", "*.yml", "*.yaml"],
    },
    project_urls={
        "Documentation": "https://stanford-medicine.github.io/precision-medicine-modules/",
        "Source": "https://github.com/stanford-medicine/precision-medicine-modules",
        "Tracker": "https://github.com/stanford-medicine/precision-medicine-modules/issues",
    },
)
