from setuptools import setup, find_packages

setup(
    name="introduction-to-epigenetics",
    version="1.0.0",
    description="Introduction to Epigenetics - Stanford Data Ocean Module",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Stanford Data Ocean",
    author_email="dataocean@stanford.edu",
    url="https://github.com/stanford-dataocean/data-science-precision-medicine",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Education",
    ],
    python_requires=">=3.7",
    install_requires=[
        "pandas>=1.3.0",
        "numpy>=1.21.0",
        "matplotlib>=3.4.0",
        "seaborn>=0.11.0",
        "scipy>=1.7.0",
        "scikit-learn>=1.0.0",
        "jupyter>=1.0.0",
        "notebook>=6.0.0",
        "ipywidgets>=7.6.0",
    ],
    extras_require={
        "advanced": [
            "pybedtools>=0.8.0",
            "pysam>=0.16.0",
            "biopython>=1.78",
            "openpyxl>=3.0.0",
        ],
        "visualization": [
            "plotly>=5.0.0",
            "bokeh>=2.3.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "epigenetics-demo=demo:main",
        ],
    },
    keywords="epigenetics bioinformatics genomics data-science education",
    project_urls={
        "Documentation": "https://github.com/stanford-dataocean/data-science-precision-medicine",
        "Source": "https://github.com/stanford-dataocean/data-science-precision-medicine",
        "Tracker": "https://github.com/stanford-dataocean/data-science-precision-medicine/issues",
    },
) 