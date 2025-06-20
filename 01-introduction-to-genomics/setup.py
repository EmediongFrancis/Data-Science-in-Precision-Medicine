#!/usr/bin/env python3
"""
Setup script for Introduction to Genomic Data Analysis
Stanford Data Ocean Precision Medicine Project

Author: Emediong Francis
Date: 2024
"""

import os
import sys
import subprocess

def install_requirements():
    """Install required packages from requirements.txt"""
    print("📦 Installing required packages...")
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "-r", "requirements.txt"])
        print("✅ Requirements installed successfully!")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error installing requirements: {e}")
        return False
    return True

def setup_directories():
    """Create necessary directories if they don't exist"""
    directories = [
        "data/raw",
        "data/processed", 
        "results/plots",
        "results/tables",
        "images"
    ]
    
    print("📁 Setting up project directories...")
    for directory in directories:
        os.makedirs(directory, exist_ok=True)
        print(f"   Created: {directory}")
    
    print("✅ Directories created successfully!")

def create_gitignore():
    """Create .gitignore file for the project"""
    gitignore_content = """
# Data files
*.csv
*.tsv
*.vcf
*.vcf.gz

# Jupyter notebook checkpoints
.ipynb_checkpoints/

# Python cache
__pycache__/
*.pyc
*.pyo

# Results and outputs
results/
*.png
*.jpg
*.pdf

# AWS credentials
.aws/
credentials

# Environment files
.env
.venv/
"""
    
    with open(".gitignore", "w") as f:
        f.write(gitignore_content.strip())
    
    print("✅ .gitignore file created!")

def main():
    """Main setup function"""
    print("🧬 Setting up Introduction to Genomic Data Analysis Project")
    print("=" * 60)
    
    # Check Python version
    if sys.version_info < (3, 8):
        print("❌ Python 3.8 or higher is required!")
        sys.exit(1)
    
    print(f"✅ Python version: {sys.version}")
    
    # Setup project
    setup_directories()
    create_gitignore()
    
    # Install requirements
    if install_requirements():
        print("\n🎉 Project setup completed successfully!")
        print("\n📚 Next steps:")
        print("1. Configure your AWS credentials for Athena access")
        print("2. Open the Jupyter notebook: jupyter lab notebooks/genomic_data_analysis.ipynb")
        print("3. Start analyzing genomic data!")
    else:
        print("\n❌ Setup failed. Please check the error messages above.")
        sys.exit(1)

if __name__ == "__main__":
    main() 