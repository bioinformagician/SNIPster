"""
Example usage of the portable path system for genomics project.

This file demonstrates how to use the updated path system to load various data files.
"""

import sys
from pathlib import Path

# Add the web_app directory to the path
sys.path.append(str(Path(__file__).parent / "web_app"))

from utilities import (
    get_genome_file_path,
    list_available_genome_files,
    get_background_data_path,
    gwas_data_path
)
import pandas as pd

def example_load_genome_file():
    """Example of how to load a genome file."""
    print("Available genome files:")
    genome_files = list_available_genome_files()
    
    for genome_file in genome_files:
        print(f"  - {genome_file}")
    
    if genome_files:
        # Load the first available genome file as an example
        first_genome = genome_files[0]
        genome_path = get_genome_file_path(first_genome)
        print(f"\nLoading genome file: {genome_path}")
        
        # Here you would load the genome file (example with pandas)
        # genome_data = pd.read_csv(genome_path, sep='\t')  # Adjust separator as needed
        # print(f"Genome data shape: {genome_data.shape}")
        print("Genome file path resolved successfully!")

def example_load_gwas_data():
    """Example of how to load GWAS data."""
    print(f"\nGWAS data path: {gwas_data_path}")
    
    if gwas_data_path.exists():
        # gwas_data = pd.read_csv(gwas_data_path)
        # print(f"GWAS data shape: {gwas_data.shape}")
        print("GWAS data path resolved successfully!")
    else:
        print("GWAS data file not found!")

def example_load_background_data():
    """Example of how to load background data."""
    sports_snps_path = get_background_data_path("sports_snps.txt")
    print(f"\nSports SNPs path: {sports_snps_path}")
    
    if sports_snps_path.exists():
        # sports_data = pd.read_csv(sports_snps_path, sep='\t')
        # print(f"Sports data shape: {sports_data.shape}")
        print("Sports SNPs data path resolved successfully!")
    else:
        print("Sports SNPs file not found!")

if __name__ == "__main__":
    print("🔬 Genomics Project - Example Usage")
    print("=" * 40)
    
    example_load_genome_file()
    example_load_gwas_data()
    example_load_background_data()
    
    print("\n" + "=" * 40)
    print("✅ All path examples completed!")
    print("\nNow your code is portable and will work regardless of where")
    print("someone clones your repository!")
