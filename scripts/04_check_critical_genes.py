# scripts/01_load_and_explore_data.py
"""
First analysis script: Load and explore thrombophilia gene expression in breast cancer
Data source: TCGA-BRCA RNA-seq dataset
"""

import pandas as pd
import os

print("=== Starting Analysis: Thrombophilia Genes in Breast Cancer ===")
print("=" * 60)

# 1. Setup paths
print("\n[SETUP] Setting up file paths...")
data_dir = "data"
raw_data_dir = os.path.join(data_dir, "raw")
gene_list_path = os.path.join(data_dir, "thrombophilia_genes.txt")

# 2. Load thrombophilia gene list
print("[LOAD] Loading thrombophilia gene list...")
try:
    with open(gene_list_path, 'r') as f:
        thrombophilia_genes = [line.strip() for line in f if line.strip()]
    print(f"[SUCCESS] Found {len(thrombophilia_genes)} thrombophilia genes:")
    for i, gene in enumerate(thrombophilia_genes, 1):
        print(f"   {i:2d}. {gene}")
except FileNotFoundError:
    print("[ERROR] Gene list file not found!")
    print(f"   Make sure this file exists: {gene_list_path}")
    exit()

# 3. Check if raw data file exists
expression_file_path = os.path.join(raw_data_dir, "TCGA_BRCA_RNAseq_Expression.txt")
print(f"\n[CHECK] Checking for expression data file...")
print(f"   Looking for: {expression_file_path}")

if not os.path.exists(expression_file_path):
    print("[ERROR] Expression data file not found!")
    print("   Please download TCGA data from GDC portal")
    exit()
else:
    print("[SUCCESS] Expression data file found!")

# 4. Load and explore data structure
print("\n[LOAD] Loading expression data (preview)...")
print("   This might take a moment for large files...")

try:
    # Preview file structure
    with open(expression_file_path, 'r') as f:
        first_lines = [f.readline() for _ in range(5)]
    
    print("\n[PREVIEW] File structure:")
    for i, line in enumerate(first_lines):
        print(f"   Line {i}: {line[:100]}..." if len(line) > 100 else f"   Line {i}: {line.strip()}")
    
    # Load full dataset
    print("\n[LOAD] Loading full dataset...")
    expression_df = pd.read_csv(expression_file_path, sep='\t', index_col=0)
    
    print("[SUCCESS] Data loaded successfully!")
    print(f"   Dataset shape: {expression_df.shape} (genes × samples)")
    print(f"   Number of genes: {expression_df.shape[0]}")
    print(f"   Number of samples: {expression_df.shape[1]}")
    
    # 5. Check thrombophilia genes in dataset
    print(f"\n[CHECK] Finding thrombophilia genes in dataset...")
    genes_found = [gene for gene in thrombophilia_genes if gene in expression_df.index]
    genes_not_found = [gene for gene in thrombophilia_genes if gene not in expression_df.index]
    
    print(f"[SUCCESS] Found {len(genes_found)} genes in dataset:")
    for gene in genes_found:
        print(f"   FOUND: {gene}")
    
    if genes_not_found:
        print(f"[WARNING] {len(genes_not_found)} genes not found:")
        for gene in genes_not_found:
            print(f"   MISSING: {gene}")

    # 6. Basic statistics for thrombophilia genes
    if genes_found:
        print(f"\n[STATS] Basic statistics for thrombophilia genes:")
        thrombophilia_data = expression_df.loc[genes_found]
        
        print("   Mean expression across samples:")
        for gene in genes_found:
            mean_expr = thrombophilia_data.loc[gene].mean()
            print(f"      {gene}: {mean_expr:.2f}")
            
        # Save processed subset
        output_dir = os.path.join(data_dir, "processed")
        os.makedirs(output_dir, exist_ok=True)
        
        output_file = os.path.join(output_dir, "thrombophilia_expression_subset.csv")
        thrombophilia_data.to_csv(output_file)
        print(f"\n[SAVE] Saved expression data to: {output_file}")

except Exception as e:
    print(f"[ERROR] Loading data: {e}")
    print("   Check file format and memory availability")

print("\n" + "=" * 60)
print("[COMPLETE] Script finished successfully!")
print("   Next: Run statistical analysis scripts")