# scripts/01_load_and_explore_data.py
import pandas as pd
import os

print("=== Starting Analysis: Thrombophilia Genes in Breast Cancer ===")

# Safe path - put file in same folder as script
script_dir = os.path.dirname(__file__)
gene_list_path = os.path.join(script_dir, "thrombophilia_genes.txt")

print(f"[INFO] Looking for gene list at: {gene_list_path}")

# 1. First check if file exists
if not os.path.exists(gene_list_path):
    print("[ERROR] Gene list file not found!")
    print("[SOLUTION] Create thrombophilia_genes.txt in scripts folder")
    print("File content:")
    print("F2")
    print("F3") 
    print("F5")
    print("SERPINC1")
    print("PROC")
    print("PROS1")
    print("THBD")
    print("FGA")
    print("FGB")
    print("FGG")
    print("PLAT")
    print("PLAU")
    print("SERPINE1")
    exit()

# 2. If file exists, continue
print("[SUCCESS] Gene list file found!")
with open(gene_list_path, 'r') as f:
    thrombophilia_genes = [line.strip() for line in f if line.strip()]
    
print(f"Found {len(thrombophilia_genes)} genes:")
for gene in thrombophilia_genes:
    print(f"  - {gene}")

print("\n[COMPLETE] Script finished successfully!")