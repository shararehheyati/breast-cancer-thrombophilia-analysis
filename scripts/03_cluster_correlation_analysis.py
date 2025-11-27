# scripts/03_cluster_correlation_analysis.py
"""
Advanced clustering and correlation analysis for thrombophilia genes
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

print("=== Starting Cluster & Correlation Analysis ===")
print("=" * 60)

# Professional settings
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.family'] = 'Arial'
os.makedirs('results/figures', exist_ok=True)
os.makedirs('results/tables', exist_ok=True)

# Load processed data
print("\n[LOAD] Loading processed data...")
try:
    expression_data = pd.read_csv("data/processed/thrombophilia_expression_subset.csv", index_col=0)
    print(f"[SUCCESS] Data loaded! {expression_data.shape[0]} genes, {expression_data.shape[1]} samples")
except FileNotFoundError:
    print("[ERROR] Processed data not found! Run 01_load_and_explore_data.py first")
    exit()

# 1. Clustered Heatmap
print("\n[PLOT] 1. Creating Clustered Heatmap...")
plt.figure(figsize=(12, 8))
clustered_grid = sns.clustermap(expression_data,
                               cmap='RdBu_r',
                               center=0,
                               yticklabels=True,
                               xticklabels=False)
plt.suptitle('Clustered Heatmap of Thrombophilia Genes', fontsize=14)
plt.savefig('results/figures/Figure5_clustered_heatmap.png', dpi=300, bbox_inches='tight')
print("[SUCCESS] Figure 5: Clustered Heatmap saved!")

# 2. Correlation Heatmap
print("\n[PLOT] 2. Creating Correlation Heatmap...")
correlation_matrix = expression_data.T.corr()
plt.figure(figsize=(10, 8))
sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', center=0, fmt='.2f')
plt.title('Correlation Matrix of Thrombophilia Genes', fontsize=14)
plt.tight_layout()
plt.savefig('results/figures/Figure6_correlation_heatmap.png', dpi=300, bbox_inches='tight')
print("[SUCCESS] Figure 6: Correlation Heatmap saved!")

print("\n[COMPLETE] Analysis Completed! Check results/figures/ folder.")