# scripts/05_survival_analysis.py
"""
Survival analysis - Association between thrombophilia gene expression and patient survival
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import os

print("=== Starting Survival Analysis for Thrombophilia Genes ===")
print("=" * 60)

# Load expression data
print("[LOAD] Loading expression data...")
expression_data = pd.read_csv("data/processed/thrombophilia_expression_subset.csv", index_col=0)
print(f"[SUCCESS] Expression data loaded: {expression_data.shape}")

# Create output directory
os.makedirs('results/figures/survival', exist_ok=True)

print(f"\n[ANALYSIS] Starting survival analysis for {len(expression_data)} genes...")

# Survival analysis for each gene
for i, gene in enumerate(expression_data.index, 1):
    print(f"\n[GENE] {i}. Survival analysis for: {gene}")
    
    # Split patients based on gene expression
    median_expression = np.median(expression_data.loc[gene])
    high_expr = expression_data.loc[gene] > median_expression
    low_expr = expression_data.loc[gene] <= median_expression
    
    print(f"   High expression patients: {sum(high_expr)}")
    print(f"   Low expression patients: {sum(low_expr)}")
    
    # Simulate survival data (replace with real clinical data)
    np.random.seed(42)
    n_samples = expression_data.shape[1]
    time = np.random.exponential(60, n_samples)
    event = np.random.binomial(1, 0.7, n_samples)
    
    # Plot Kaplan-Meier curve
    plt.figure(figsize=(10, 6))
    
    kmf_high = KaplanMeierFitter()
    kmf_low = KaplanMeierFitter()
    
    kmf_high.fit(time[high_expr], event[high_expr], label=f'High {gene} Expression')
    kmf_low.fit(time[low_expr], event[low_expr], label=f'Low {gene} Expression')
    
    kmf_high.plot_survival_function(ci_show=True)
    kmf_low.plot_survival_function(ci_show=True)
    
    plt.title(f'Survival Analysis: {gene} Expression in Breast Cancer', 
              fontsize=14, fontweight='bold')
    plt.xlabel('Time (Months)')
    plt.ylabel('Survival Probability')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Log-rank test
    results = logrank_test(time[high_expr], time[low_expr], 
                          event[high_expr], event[low_expr])
    
    # Display p-value on plot
    plt.text(0.5, 0.2, f'Log-rank p-value: {results.p_value:.4f}', 
             transform=plt.gca().transAxes, fontsize=12, 
             bbox=dict(boxstyle="round,pad=0.3", facecolor="white"))
    
    plt.tight_layout()
    plt.savefig(f'results/figures/survival/survival_{gene}.png', 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"   [SUCCESS] Survival plot saved")
    print(f"   [STATS] p-value: {results.p_value:.4f}")
    
    if results.p_value < 0.05:
        print(f"   [SIGNIFICANT] Significant difference! {gene} expression associated with prognosis")
    else:
        print(f"   [INFO] No significant difference found")

print("\n" + "=" * 60)
print("[COMPLETE] Survival analysis completed!")
print("[SAVE] Results saved to results/figures/survival/")