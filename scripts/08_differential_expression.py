# scripts/08_differential_expression.py
"""
Differential expression analysis of thrombophilia genes across breast cancer subtypes
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os

print("=== Differential Expression Analysis Across Breast Cancer Subtypes ===")
print("=" * 60)

# Load expression data
print("[LOAD] Loading expression data...")
expression_data = pd.read_csv("data/processed/thrombophilia_expression_subset.csv", index_col=0)
print(f"[SUCCESS] Expression data loaded: {expression_data.shape}")

# Simulate breast cancer subtypes
np.random.seed(42)
n_samples = expression_data.shape[1]
subtypes = np.random.choice(['Luminal_A', 'Luminal_B', 'HER2+', 'Triple_Negative'], 
                           n_samples, p=[0.5, 0.2, 0.15, 0.15])

print(f"\n[DISTRIBUTION] Simulated subtype distribution:")
subtype_counts = pd.Series(subtypes).value_counts()
for subtype, count in subtype_counts.items():
    print(f"   {subtype}: {count} samples ({count/n_samples*100:.1f}%)")

# Create output directories
os.makedirs('results/figures/differential_expression', exist_ok=True)
os.makedirs('results/tables/differential_expression', exist_ok=True)

print(f"\n[ANALYSIS] Starting differential analysis for {len(expression_data)} genes...")

# Store results
results = []

# Analysis for each gene
for gene in expression_data.index:
    print(f"\n[GENE] Analyzing: {gene}")
    
    # Expression data for this gene
    gene_expression = expression_data.loc[gene]
    
    # Create DataFrame for analysis
    df = pd.DataFrame({
        'expression': gene_expression.values,
        'subtype': subtypes
    })
    
    # ANOVA test for comparison between all groups
    groups = [df[df['subtype'] == subtype]['expression'] for subtype in df['subtype'].unique()]
    f_stat, p_value = stats.f_oneway(*groups)
    
    print(f"   ANOVA p-value: {p_value:.4f}")
    
    # If significant, perform pairwise comparisons
    if p_value < 0.05:
        print(f"   [SIGNIFICANT] Significant differences between subtypes!")
        
        # Pairwise comparisons
        significant_pairs = []
        subtype_list = df['subtype'].unique()
        
        # Compare each subtype pair
        for i in range(len(subtype_list)):
            for j in range(i+1, len(subtype_list)):
                group1 = df[df['subtype'] == subtype_list[i]]['expression']
                group2 = df[df['subtype'] == subtype_list[j]]['expression']
                
                # Independent t-test
                t_stat, pair_p_value = stats.ttest_ind(group1, group2)
                
                if pair_p_value < 0.05:
                    pair = f"{subtype_list[i]}_vs_{subtype_list[j]}"
                    significant_pairs.append(pair)
                    print(f"      FOUND: {pair}: p = {pair_p_value:.4f}")
        
        print(f"   Significant differences: {significant_pairs}")
    else:
        significant_pairs = []
        print(f"   [INFO] No significant differences between subtypes")
    
    # Calculate mean expression in each subtype
    means = df.groupby('subtype')['expression'].mean()
    
    # Store results
    results.append({
        'Gene': gene,
        'ANOVA_P_Value': p_value,
        'Significant': p_value < 0.05,
        'Luminal_A_Mean': means.get('Luminal_A', np.nan),
        'Luminal_B_Mean': means.get('Luminal_B', np.nan),
        'HER2+_Mean': means.get('HER2+', np.nan),
        'Triple_Negative_Mean': means.get('Triple_Negative', np.nan),
        'Significant_Pairs': ', '.join(significant_pairs) if significant_pairs else 'None'
    })
    
    # Plot boxplot for significant genes
    if p_value < 0.05:
        plt.figure(figsize=(10, 6))
        sns.boxplot(data=df, x='subtype', y='expression', palette='Set2')
        plt.title(f'Differential Expression of {gene} in Breast Cancer Subtypes\n(ANOVA p = {p_value:.4f})')
        plt.xlabel('Breast Cancer Subtype')
        plt.ylabel('Expression Level')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(f'results/figures/differential_expression/{gene}_expression_by_subtype.png', 
                   dpi=300, bbox_inches='tight')
        plt.close()
        print(f"   [PLOT] Differential expression plot saved")

# Convert to DataFrame
results_df = pd.DataFrame(results)

print("\n" + "=" * 60)
print("[SUMMARY] Differential Expression Results Summary:")
print("=" * 60)

# Significant genes
significant_genes = results_df[results_df['Significant'] == True]
print(f"\n[SIGNIFICANT] Genes with significant differences between subtypes:")
if len(significant_genes) > 0:
    for _, row in significant_genes.iterrows():
        print(f"   FOUND: {row['Gene']}: p = {row['ANOVA_P_Value']:.4f}")
        if row['Significant_Pairs'] != 'None':
            print(f"      Differences: {row['Significant_Pairs']}")
else:
    print("   [INFO] No genes show significant differences between subtypes")

print(f"\n[STATS] Overall statistics:")
print(f"   Significant genes: {len(significant_genes)}/{len(results_df)}")
print(f"   Percentage significant: {len(significant_genes)/len(results_df)*100:.1f}%")

# Save results
results_df.to_csv('results/tables/differential_expression_results.csv', index=False)
print(f"\n[SAVE] Results saved to 'results/tables/differential_expression_results.csv'")

# Create summary plot
plt.figure(figsize=(12, 8))
significant = results_df['Significant'].values
genes = results_df['Gene'].values
p_values = results_df['ANOVA_P_Value'].values

colors = ['red' if sig else 'gray' for sig in significant]
plt.barh(genes, -np.log10(p_values), color=colors)
plt.axvline(-np.log10(0.05), color='red', linestyle='--', alpha=0.7, label='p=0.05')
plt.xlabel('-log10(p-value)')
plt.title('Differential Expression Analysis - Thrombophilia Genes\nAcross Breast Cancer Subtypes')
plt.legend()
plt.tight_layout()
plt.savefig('results/figures/differential_expression_summary.png', dpi=300, bbox_inches='tight')
print("[PLOT] Differential expression summary plot saved")

print("\n" + "=" * 60)
print("[NEXT] Next step: Deep analysis of significant genes")
print("=" * 60)