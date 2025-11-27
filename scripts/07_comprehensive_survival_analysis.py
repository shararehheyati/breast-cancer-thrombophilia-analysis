# scripts/07_comprehensive_survival_analysis.py
"""
Comprehensive survival results analysis - Determine if high or low expression is favorable
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

print("=== Comprehensive Survival Analysis for All Genes ===")
print("=" * 60)

# Pre-existing biological knowledge
biological_knowledge = {
    "PLAU": {"role": "Harmful - Promotes metastasis", "expected": "Low expression better"},
    "SERPINE1": {"role": "Harmful - Fibrinolysis inhibitor", "expected": "Low expression better"},
    "F2": {"role": "Neutral - Prothrombin", "expected": "Unclear"},
    "F5": {"role": "Harmful - Factor V Leiden", "expected": "Low expression better"},
    "SERPINC1": {"role": "Beneficial - Antithrombin", "expected": "High expression better"},
    "PROC": {"role": "Beneficial - Protein C", "expected": "High expression better"},
    "PROS1": {"role": "Beneficial - Protein S", "expected": "High expression better"},
    "FGG": {"role": "Neutral - Fibrinogen", "expected": "Unclear"},
    "FGA": {"role": "Neutral - Fibrinogen", "expected": "Unclear"},
    "FGB": {"role": "Neutral - Fibrinogen", "expected": "Unclear"},
    "F3": {"role": "Harmful - Tissue factor", "expected": "Low expression better"},
    "PLAT": {"role": "Beneficial - tPA (fibrinolytic)", "expected": "High expression better"},
    "THBD": {"role": "Beneficial - Thrombomodulin", "expected": "High expression better"}
}

# Survival results from previous analysis (p-values)
survival_results = {
    "F5": 0.4973, "F2": 0.6352, "SERPINC1": 0.2565, "PROC": 0.2933,
    "PROS1": 0.9075, "FGG": 0.1644, "FGA": 0.5633, "FGB": 0.4707,
    "F3": 0.8993, "PLAT": 0.3651, "PLAU": 0.0111, "SERPINE1": 0.6293,
    "THBD": 0.5223
}

print("[ANALYSIS] Gene-by-gene analysis:")
print("=" * 50)

results_summary = []

for gene, p_value in survival_results.items():
    bio_info = biological_knowledge[gene]
    
    # Determine significance
    significance = "Significant" if p_value < 0.05 else "Not significant"
    
    # Symbol based on significance
    symbol = "[SIGNIFICANT]" if p_value < 0.05 else "[INFO]"
    
    print(f"\n{symbol} {gene}:")
    print(f"   Biological role: {bio_info['role']}")
    print(f"   Expected: {bio_info['expected']}")
    print(f"   p-value: {p_value:.4f} ({significance})")
    
    # Save for final summary
    results_summary.append({
        'Gene': gene,
        'Biological_Role': bio_info['role'],
        'Expected': bio_info['expected'],
        'P_Value': p_value,
        'Significance': significance
    })

# Convert to DataFrame for better display
summary_df = pd.DataFrame(results_summary)

print("\n" + "=" * 60)
print("[SUMMARY] Results summary:")
print("=" * 60)

# Significant genes
significant_genes = summary_df[summary_df['Significance'] == 'Significant']
print(f"\n[SIGNIFICANT] Genes with p < 0.05:")
if len(significant_genes) > 0:
    for _, row in significant_genes.iterrows():
        print(f"   FOUND: {row['Gene']}: p = {row['P_Value']:.4f} - {row['Biological_Role']}")
else:
    print("   [INFO] No significant genes found")

# Genes with trend toward significance
trend_genes = summary_df[(summary_df['P_Value'] < 0.2) & (summary_df['P_Value'] >= 0.05)]
print(f"\n[TREND] Genes with trend (0.05 < p < 0.2):")
if len(trend_genes) > 0:
    for _, row in trend_genes.iterrows():
        print(f"   TREND: {row['Gene']}: p = {row['P_Value']:.4f} - {row['Biological_Role']}")
else:
    print("   [INFO] No genes with trend found")

print(f"\n[STATS] Overall statistics:")
print(f"   Total genes: {len(summary_df)}")
print(f"   Significant genes: {len(significant_genes)}")
print(f"   Trend genes: {len(trend_genes)}")
print(f"   Non-significant genes: {len(summary_df) - len(significant_genes) - len(trend_genes)}")

# Save results
summary_df.to_csv('results/tables/survival_analysis_summary.csv', index=False)
print(f"\n[SAVE] Results saved to 'results/tables/survival_analysis_summary.csv'")

print("\n" + "=" * 60)
print("[GUIDE] Instructions for reviewing survival plots:")
print("=" * 60)

print("""
For final interpretation, examine the survival plots:

[ANALYSIS] For each gene answer these questions:
1. Which line is higher? High Expression or Low Expression?
2. Does the pattern match biological expectations?
3. How large is the gap between the two lines?

[EXAMPLE] For PLAU:
   - Expected: Low expression better (harmful gene)
   - In plot: Low Expression line should be higher
   - Conclusion: If confirmed, findings match expectations

[NOTE] Even if p-value is not significant, trends can be interesting!
""")

# Create summary plot
plt.figure(figsize=(12, 8))
genes = summary_df['Gene']
p_values = summary_df['P_Value']
colors = ['red' if p < 0.05 else 'orange' if p < 0.2 else 'gray' for p in p_values]

plt.barh(genes, -np.log10(p_values), color=colors)
plt.axvline(-np.log10(0.05), color='red', linestyle='--', alpha=0.7, label='p=0.05')
plt.axvline(-np.log10(0.2), color='orange', linestyle='--', alpha=0.7, label='p=0.2')
plt.xlabel('-log10(p-value)')
plt.title('Survival Analysis Results - Thrombophilia Genes')
plt.legend()
plt.tight_layout()
plt.savefig('results/figures/survival_summary_plot.png', dpi=300, bbox_inches='tight')
print("\n[PLOT] Summary plot saved to 'results/figures/survival_summary_plot.png'")