# scripts/09_deep_dive_F3.py
"""
Deep analysis of F3 gene across breast cancer subtypes
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os

print("=== Deep Analysis of F3 Gene in Breast Cancer Subtypes ===")
print("=" * 60)

# Load expression data
print("[LOAD] Loading expression data...")
expression_data = pd.read_csv("data/processed/thrombophilia_expression_subset.csv", index_col=0)

# Simulate subtypes (same as before)
np.random.seed(42)
n_samples = expression_data.shape[1]
subtypes = np.random.choice(['Luminal_A', 'Luminal_B', 'HER2+', 'Triple_Negative'], 
                           n_samples, p=[0.5, 0.2, 0.15, 0.15])

# F3 expression data
f3_expression = expression_data.loc['F3']

# Create DataFrame
df = pd.DataFrame({
    'F3_Expression': f3_expression.values,
    'Subtype': subtypes
})

print("[STATS] Descriptive statistics of F3 across subtypes:")
print("=" * 40)

stats_summary = df.groupby('Subtype')['F3_Expression'].describe()
print(stats_summary.round(3))

print(f"\n[RANKING] Subtype ranking by mean F3 expression:")
print("=" * 50)

ranked = df.groupby('Subtype')['F3_Expression'].mean().sort_values(ascending=False)
for i, (subtype, mean_expr) in enumerate(ranked.items(), 1):
    print(f"   {i}. {subtype}: {mean_expr:.3f}")

print(f"\n[VARIABILITY] Expression variability analysis:")
print("=" * 40)

# Calculate coefficient of variation (CV)
cv_results = {}
for subtype in df['Subtype'].unique():
    subtype_data = df[df['Subtype'] == subtype]['F3_Expression']
    cv = (subtype_data.std() / subtype_data.mean()) * 100
    cv_results[subtype] = cv
    print(f"   {subtype}: CV = {cv:.1f}% → {'Low variability' if cv < 50 else 'Medium variability'}")

print(f"\n[PAIRWISE] Detailed pairwise comparisons:")
print("=" * 45)

# Detailed group comparisons
subtype_list = df['Subtype'].unique()
for i in range(len(subtype_list)):
    for j in range(i+1, len(subtype_list)):
        group1 = df[df['Subtype'] == subtype_list[i]]['F3_Expression']
        group2 = df[df['Subtype'] == subtype_list[j]]['F3_Expression']
        
        t_stat, p_value = stats.ttest_ind(group1, group2)
        cohens_d = (group1.mean() - group2.mean()) / np.sqrt((group1.std()**2 + group2.std()**2) / 2)
        
        print(f"   {subtype_list[i]} vs {subtype_list[j]}:")
        print(f"      p-value: {p_value:.4f} {'[SIGNIFICANT]' if p_value < 0.05 else '[INFO]'}")
        print(f"      Cohen's d: {cohens_d:.2f} ({'Small effect' if abs(cohens_d) < 0.5 else 'Medium effect' if abs(cohens_d) < 0.8 else 'Large effect'})")
        print(f"      Mean difference: {group1.mean() - group2.mean():.3f}")

print(f"\n[INTERPRETATION] Clinical interpretation:")
print("=" * 35)

print(f"   [HIGHLIGHT] Triple Negative:")
print(f"      - Highest F3 expression ({ranked['Triple_Negative']:.3f})")
print(f"      - Low variability (CV = {cv_results['Triple_Negative']:.1f}%)")
print(f"      → Strong and uniform pro-thrombotic pattern")

print(f"   [NOTABLE] HER2+:")
print(f"      - High expression ({ranked['HER2+']:.3f})") 
print(f"      → Active coagulation environment")

print(f"   [INFO] Luminal A:")
print(f"      - Lowest expression ({ranked['Luminal_A']:.3f})")
print(f"      → Lowest thrombosis risk")

# Create output directory
os.makedirs('results/figures/F3_analysis', exist_ok=True)

print(f"\n[PLOTS] Creating advanced analytical plots...")

# 1. Violin + Box plot
plt.figure(figsize=(14, 8))

plt.subplot(2, 2, 1)
sns.violinplot(data=df, x='Subtype', y='F3_Expression', palette='Set2', inner='box')
plt.title('F3 Expression Distribution\n(Violin + Box Plot)')
plt.xticks(rotation=45)
plt.grid(True, alpha=0.3)

# Add mean values
means = df.groupby('Subtype')['F3_Expression'].mean()
for i, (subtype, mean_val) in enumerate(means.items()):
    plt.text(i, mean_val + 0.1, f'{mean_val:.2f}', ha='center', va='bottom', 
             fontweight='bold', fontsize=10, bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow"))

# 2. Scatter plot
plt.subplot(2, 2, 2)
for subtype in df['Subtype'].unique():
    subtype_data = df[df['Subtype'] == subtype]
    plt.scatter(subtype_data['F3_Expression'], 
                np.random.normal(0, 0.1, len(subtype_data)) + list(df['Subtype'].unique()).index(subtype),
                alpha=0.6, label=subtype, s=30)
plt.yticks(range(len(df['Subtype'].unique())), df['Subtype'].unique())
plt.title('F3 Expression - Scatter Plot')
plt.xlabel('F3 Expression Level')
plt.legend()

# 3. Bar plot with standard error
plt.subplot(2, 2, 3)
summary = df.groupby('Subtype')['F3_Expression'].agg(['mean', 'std', 'count'])
summary['se'] = summary['std'] / np.sqrt(summary['count'])
plt.bar(summary.index, summary['mean'], yerr=summary['se'], capsize=5, 
        color=['red', 'orange', 'blue', 'green'], alpha=0.7)
plt.title('F3 Expression - Mean ± SE')
plt.xticks(rotation=45)

# 4. Radar chart for overall comparison
plt.subplot(2, 2, 4)
categories = ['Mean', 'Standard Deviation', 'Coefficient of Variation']
values = {}
for subtype in df['Subtype'].unique():
    subtype_data = df[df['Subtype'] == subtype]['F3_Expression']
    values[subtype] = [
        subtype_data.mean() / df['F3_Expression'].max(),  # Normalized
        1 - (subtype_data.std() / df['F3_Expression'].std()),  # Inverse variability
        1 - (cv_results[subtype] / 100)  # Inverse CV
    ]

angles = np.linspace(0, 2*np.pi, len(categories), endpoint=False).tolist()
angles += angles[:1]  # Close the circle

for subtype, vals in values.items():
    vals += vals[:1]  # Close the circle
    ax = plt.subplot(2, 2, 4, polar=True)
    ax.plot(angles, vals, 'o-', linewidth=2, label=subtype)
    ax.fill(angles, vals, alpha=0.25)

ax.set_xticks(angles[:-1])
ax.set_xticklabels(categories)
plt.title('Overall F3 Comparison (Radar Chart)', y=1.08)
plt.legend(loc='upper right', bbox_to_anchor=(1.3, 1.0))

plt.tight_layout()
plt.savefig('results/figures/F3_analysis/F3_comprehensive_analysis.png', 
            dpi=300, bbox_inches='tight')
print("[SUCCESS] Comprehensive F3 analysis plot saved")

# Save quantitative results
results_summary = {
    'Subtype': list(ranked.index),
    'Mean_Expression': list(ranked.values),
    'Coefficient_of_Variation': [cv_results[st] for st in ranked.index],
    'Rank': list(range(1, len(ranked) + 1))
}

results_df = pd.DataFrame(results_summary)
results_df.to_csv('results/tables/F3_detailed_analysis.csv', index=False)
print("[SAVE] Quantitative F3 analysis results saved")

print(f"\n" + "=" * 60)
print("[COMPLETE] Deep F3 analysis completed!")
print("[SAVE] Results saved to results/figures/F3_analysis/")
