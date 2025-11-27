# generate_figure1.py
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

print("[FIGURE] Starting Figure 1 Generation...")

# Figure settings
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.family'] = 'Arial'

# Create output directory
os.makedirs('../manuscript_figures/Figure1_Expression_Heatmap', exist_ok=True)

# Sample data for testing
genes = ['F2', 'F3', 'F5', 'SERPINC1', 'PROC', 'PROS1', 'THBD', 
         'FGA', 'FGB', 'FGG', 'PLAT', 'PLAU', 'SERPINE1']

# Create sample heatmap
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# Panel A: Heatmap
np.random.seed(42)
data = np.random.randn(len(genes), 50)
sns.heatmap(data, xticklabels=False, yticklabels=genes, 
            cmap='RdBu_r', center=0, ax=ax1)
ax1.set_title('A) Gene Expression Heatmap')
ax1.set_ylabel('Thrombophilia Genes')

# Panel B: Boxplot
subtypes = ['Luminal A', 'Luminal B', 'HER2-enriched', 'Triple-negative']
box_data = [np.random.normal(i, 0.5, 100) for i in range(4)]
ax2.boxplot(box_data, labels=subtypes)
ax2.set_title('B) PLAU Expression by Subtype')
ax2.set_ylabel('Expression Level')

plt.tight_layout()

# Save figure
output_path = '../manuscript_figures/Figure1_Expression_Heatmap/Figure1.tiff'
plt.savefig(output_path, format='tiff', bbox_inches='tight', dpi=300)
plt.close()

print(f"[SUCCESS] Figure 1 created: {output_path}")
print("[OUTPUT] Check the manuscript_figures folder!")