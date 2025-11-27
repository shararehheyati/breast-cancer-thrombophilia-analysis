# scripts/10_multivariate_survival.py
"""
Multivariate survival analysis for PLAU gene
Objective: Determine if PLAU is an independent biomarker
Learning: Cox Proportional Hazards Model and Hazard Ratio
"""

# Import libraries
import pandas as pd  # For tabular data manipulation
import numpy as np   # For numerical computations
from lifelines import CoxPHFitter  # For Cox survival analysis
import matplotlib.pyplot as plt  # For plotting
import os

print("=== Starting Multivariate Survival Analysis for PLAU ===")
print("=" * 60)

# Explanation of analysis importance
print("[OBJECTIVE] Analysis goals:")
print("   - Determine if PLAU predicts prognosis independently of clinical factors")
print("   - If yes: PLAU is a strong independent biomarker")
print("   - If no: PLAU effect might be due to correlation with stage or grade")

# Load expression data
print("\n[LOAD] Loading expression data...")
expression_data = pd.read_csv("data/processed/thrombophilia_expression_subset.csv", index_col=0)
plau_expression = expression_data.loc['PLAU']  # Only PLAU gene

print(f"[SUCCESS] PLAU data loaded: {len(plau_expression)} samples")

# Create simulated clinical data for demonstration
print("\n[DATA] Creating simulated clinical data...")
np.random.seed(42)  # For reproducible results
n_samples = len(plau_expression)

# Simulated clinical data for training
clinical_data = pd.DataFrame({
    'patient_id': range(n_samples),
    'PLAU_expression': plau_expression.values,
    
    # Important clinical factors in breast cancer:
    'age': np.random.normal(58, 12, n_samples),  # Patient age - mean 58 years
    'stage': np.random.choice(['I', 'II', 'III', 'IV'], n_samples, p=[0.2, 0.5, 0.25, 0.05]),
    'grade': np.random.choice([1, 2, 3], n_samples, p=[0.3, 0.5, 0.2]),
    
    # Survival data
    'time': np.random.exponential(60, n_samples),  # Survival time in months
    'event': np.random.binomial(1, 0.7, n_samples)  # 1 = death, 0 = censored
})

print("[PREVIEW] Sample of simulated data:")
print(clinical_data.head())

print("\n[STATS] Descriptive statistics:")
print(f"   - Number of patients: {len(clinical_data)}")
print(f"   - Mean age: {clinical_data['age'].mean():.1f} years")
print(f"   - Stage distribution: {clinical_data['stage'].value_counts().to_dict()}")
print(f"   - Grade distribution: {clinical_data['grade'].value_counts().to_dict()}")

# Convert categorical variables to numerical for model
print("\n[PREP] Converting categorical variables to numerical...")
clinical_data['stage_num'] = clinical_data['stage'].map({'I': 1, 'II': 2, 'III': 3, 'IV': 4})
clinical_data['grade_num'] = clinical_data['grade']

print("[SUCCESS] Stage converted to numeric: I=1, II=2, III=3, IV=4")
print("[SUCCESS] Grade converted to numeric: 1, 2, 3")

print(f"\n[ANALYSIS] Starting statistical analyses...")

# 1. Univariate analysis - PLAU only
print("\n" + "="*50)
print("1. Univariate Analysis (PLAU only)")
print("="*50)
print("[QUESTION] Does PLAU alone associate with survival?")

cph_univariate = CoxPHFitter()
cph_univariate.fit(clinical_data[['time', 'event', 'PLAU_expression']], 
                  duration_col='time', event_col='event')

print("[RESULTS] Univariate analysis results:")
print(cph_univariate.summary)

# 2. Multivariate analysis - All factors
print("\n" + "="*50)
print("2. Multivariate Analysis (PLAU + age + stage + grade)")
print("="*50)
print("[QUESTION] Does PLAU have independent effect after adjusting for age, stage, grade?")

cph_multivariate = CoxPHFitter()
cph_multivariate.fit(clinical_data[['time', 'event', 'PLAU_expression', 'age', 'stage_num', 'grade_num']], 
                    duration_col='time', event_col='event')

print("[RESULTS] Multivariate analysis results:")
print(cph_multivariate.summary)

# Compare the two models
print("\n" + "="*50)
print("3. Comparison: Univariate vs Multivariate Models")
print("="*50)

hr_uni = np.exp(cph_univariate.params_['PLAU_expression'])
hr_multi = np.exp(cph_multivariate.params_['PLAU_expression'])

print(f"[COMPARISON] Univariate Hazard Ratio (PLAU): {hr_uni:.3f}")
print(f"[COMPARISON] Multivariate Hazard Ratio (PLAU): {hr_multi:.3f}")
print(f"[COMPARISON] HR change: {((hr_multi - hr_uni) / hr_uni * 100):+.1f}%")

# Create plots
print("\n[PLOTS] Generating analysis plots...")
plt.figure(figsize=(15, 6))

plt.subplot(1, 2, 1)
cph_univariate.plot()
plt.title('Univariate Analysis - PLAU Only\n(HR = {:.3f})'.format(hr_uni))

plt.subplot(1, 2, 2)
cph_multivariate.plot()
plt.title('Multivariate Analysis - PLAU + Clinical Factors\n(HR = {:.3f})'.format(hr_multi))

plt.tight_layout()
plt.savefig('results/figures/multivariate_survival_analysis.png', dpi=300, bbox_inches='tight')
print("[SUCCESS] Analysis plots saved")

# Results interpretation
print("\n" + "="*50)
print("4. Results Interpretation for Manuscript")
print("="*50)

plau_p_multivariate = cph_multivariate.summary['p'][0]

if plau_p_multivariate < 0.05:
    print("[SIGNIFICANT] **Key finding:** PLAU is an independent biomarker!")
    print(f"   - Even after adjusting for age, stage, grade")
    print(f"   - Hazard Ratio: {hr_multi:.3f} (p = {plau_p_multivariate:.4f})")
    print(f"   - Patients with high PLAU expression have {((hr_multi-1)*100):+.1f}% higher risk")
    
    print(f"\n[CLINICAL] **Clinical interpretation:**")
    print(f"   - PLAU can be used for risk stratification")
    print(f"   - High PLAU patients may need more aggressive treatment")
else:
    print("[INFO] PLAU effect may depend on other factors")
    print("   - Observed effect might be due to correlation with stage or grade")



print(f"\n[OUTPUT] **Saved files:**")
print(f"   - Complete results: results/tables/multivariate_analysis_data.csv")
print(f"   - Plots: results/figures/multivariate_survival_analysis.png")

# Save results
os.makedirs('results/tables', exist_ok=True)
clinical_data.to_csv('results/tables/multivariate_analysis_data.csv', index=False)
print(f"\n[SAVE] Analysis data saved")

print(f"\n[COMPLETE] Multivariate analysis completed!")
print("[NEXT] Next step: Pathway Analysis")
print("\n" + "="*60)