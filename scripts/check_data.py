import os
import pandas as pd
import numpy as np

print("[CHECK] Checking project data structure...")

# 1. Check folder structure
print("\n" + "="*50)
print("[STRUCTURE] Project folder structure:")
print("="*50)

for root, dirs, files in os.walk("."):
    # Only show important folders
    if 'data' in root or 'results' in root or 'scripts' in root:
        level = root.replace('.', '').count(os.sep)
        indent = ' ' * 2 * level
        print(f"{indent}[FOLDER] {os.path.basename(root)}/")
        subindent = ' ' * 2 * (level + 1)
        for file in files:
            if file.endswith(('.csv', '.tsv', '.txt', '.py')):
                print(f"{subindent}[FILE] {file}")

# 2. Check gene expression data
print("\n" + "="*50)
print("[EXPRESSION] Checking gene expression data:")
print("="*50)

try:
    expression_files = []
    for root, dirs, files in os.walk("data"):
        for file in files:
            if 'expression' in file.lower() or 'gene' in file.lower():
                expression_files.append(os.path.join(root, file))
    
    if expression_files:
        print("[SUCCESS] Gene expression files found:")
        for file in expression_files:
            print(f"   [FILE] {file}")
            
        # Load first expression file
        expression_data = pd.read_csv(expression_files[0], index_col=0)
        print(f"   Data shape: {expression_data.shape}")
        print(f"   Number of genes: {len(expression_data)}")
        print(f"   Number of samples: {len(expression_data.columns)}")
        print("\n   Data sample:")
        print(expression_data.iloc[:5, :5])
    else:
        print("[ERROR] No gene expression files found")

except Exception as e:
    print(f"[ERROR] Loading expression data: {e}")

# 3. Check clinical data
print("\n" + "="*50)
print("[CLINICAL] Checking clinical data:")
print("="*50)

try:
    clinical_files = []
    for root, dirs, files in os.walk("data"):
        for file in files:
            if 'clinical' in file.lower() or 'survival' in file.lower() or 'stage' in file.lower():
                clinical_files.append(os.path.join(root, file))
    
    if clinical_files:
        print("[SUCCESS] Clinical files found:")
        for file in clinical_files:
            print(f"   [FILE] {file}")
            
        # Load first clinical file
        clinical_data = pd.read_csv(clinical_files[0])
        print(f"   Data shape: {clinical_data.shape}")
        print(f"   Available columns: {list(clinical_data.columns)}")
        print("\n   Data sample:")
        print(clinical_data.head())
    else:
        print("[ERROR] No clinical files found")

except Exception as e:
    print(f"[ERROR] Loading clinical data: {e}")

# 4. Check PLAU gene in data
print("\n" + "="*50)
print("[PLAU] Checking PLAU gene:")
print("="*50)

try:
    if 'expression_data' in locals():
        if 'PLAU' in expression_data.index:
            plau_expression = expression_data.loc['PLAU']
            print(f"[SUCCESS] PLAU gene found!")
            print(f"   Mean expression: {plau_expression.mean():.2f}")
            print(f"   Expression range: {plau_expression.min():.2f} - {plau_expression.max():.2f}")
        else:
            print("[ERROR] PLAU gene not found in data")
            print("   Available genes (first 10):")
            print(expression_data.index[:10].tolist())
except Exception as e:
    print(f"[ERROR] Checking PLAU gene: {e}")

print("\n" + "="*50)
print("[COMPLETE] Data check completed!")
print("="*50)