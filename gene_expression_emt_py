#Install & Import Libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
from scipy.integrate import odeint
from sklearn.linear_model import RidgeCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix
import networkx as nx
#Load GEO Dataset 
try:
    gse = GEOparse.get_GEO("GSE45827", destdir="./")
    data = gse.pivot_samples("VALUE")
    expr = data.values
    
    meta = gse.phenotype_data
    print("Available columns:", meta.columns)

    # a usable column 
    col = None
    for c in meta.columns:
        if "characteristics" in c.lower():
            col = c
            break

    if col:
        meta_col = meta[col].astype(str)
    else:
        print(" No characteristics column found")
        meta_col = pd.Series([''] * len(meta))
group = np.where(
        meta_col.str.contains('basal|her2', case=False, na=False),
        1, 0
    )
 print(" GEO data loaded successfully")
except Exception as e:
    print(" GEO load failed:", e)
    
    np.random.seed(0)
    expr = np.random.randn(1000, 50)
    group = np.random.randint(0, 2, 50)
print(meta_col.unique()[:50])

for col in meta.columns:
    print("\nCOLUMN:", col)
    print(meta[col].astype(str).unique()[:10])
new one
from scipy.stats import ttest_ind
import numpy as np

# Use correct column
meta_col = meta['characteristics_ch1.1.tumor subtype'].astype(str)

# Filter only Basal & HER2
mask = meta_col.str.contains('Basal|Her2', case=False, na=False)

expr_filtered = expr[:, mask]
meta_filtered = meta_col[mask]

# Create groups
group = np.where(
    meta_filtered.str.contains('Basal', case=False),
    1, 0
)

print("Group counts:", np.bincount(group))

3. Differential Expression # Split
group1 = expr_filtered[:, group == 1]
group0 = expr_filtered[:, group == 0]

print("Shapes:", group1.shape, group0.shape)

# T-test
t_stat, p_vals = ttest_ind(group1, group0, axis=1, equal_var=False)

# LogFC
logFC = group1.mean(axis=1) - group0.mean(axis=1)

# Significant genes
sig_idx = (p_vals < 0.05) & (np.abs(logFC) > 1)

# Top genes
top_n = min(30, np.sum(sig_idx))
top_genes_idx = np.where(sig_idx)[0][:top_n]

# Final matrix
X = expr_filtered[top_genes_idx].T

from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

pca = PCA(n_components=2)
X_pca = pca.fit_transform(X)

plt.scatter(X_pca[:,0], X_pca[:,1], c=group)
plt.title("PCA: Basal vs HER2")
plt.show()


from sklearn.ensemble import RandomForestClassifier

model = RandomForestClassifier()
model.fit(X, group)

print("Training accuracy:", model.score(X, group))

import seaborn as sns

sns.heatmap(X[:20], cmap='viridis')
plt.title("Top Differential Genes")
plt.show()
4. Ridge-Based Network Inference
A = np.zeros((X.shape[1], X.shape[1]))

for i in range(X.shape[1]):
    y = X[:, i]
    X_pred = np.delete(X, i, axis=1)
    model = RidgeCV(alphas=[0.1, 1, 10]).fit(X_pred, y)
    coefs = model.coef_
    A[i, np.arange(X.shape[1]) != i] = coefs
#Gene Co-expression Heatmap
cor_matrix = np.corrcoef(X, rowvar=False)
sns.heatmap(cor_matrix)
plt.title("Gene Co-expression")
plt.show()

#Network Graph (Top 10% edges)
threshold = np.quantile(np.abs(cor_matrix), 0.9)
adj = np.where(np.abs(cor_matrix) > threshold, cor_matrix, 0)
G = nx.from_numpy_array(adj)
nx.draw(G, node_size=50)
plt.title("Gene Network")
plt.show()

#Random Forest Classification
X_train, X_test, y_train, y_test = train_test_split(X, group, test_size=0.3)
rf = RandomForestClassifier(n_estimators=100)
rf.fit(X_train, y_train)
pred = rf.predict(X_test)
print("Confusion Matrix:")
print(confusion_matrix(y_test, pred))

# Feature importance
importances = rf.feature_importances_
plt.bar(range(len(importances)), importances)
plt.title("Feature Importance")
plt.show()

#Mechanistic ODE Model (SMAD, SNAIL, ZEB, E-cadherin)
def emt_model(y, t, params):
    SMAD, SNAIL, ZEB, ECAD = y
    k1, k2, k3, k4, d1, d2, d3, d4 = params
    
    dSMAD = k1 - d1 * SMAD
    dSNAIL = k2 * (SMAD**2 / (1 + SMAD**2)) - d2 * SNAIL
    dZEB = k3 * (SNAIL**2 / (1 + SNAIL**2)) - d3 * ZEB
    dECAD = k4 * (1 / (1 + SNAIL**2 + ZEB**2)) - d4 * ECAD
    
    return [dSMAD, dSNAIL, dZEB, dECAD]

# Initial conditions
y0 = [0.1, 0.1, 0.1, 1.0]
t = np.linspace(0, 100, 200)
params = [1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5]

solution = odeint(emt_model, y0, t, args=(params,))

# Plot dynamics
plt.plot(t, solution[:,0], label='SMAD')
plt.plot(t, solution[:,1], label='SNAIL')
plt.plot(t, solution[:,2], label='ZEB')
plt.plot(t, solution[:,3], label='E-cadherin')
plt.xlabel("Time")
plt.ylabel("Expression")
plt.title("TGF-beta Induced EMT Dynamics")
plt.legend()
plt.show()
#Decision-Making Simulation (Vary TGF-beta)
k1_values = [0.5, 1, 2, 3]
for k1 in k1_values:
    params = [k1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5]
    sol = odeint(emt_model, y0, t, args=(params,))
    plt.plot(t, sol[:,3], label=f'TGF={k1}')  # E-cadherin

plt.xlabel("Time")
plt.ylabel("E-cadherin")
plt.title("Cell Fate Switch under TGF-beta Levels")
plt.legend()
plt.show()

