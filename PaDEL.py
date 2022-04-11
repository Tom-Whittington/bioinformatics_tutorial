import padelpy
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
import os

df_full=pd.read_csv('CSV/coronavirus processed_data.csv')
df=df_full[["canonical_smiles", "molecule_chembl_id"]]
df.to_csv("molecule.smi", sep='\t', index=False, header=False)

padelpy.padeldescriptor(mol_dir='molecule.smi',
                        fingerprints=True,
                        removesalt=True,
                        standardizenitro=True,
                        descriptortypes='References/PubchemFingerprinter.xml',
                        d_file='descriptors.csv')

df_X = pd.read_csv('descriptors.csv')

files = ['molecule.smi', 'descriptors.csv']
for file in files:
    os.remove(file)

df_X = df_X.drop(columns=["Name"])
df_y=df_full["pIC50"]
df=pd.concat([df_X,df_y], axis=1)
X = df.drop('pIC50', axis=1)
y = df['pIC50']

print(X.shape, y.shape)

from sklearn.feature_selection import VarianceThreshold
selection=VarianceThreshold(threshold=(.8*(1-.8)))
X=selection.fit_transform(X)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

import numpy as np
np.random.seed(100)
model = RandomForestRegressor(n_estimators=100)
model.fit(X_train,y_train)
r2 = model.score(X_test, y_test)
y_pred = model.predict(X_test)

sns.set(color_codes=True)
sns.set_style("white")

ax = sns.regplot(y_test, y_pred, scatter_kws={'alpha':0.4})
ax.set_xlabel('Experimental pIC50', fontsize='large', fontweight='bold')
ax.set_ylabel('Predicted pIC50', fontsize='large', fontweight='bold')
ax.set_xlim(0, 12)
ax.set_ylim(0, 12)
ax.figure.set_size_inches(5,5)
plt.show()






