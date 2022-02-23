#!/usr/bin/env python
# coding: utf-8

# In[30]:


# Install RDKit.
#!pip install rdkit-pypi


# In[31]:


from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Descriptors, Lipinski
from rdkit.Chem import AllChem
from rdkit import DataStructs
import numpy as np
import pandas as pd


# In[32]:


df=pd.read_csv('coronavirus_CHEMBL3927_preprocessed.csv')


# def lipinksi(smiles, verbose = False)
# 
#     moldata=[]
#     for elem in smiles:
#         mol=Chem.MolFromSmiles(elem)
#         moldata.append(mol)
#     
#     baseData = np.arrange(1,1)
#     i=0
#     for mol in moldata:
#         
#         desc_MolWt = Descriptors.MolWt(mol)
#         desc_MolLogP = Descriptors.MolLogP(mol)
#         desc_NumHDonors = Lipinksi.NumHdonors(mol)
#         desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)

# In[35]:


mol_list = []
for elem in df["canonical_smiles"]:
    mol=Chem.MolFromSmiles(elem)
    mol_list.append(mol)


# In[38]:


mol_list[4]


# In[ ]:


desc_MolWt


# In[ ]:


desc_MolLogP


# In[ ]:


desc_NumHDonors

