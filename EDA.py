
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors, Lipinski
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import PandasTools
import numpy as np
import pandas as pd

df=pd.read_csv('CSV/coronavirus_CHEMBL3927_preprocessed.csv')
print(df)

def lipinksi(smiles, verbose = False):

    moldata=[]
    for elem in smiles:
        mol=Chem.MolFromSmiles(elem)
        moldata.append(mol)

    baseData = np.arrange(1,1)
    i=0
    for mol in moldata:

        desc_MolWt = Descriptors.MolWt(mol)
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_NumHDonors = Lipinksi.NumHdonors(mol)
        desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)

# PandasTools.AddMoleculeColumnToFrame(df, smilesCol='canonical_smiles')

df["structure"]=df["canonical_smiles"].map(lambda x: Chem.MolFromSmiles(x))
df["molwt"]=df["structure"].map(lambda x: Chem.Descriptors.MolWt(x))
df["logp"]=df["structure"].map(lambda x: Chem.Descriptors.MolLogP(x))
# df[["molwt", "logp", "numHdonor", "numHacceptor"]]=df["structure"].map(lambda x: [Chem.Descriptors.MolWt(x),)


# mol_list = []
# for elem in df["canonical_smiles"]:
#     mol=Chem.MolFromSmiles(elem)
#     mol_list.append(mol)

# df["structure"]=mol_list

print(df)

# MolWt_list = []
# molLogP_list = []
# molNumHDonor_list = []
# molNumHAcceptor_list = []
#
# for mol in df["structure"]:
#     molWt=Descriptors.MolWt(mol)
#     molLogP = Descriptors.MolLogP(mol)
#     molNumHDonor = Lipinksi.NumHdonors(mol)
#     molNumHAcceptor = Lipinski.NumHAcceptors(mol)

#
#
# # In[38]:
#
#
# mol_list[4]
#
#
# # In[ ]:
#
#
# desc_MolWt
#
#
# # In[ ]:
#
#
# desc_MolLogP
#
#
# # In[ ]:
#
#
# desc_NumHDonors

