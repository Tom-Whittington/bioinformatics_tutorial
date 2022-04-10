
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors, Lipinski
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

pd.set_option('display.max_columns', 8)


def generate_descriptors(input_df):
    """Generates structure from smiles and then uses structure to generate molecular descriptors"""
    input_df["structure"] = input_df["canonical_smiles"].map(lambda x: Chem.MolFromSmiles(x))
    input_df["mw"] = input_df["structure"].map(lambda x: Chem.Descriptors.MolWt(x))
    input_df["logp"] = input_df["structure"].map(lambda x: Chem.Descriptors.MolLogP(x))
    input_df["numHdonors"] = input_df["structure"].map(lambda x: Chem.Lipinski.NumHDonors(x))
    input_df["numHacceptors"] = input_df["structure"].map(lambda x: Chem.Lipinski.NumHAcceptors(x))
    return input_df


def pIC50(input_df):
    """Caps IC50 at 100000000 before converting to pIC50 to increase uniformity"""
    input_df["standard_value"][input_df["standard_value"] > 100000000] = 100000000
    input_df["pIC50"] = -np.log10(df["standard_value"]*(10**-9))
    input_df.drop(columns='standard_value', inplace=True)
    df.to_csv()
    return input_df

condition = 'coronavirus'

df = pd.read_csv('CSV/coronavirus_CHEMBL3927_preprocessed.csv')

generate_descriptors(df)

pIC50(df)

df = df[df["act_class"] != 'intermediate']

sns.countplot(x="act_class", data=df, edgecolor="black")
plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('Frequency', fontsize=14, fontweight='bold')
plt.savefig(condition + ' class frequency plot.svg')
plt.show()
