from rdkit import Chem
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
    input_df.to_csv('CSV/' + condition + ' processed_data.csv')
    return input_df


def plots(input_df):
    """Generates plot of data"""
    input_df = input_df[input_df["act_class"] != 'intermediate']

# frequency plot of activity class
    sns.countplot(x="act_class", data=input_df, edgecolor="black")
    plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
    plt.ylabel('Frequency', fontsize=14, fontweight='bold')
    # plt.savefig(condition + ' class frequency plot.svg')
    plt.show()

# scatter plot of molecular weight vs logP
    sns.scatterplot(x="mw", y="logp", data=input_df, hue="act_class", size="pIC50", edgecolor="black", alpha=0.7)
    plt.xlabel('MW', fontsize=14, fontweight='bold')
    plt.ylabel('LogP', fontsize=14, fontweight='bold')
    # plt.savefig(condition + ' MW vs LogP plot.svg')
    # plt.show()


def mannwhitneyu_boxplot(input_df):
    """"Generates box plots for each of the lipinkski descriptors split via class and then tests significance"""
    from numpy.random import seed
    from scipy.stats import mannwhitneyu

    input_df = input_df[input_df["act_class"] != 'intermediate']

    lipinkski_descriptors = ['mw', 'logp', 'numHdonors', 'numHacceptors']
    for descriptor in lipinkski_descriptors:

        # box plot
        sns.boxplot(x='act_class', y=descriptor, data=input_df)
        plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
        plt.ylabel(descriptor, fontsize=14, fontweight='bold')
        # plt.savefig(condition + ' ' + descriptor + ' box plot.svg')

        # seed random number generator
        seed(1)

        # filters df based on activity class
        active = input_df[input_df["act_class"] == 'active'][descriptor]
        inactive = input_df[input_df["act_class"] == 'inactive'][descriptor]

        stat, p = mannwhitneyu(active, inactive)

        alpha = 0.05
        if p > alpha:
            interpretation = 'Same distribution (fail to reject H0)'
        else:
            interpretation = 'Different distribution (reject H0)'

        results = pd.DataFrame({'Descriptor': descriptor,
                                'Statistics': stat,
                                'p': p,
                                'alpha': alpha,
                                'Interpretation': interpretation}, index=[0])
        print(results)
        # results.to_csv(condition + ' mannwhitneyu ' + descriptor + '.csv')


condition = 'coronavirus'

df = pd.read_csv('CSV/coronavirus_CHEMBL3927_preprocessed.csv')

generate_descriptors(df)

pIC50(df)

plots(df)

mannwhitneyu_boxplot(df)
