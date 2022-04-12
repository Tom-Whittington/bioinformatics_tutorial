from chembl_webresource_client.new_client import new_client
from sklearn.feature_selection import VarianceThreshold
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from numpy.random import seed
from scipy.stats import mannwhitneyu
from rdkit.Chem import Descriptors, Lipinski
from rdkit import Chem
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import padelpy
import os
import tk
import easygui




def search(condition):
    """Target is selected and searched within the Chembl db before being turned into a dataframe from dict"""

    # TODO: Add input box to allow selecting conditions

    target = new_client.target
    target_query = target.search(condition)
    df = pd.DataFrame.from_dict(target_query)

    selected_target = (df["target_chembl_id"].iloc[4])

    # Now using the target chembl id the database is searched for activity data but only returning results that have
    # activity using IC50.

    activity = new_client.activity
    act_query = activity.filter(target_chembl_id=selected_target).filter(standard_type="IC50")
    df = pd.DataFrame.from_dict(act_query)
    # df.to_csv(condition + '_' + selected_target + '_' + 'query_raw.csv', index=False)

    return df


def search_molecules(targets):
    # TODO: Add way of choosing target and saving target name
    selected_target = targets["target_chembl_id"][4]

    # Now using the target chembl id the database is searched for activity data but only returning results that have
    # activity using IC50.

    activity = new_client.activity
    act_query = activity.filter(target_chembl_id=selected_target).filter(standard_type="IC50")
    df = pd.DataFrame.from_dict(act_query)
    # df.to_csv(condition + '_' + selected_target + '_' + 'query_raw.csv', index=False)

    return targets


def formatter(input_df):
    """Unneeded columns and NAs are dropped and standard_value is renamed to IC50"""
    input_df = input_df[["canonical_smiles", "molecule_chembl_id", "standard_value"]]
    input_df.rename(columns={'standard_value': 'IC50'}, inplace=True)
    input_df = input_df[input_df["IC50"].notna()]
    input_df["IC50"] = input_df["IC50"].astype("float")
    input_df["act_class"] = input_df["IC50"].apply(lambda x: act_classifier(x))

    return input_df

def formatter(input_df):
    """Unneeded columns and NAs are dropped and standard_value is renamed to IC50"""
    input_df = input_df[["canonical_smiles", "molecule_chembl_id", "standard_value"]]
    input_df.rename(columns={'standard_value': 'IC50'}, inplace=True)
    input_df = input_df[input_df["IC50"].notna()]
    input_df["IC50"] = input_df["IC50"].astype("float")

    return input_df


def act_classifier(x):
    """Classifies molecule activity into active, intermediate and inactive vis IC50 value"""
    if x >= 10000:
        return "inactive"
    elif x <= 1000:
        return "active"
    else:
        return "intermediate"


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
    input_df["IC50"][input_df["IC50"] > 100000000] = 100000000
    input_df["pIC50"] = -np.log10(input_df["IC50"] * (10 ** -9))
    input_df.drop(columns='IC50', inplace=True)
    # input_df.to_csv('CSV/' + condition + ' processed_data.csv')
    return input_df


def plots(input_df, condition):
    """Generates plot of data"""
    plot_df = input_df[input_df["act_class"] != 'intermediate']

    # frequency plot of activity class
    sns.countplot(x="act_class", data=plot_df, edgecolor="black")
    plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
    plt.ylabel('Frequency', fontsize=14, fontweight='bold')
    # plt.savefig(condition + ' class frequency plot.svg')
    plt.show()

    # scatter plot of molecular weight vs logP
    sns.scatterplot(x="mw", y="logp", data=plot_df, hue="act_class", size="pIC50", edgecolor="black", alpha=0.7)
    plt.xlabel('MW', fontsize=14, fontweight='bold')
    plt.ylabel('LogP', fontsize=14, fontweight='bold')
    # plt.savefig(condition + ' MW vs LogP plot.svg')
    plt.show()


def mannwhitneyu_boxplot(input_df):
    """"Generates box plots for each of the lipinkski descriptors split via class and then tests significance"""


    plot_df = input_df[input_df["act_class"] != 'intermediate']

    lipinkski_descriptors = ['mw', 'logp', 'numHdonors', 'numHacceptors']
    for descriptor in lipinkski_descriptors:

        # box plot
        sns.boxplot(x='act_class', y=descriptor, data=plot_df)
        plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
        plt.ylabel(descriptor, fontsize=14, fontweight='bold')
        # plt.savefig(condition + ' ' + descriptor + ' box plot.svg')
        plt.show()

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


# condition = easygui.enterbox('What condition are you interested in')

def padel(df):
    df = df[["canonical_smiles", "molecule_chembl_id"]]
    df.to_csv("molecule.smi", sep='\t', index=False, header=False)

    padelpy.padeldescriptor(mol_dir='molecule.smi',
                            fingerprints=True,
                            removesalt=True,
                            standardizenitro=True,
                            descriptortypes='References/PubchemFingerprinter.xml',
                            d_file=condition + ' descriptors.csv')

    df_desc = pd.read_csv('descriptors.csv')
    os.remove('molecule.smi')
    df_X = df_desc.drop(columns=["Name"])
    df_y = df["pIC50"]
    df_final = pd.concat([df_X, df_y], axis=1)
    X = df_final.drop('pIC50', axis=1)
    y = df_final['pIC50']

    return X, y


def random_forest(X, y):
    selection = VarianceThreshold(threshold=(.8 * (1 - .8)))
    X = selection.fit_transform(X)

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
    np.random.seed(100)
    model = RandomForestRegressor(n_estimators=100)
    model.fit(X_train, y_train)

    r2 = model.score(X_test, y_test)
    y_pred = model.predict(X_test)

    return y_test, y_pred


def random_forest_plot(y_test, y_pred):
    sns.set(color_codes=True)
    sns.set_style("white")

    ax = sns.regplot(y_test, y_pred, scatter_kws={'alpha': 0.4})
    ax.set_xlabel('Experimental pIC50', fontsize='large', fontweight='bold')
    ax.set_ylabel('Predicted pIC50', fontsize='large', fontweight='bold')
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 12)
    ax.figure.set_size_inches(5, 5)
    plt.show()


def main():
    condition = 'coronavirus'
    # search(condition)
    # search_molecules(targets)
    df=search(condition)
    #search_molecules(df)
    df=formatter(df)

    df["act_class"] = df["IC50"].apply(lambda x: act_classifier(x))
    generate_descriptors(df)
    pIC50(df)
    plots(df, condition)
    mannwhitneyu_boxplot(df)
    padel(df)
    random_forest(X, y)
    random_forest_plot(y_test, y_pred)

if __name__ == "__main__":
    main()
