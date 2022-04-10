import pandas as pd
from chembl_webresource_client.new_client import new_client
import tk
import easygui

def classifier(x):
    """Classifies molecule activity into active, intermediate and inactive vis IC50 value"""
    if x >= 10000:
        return "inactive"
    elif x <= 1000:
        return "active"
    else:
        return "intermediate"


# TODO: Add input box to allow selecting conditons

# condition = easygui.enterbox('What condition are you interested in')

# Target is selected and searched within the Chembl db before being turned into a dataframe from dict

condition = 'coronavirus'
target = new_client.target
target_query = target.search(condition)
targets = pd.DataFrame.from_dict(target_query)

# Target row is subsetted from dataframe and it's chembl id is extracted. Any target can be chosen and searched

# TODO: Add way of choosing target and saving target name

selected_target = targets["target_chembl_id"][4]

# Now using the target chembl id the database is searched for activity data but only returning results that have
# activity using IC50.

activity = new_client.activity
act_query = activity.filter(target_chembl_id=selected_target).filter(standard_type="IC50")
df = pd.DataFrame.from_dict(act_query)

# df.to_csv(condition + '_' + selected_target + '_' + 'query_raw.csv', index=False)

# Unneeded columns and NAs are dropped and standard_value is renamed to IC50
df = df[["canonical_smiles", "molecule_chembl_id", "standard_value"]]
df.rename(columns={'standard_value' : 'IC50'}, inplace=True)
df = df[df["IC50"].notna()]
df["IC50"] = df["IC50"].astype("float")

# Molecules are classified via IC50
df["act_class"] = df["IC50"].apply(lambda x: classifier(x))

print(df)
# df.to_csv(condition + '_' + selected_target + '_' + 'preprocessed.csv', index=False)
