import pandas as pd
from chembl_webresource_client.new_client import new_client
import tk
import easygui


# Input box will likely be used in future

# condition = easygui.enterbox('What condition are you interested in')

# Target is selected and searched within the Chembl db before being turned into a dataframe from dict

condition='coronavirus'
target = new_client.target
target_query=target.search(condition)
targets=pd.DataFrame.from_dict(target_query)

# Target row is subsetted from dataframe and it's chembl id is extracted. Any target can be chosen and searched

selected_target=targets["target_chembl_id"][4]

# Now using the target chembl id the database is searched for activity data but only returning results that have activity using IC50.

activity=new_client.activity
act_query=activity.filter(target_chembl_id=selected_target).filter(standard_type="IC50")
df=pd.DataFrame.from_dict(act_query)


df.to_csv(condition + '_' + selected_target + '_' + 'query_raw.csv', index=False)

df=df[["canonical_smiles", "molecule_chembl_id", "standard_value"]]

df=df[df["standard_value"].notna()]

df["standard_value"] = df["standard_value"].astype("float")

act_class = []
for i in df["standard_value"]:
    if float(i) >=10000:
        act_class.append("inactive")
    elif float(i) <= 1000:
        act_class.append("active")
    else:
        act_class.append("intermediate")

df["act_class"]=act_class

df.to_csv(condition + '_' + selected_target + '_' + 'preprocessed.csv', index=False)
print(df)
