#!/usr/bin/env python
# coding: utf-8

# ## Data Import and treatment

# In[1]:


import pandas as pd
from chembl_webresource_client.new_client import new_client
import tk
import easygui


# Input box will likely be used in future

# In[ ]:


# condition = easygui.enterbox('What condition are you interested in')


# Target is selected and searched within the Chembl db before being turned into a dataframe from dict

# In[5]:


condition='coronavirus'
target = new_client.target
target_query=target.search(condition)
targets=pd.DataFrame.from_dict(target_query)
targets


# Target row is subsetted from dataframe and it's chembl id is extracted. Any target can be chosen and searched

# In[8]:


selected_target=targets["target_chembl_id"][4]
selected_target


# Now using the target chembl id the database is searched for activity data but only returning results that have activity using IC50.

# In[23]:


activity=new_client.activity
act_query=activity.filter(target_chembl_id=selected_target).filter(standard_type="IC50")
df=pd.DataFrame.from_dict(act_query)


# Raw query results are saved to CSV

# In[34]:


df.to_csv(condition + '_' + selected_target + '_' + 'query_raw.csv', index=False)


# In[20]:


df=df[["canonical_smiles", "molecule_chembl_id", "standard_value"]]
df


# In[ ]:




