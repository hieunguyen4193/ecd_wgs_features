#!/usr/bin/env python3
# coding: utf-8

# In[4]:


import pandas as pd
import os, re
import sys

input_dir = str(sys.argv[1])
output_dir = str(sys.argv[2])
example_samples = str(sys.argv[3])

sample_list = []
for file in os.listdir(input_dir):
    if file[-4:] == '.csv':
        sample_list.append(file.split('.')[0])

sample = sample_list[0]
df = pd.read_csv('{}/{}.csv'.format(input_dir, sample))
df['Sample'][0] = sample
df.to_csv('{}/{}.csv'.format(input_dir, sample), index = None)

count = 0
for sample in sample_list[1:]:
    df_temp = pd.read_csv('{}/{}.csv'.format(input_dir, sample))
    df_temp['Sample'][0] = sample
    df_temp.to_csv('{}/{}.csv'.format(input_dir, sample), index = None)
    df = pd.concat([df, df_temp], axis = 0)
    count += 1
    print(count)
df.columns = [str(c) for c in df.columns]
df.rename(columns = {"Sample": "SampleID"}, inplace = True)

cols = pd.read_csv('{}/NUCLEOSOME.csv'.format(example_samples), index_col = [0]).columns.tolist()
df = df[cols]
df.to_csv('{}/NUCLEOSOME.csv'.format(output_dir))


# In[ ]:




