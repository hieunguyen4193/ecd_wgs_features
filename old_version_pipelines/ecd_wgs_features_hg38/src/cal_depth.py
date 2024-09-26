#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import os, re, glob, sys
import subprocess

sample_id = sys.argv[1]
input_file = sys.argv[2]
output_dir = sys.argv[3]
n_thread = sys.argv[4]


script = 'samtools view -@ {} -c {}'.format(n_thread, input_file)


result = subprocess.Popen(script.split(' '), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
data, _ = result.communicate()
data = int(data[:-1])
data

df = pd.DataFrame([sample_id, data]).T
df.columns = ['SampleID', 'N_read']
df['Depth'] = [int(n) * 85 / 3e9 for n in df['N_read']]
df.to_csv('{}/{}_depth_info.csv'.format(output_dir, sample_id), index = None)
df

