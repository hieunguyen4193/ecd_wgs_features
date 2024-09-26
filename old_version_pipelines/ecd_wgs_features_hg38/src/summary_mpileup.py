import pandas as pd
import os
import glob
import re
import sys

input_dir = sys.argv[1]
output_dir = sys.argv[2]

file_list = glob.glob('{}/*.csv'.format(input_dir))

df = pd.DataFrame()
for file in file_list:
    df_temp = pd.read_csv(file, header = None)
    df = pd.concat([df, df_temp])
    
    
df.columns = ['SampleID', 'Mpileup']
df['Breadth of coverage'] = df['Mpileup'] / 3e9
df.to_csv('{}/Mpileup.csv'.format(output_dir), index = None)