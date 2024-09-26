import pandas as pd
import numpy as np
import os, re, sys

input_file = str(sys.argv[1])
output_dir = str(sys.argv[2])
    
if os.path.getsize(input_file) == 0:
    
    sample = input_file.split('/')[-1].split('.')[0]
    sample = re.sub('chrchr', 'chr', sample)
    cols = ['Sample'] + [int(i) for i in range(-300,301,1)]
    value = [sample] + [None] * 601
    df = pd.DataFrame([value])
    df.columns = cols
    df.to_csv("{}/{}.csv".format(output_dir, sample), index = None)
    
else:
    df = pd.read_csv(input_file, sep = '\t', header = None)
    df = df[[3]]
    df.columns = ['value']
    df = df[(df['value'] >= -300) & (df['value'] <= 300)]

    group = df['value'].value_counts()
    a = group.index
    b = group.values

    df = pd.DataFrame([a, b]).T
    df.columns = ['value', 'count']

    sample = input_file.split('/')[-1].split('.')[0]
    sample = re.sub('chrchr', 'chr', sample)
    df = df.T

    df.columns = [int(i) for i in df.iloc[0]]

    df.columns.name = None
    df.drop('value', axis = 0, inplace = True)
    df.index = [sample]

    for c in range(-300, 301):
        if int(c) not in df.columns.tolist():
            df[int(c)] = 0

    df = df/df.sum().sum()
    df = df[range(-300, 301)]
    df.reset_index(inplace = True)
    df.rename(columns = {'index': 'Sample'}, inplace = True)
    df.to_csv("{}/{}.csv".format(output_dir, sample), index = None)