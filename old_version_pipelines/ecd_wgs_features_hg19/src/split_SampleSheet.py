import pandas as pd 
import os , re, sys

path = sys.argv[1]

df = pd.read_csv(path)

c = 1
n_sample = 100
for i in range(0, len(df), n_sample):
    start = i
    end = start + n_sample
    df.iloc[start:end].to_csv('{}_{}.csv'.format(path[:-4], c), index = None)
    print(c)
    c += 1
