#!/usr/bin/env python3
# coding: utf-8

import pandas as pd
import numpy as np
import os, re
import sys

def get_tumor_fraction(path):
    f = open(path, "r")
    f.readline()
    data = f.readline().split('\t')
    return float(data[1])

input_dir = sys.argv[1]
output_dir = sys.argv[2]


save = []
for file in os.listdir(input_dir):
    if file[-11:] == '.params.txt':
        sample = file[:-11].split('/')[-1]
        tumor_faction = get_tumor_fraction(input_dir + '/' + file)
        save.append([sample, tumor_faction])

df = pd.DataFrame(save)
df.columns = ['SampleID', 'ichorCNA']

df.to_csv('{}/IchorCNA.csv'.format(output_dir))
