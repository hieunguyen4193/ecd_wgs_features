#!/usr/bin/env python3
# coding: utf-8

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
    df_temp = pd.read_csv(file)
    df = pd.concat([df, df_temp])
    
df.to_csv('{}/Depth.csv'.format(output_dir))