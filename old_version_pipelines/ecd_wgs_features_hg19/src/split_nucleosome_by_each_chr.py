import pandas as pd
import sys

input_file = sys.argv[1]
chr = sys.argv[2]
output_file = sys.argv[3]

df = pd.read_csv(input_file, header = None, sep = '\t')
df = df[df[0] == chr].sort_values(1)
df.to_csv(output_file, index = None, header = None, sep = '\t')
