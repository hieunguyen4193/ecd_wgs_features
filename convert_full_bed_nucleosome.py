import pandas as pd
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

df = pd.read_csv(input_file, header = None, sep = '\t')
df[0] = df[0].astype(str)
df_new = pd.DataFrame()

list_of_chr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
    'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
    'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22',
    'chrX', 'chrY', 'chrM']

for chr in list_of_chr:
    df_new = pd.concat([df_new, df[df[0] == chr].sort_values(1)])
df = df_new.copy()
df.to_csv(output_file, index = None, header = None, sep = '\t')