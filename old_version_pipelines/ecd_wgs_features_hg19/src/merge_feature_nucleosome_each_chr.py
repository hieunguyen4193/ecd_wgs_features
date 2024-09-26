import pandas as pd
import sys
import glob

input_dir = sys.argv[1]
sampleid = sys.argv[2]
output_dir = sys.argv[3]

df = pd.DataFrame()

for chr in range(1, 23):
    file_path = '{}/{}.chr.{}.csv'.format(input_dir, sampleid, chr)
    df_now = pd.read_csv(file_path)
    df_now['Sample'] = 'chr.{}'.format(chr)
    
    if len(df) == 0:
        df = df_now.copy()
    else:
        df_now = df_now[df.columns]
        df = pd.concat([df, df_now])

df = df.T
df.columns = df.iloc[0]
df = df.iloc[1:]
df.reset_index(inplace = True)
df.rename(columns = {'index': "size"}, inplace = True)
df.to_csv('{}/{}.csv'.format(output_dir, sampleid))