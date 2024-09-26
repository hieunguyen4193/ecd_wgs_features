import pandas as pd 
import numpy as np 
import os, re, sys, glob
import pickle

input_dir = sys.argv[1]
output_dir = sys.argv[2]
resource = sys.argv[3]


os.system('mkdir -p {}'.format(output_dir))

dir_list = glob.glob('{}/*_features'.format(input_dir))


cols = pd.read_csv('{}/rerun_samples/NUCLEOSOME.csv'.format(resource)).columns[1:].tolist()
df_nucleosome = pd.DataFrame(columns = cols)

cols = pd.read_csv('{}/rerun_samples/EM.csv'.format(resource)).columns[1:].tolist()
df_em = pd.DataFrame(columns = cols)

cols = pd.read_csv('{}/rerun_samples/FLEN.csv'.format(resource)).columns[1:].tolist()
df_flen = pd.DataFrame(columns = cols)


for dir in dir_list:
    # NUCLEOSOME
    nucleosome = pd.read_csv('{}/GWfeature_Nucleosome.csv'.format(dir))
    nucleosome.rename(columns = {"Sample": "SampleID"}, inplace = True)
    nucleosome = nucleosome[df_nucleosome.columns]
    

    # SAMPLE
    sample = nucleosome['SampleID'][0]
    

    # EM
    em = pd.read_csv('{}/GWfeature_EM.csv'.format(dir)).T
    em.columns = em.iloc[1]
    em = em.iloc[2:]
    em.reset_index(inplace = True)
    em.rename(columns = {"index": "SampleID"}, inplace = True)
    em['SampleID'][0] = sample
    em = em[df_em.columns]
    

    # FLEN
    flen = pd.read_csv('{}/GWfeature_flen.csv'.format(dir))
    flen = flen[flen['size'] != 0].T
    flen.columns = flen.iloc[1]
    flen = flen.iloc[[2]]
    flen.columns = [str(c) for c in flen.columns]
    flen.reset_index(inplace = True)
    flen.rename(columns = {"index": "SampleID"}, inplace = True)
    flen['SampleID'][0] = sample
    for c in df_flen.columns:
        if c not in flen.columns.tolist():
            flen[c] = 0
    flen = flen[df_flen.columns]

    
    df_nucleosome = pd.concat([df_nucleosome, nucleosome])
    df_em = pd.concat([df_em, em])
    df_flen = pd.concat([df_flen, flen])
    
    
df_nucleosome.reset_index(drop = True, inplace = True)
df_em.reset_index(drop = True, inplace = True)
df_flen.reset_index(drop = True, inplace = True)

df_nucleosome.to_csv('{}/NUCLEOSOME.csv'.format(output_dir))
print('--FINISHED NUCLEOSOME--', df_nucleosome.shape)
df_em.to_csv('{}/EM.csv'.format(output_dir))
print('--FINISHED EM--', df_em.shape)
df_flen.to_csv('{}/FLEN.csv'.format(output_dir))
print('--FINISHED FLEN--', df_flen.shape)