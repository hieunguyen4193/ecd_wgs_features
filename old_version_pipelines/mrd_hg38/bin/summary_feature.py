#!/usr/bin/env python3
# coding: utf-8

import pandas as pd 
import numpy as np 
import os, re, sys, glob
import pickle

input_dir = sys.argv[1]
output_dir = sys.argv[2]
feature_type = sys.argv[3]
example_samples = sys.argv[4]


os.system('mkdir -p {}'.format(output_dir))

#####################  GW_feature_v0.2_EM ########################
if feature_type == 'EM':
            file_list = glob.glob('{}/*_GWfeature_EM.csv'.format(input_dir))

            summary = pd.DataFrame()

            for file in file_list:
                b = pd.read_csv(file)
                b = b.iloc[:, 1:]
                summary = pd.concat([summary, b]) if len(summary) == 0 else pd.merge(summary, b, on = 'motif', how = 'left')

            summary = summary.T
            summary.columns = summary.iloc[0]
            summary = summary.iloc[1:]
            summary.reset_index(inplace = True)
            summary.rename(columns = {'index': "SampleID"}, inplace = True)
                
            print(summary.shape)
            cols = pd.read_csv('{}/EM.csv'.format(example_samples), index_col = [0]).columns.tolist()
            summary = summary[cols]
            summary.to_csv("{}/EM.csv".format(output_dir))
            
            
#####################  GW_feature_v0.2_FLEN ########################
if feature_type == 'FLEN':
            file_list = glob.glob('{}/*_GWfeature_FLEN.csv'.format(input_dir))

            summary = pd.DataFrame()

            for file in file_list:
                b = pd.read_csv(file)
                b = b.iloc[:, 1:]
                b = b[b['size'] != 0]
                b.reset_index(drop = True, inplace = True)
                for size in range(50, 351):
                    size = int(size)
                    if size not in b['size'].tolist():
                        b.loc[len(b)] = [size, 0]
                summary = pd.concat([summary, b]) if len(summary) == 0 else pd.merge(summary, b, on = 'size', how = 'left')

            summary.sort_values('size', ascending = True, inplace = True)
            summary = summary.T
            summary.columns = summary.iloc[0]
            summary.columns = [str(c) for c in summary.columns]
            summary = summary.iloc[1:]
            summary.reset_index(inplace = True)
            summary.rename(columns = {'index': "SampleID"}, inplace = True)
                
            print(summary.shape)
            cols = pd.read_csv('{}/FLEN.csv'.format(example_samples), index_col = [0]).columns.tolist()
            summary = summary[cols]
            summary.to_csv("{}/FLEN.csv".format(output_dir))


# #####################  GW_feature_v0.2_CNA_SHORT ########################
# elif feature_type == 'CNA_SHORT':
#             file_list = glob.glob('{}/*_GWfeature_CNA_SHORT.csv'.format(input_dir))

#             summary = pd.DataFrame()

#             for file in file_list:
#                 b = pd.read_csv(file)
#                 sample = b.columns[1].split('.short')[0]
#                 b.columns = ['region', sample]
#                 summary = pd.concat([summary, b]) if len(summary) == 0 else pd.merge(summary, b, on = 'region', how = 'left')

#             summary = summary.T
#             summary.columns = summary.iloc[0]
#             summary = summary.iloc[1:]
#             summary.reset_index(inplace = True)
#             summary.rename(columns = {'index': "SampleID"}, inplace = True)
                
#             print(summary.shape)
#             summary.to_csv("{}/GWfeature_CNA_SHORT.csv".format(output_dir), index = None)
            
            
# #####################  GW_feature_v0.2_CNA ########################
# elif feature_type == 'CNA':
#             file_list = glob.glob('{}/*_GWfeature_CNA.csv'.format(input_dir))

#             summary = pd.DataFrame()

#             for file in file_list:
#                 b = pd.read_csv(file)
#                 sample = b.columns[1].split('.sorted')[0]
#                 b.columns = ['region', sample]
#                 summary = pd.concat([summary, b]) if len(summary) == 0 else pd.merge(summary, b, on = 'region', how = 'left')

#             summary = summary.T
#             summary.columns = summary.iloc[0]
#             summary = summary.iloc[1:]
#             summary.reset_index(inplace = True)
#             summary.rename(columns = {'index': "SampleID"}, inplace = True)
                
#             print(summary.shape)
#             summary.to_csv("{}/GWfeature_CNA.csv".format(output_dir), index = None)

# # #####################  GW_feature_v0.2_FLEN_RATIO_SHORT_TOTAL ########################
# elif feature_type == 'FLEN_RATIO':
#             for name, name_to_save in [['ratio_short_long', 'FLEN_RATIO_SHORT_LONG'],
#                                         ['ratio_short_total', 'FLEN_RATIO_SHORT_TOTAL']]:
#                 file_list = glob.glob('{}/*_GWfeature_FLEN_RATIO.csv'.format(input_dir))

#                 summary = pd.DataFrame()

#                 for file in file_list:
#                     b = pd.read_csv(file)
#                     sample = b.columns[-1].split('.')[0]
#                     b = b[['region', '{}.{}'.format(sample, name)]]
#                     summary = pd.concat([summary, b]) if len(summary) == 0 else pd.merge(summary, b, on = 'region', how = 'left')

#                 summary = summary.T
#                 summary.columns = summary.iloc[0]
#                 summary = summary.iloc[1:]
#                 summary.reset_index(inplace = True)
#                 summary.rename(columns = {'index': "SampleID"}, inplace = True)
                    
#                 print(summary.shape)
#                 summary.to_csv("{}/GWfeature_{}.csv".format(output_dir, name_to_save), index = None)