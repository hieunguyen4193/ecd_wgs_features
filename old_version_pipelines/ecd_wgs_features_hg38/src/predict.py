#!/usr/bin/env python3
# coding: utf-8
import pandas as pd
import numpy as np
import pathlib
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn import metrics
from sklearn.metrics import RocCurveDisplay
from sklearn.decomposition import NMF
import warnings
warnings.filterwarnings("ignore")
from sklearn.mixture import GaussianMixture
import os
import ot
import sys
import pickle

#####---------------------------------------------------------------------------#####
##### input args
#####---------------------------------------------------------------------------#####
path_to_input_feature_folder = sys.argv[1]
path_to_model_files = sys.argv[2]
path_to_save_output = sys.argv[3]
path_to_example_samples = sys.argv[4]

os.system('mkdir -p {}'.format(path_to_save_output))
#####---------------------------------------------------------------------------#####
##### MAIN
#####---------------------------------------------------------------------------#####

print(path_to_input_feature_folder)
rerun_samples = pd.read_csv(path_to_example_samples + '/EM.csv')['SampleID'].tolist()

featuredf = dict()
for feat in ["FLEN", "IchorCNA", "NUCLEOSOME", "EM"]:
    path_to_feature = [item for item in pathlib.Path(path_to_input_feature_folder).glob("*{}*".format(feat))][0]
    featuredf[feat] = pd.read_csv(path_to_feature, index_col = [0])
    
    ############### HIEU ADD TO CONTROL ORDER COLUMNs ##############
    cols = pd.read_csv('{}/{}.csv'.format(path_to_example_samples, feat), index_col = [0]).columns.tolist()
    featuredf[feat] = featuredf[feat][cols] 
    ############### HIEU ADD TO CONTROL ORDER COLUMNs ##############

def predict(input_featuredf, path_to_save_output, path_to_model_files, nmf_flen_cancer_signal, output_filename = "output"):
    all_samples = input_featuredf["FLEN"]["SampleID"].unique()
    output = dict()
    
    #####---------------------------------------------------------------------------#####
    ##### FLEN SCORE
    #####---------------------------------------------------------------------------#####
    inputdf = input_featuredf["FLEN"].set_index("SampleID").T.reset_index()
    inputdf.columns = ["size"] + list(inputdf.columns)[1:]
    healthy_flendf = pd.read_csv(os.path.join(path_to_model_files, "healthy_flendf.csv"))
    healthy_flendf.columns = ["size", "Healthy"]
    
    healthy_flendf["size"] = healthy_flendf["size"].astype(str)
    inputdf["size"] = inputdf["size"].astype(str)
    inputdf = inputdf.merge(healthy_flendf, right_on = "size", left_on = "size")
    for sampleid in [item for item in inputdf.columns if item not in ["size", "Healthy"]]:
        inputdf[sampleid] = abs(inputdf[sampleid] - inputdf["Healthy"])
    
    flen_scoredf = inputdf.drop(["Healthy", "size"], axis = 1).sum().reset_index()
    flen_scoredf.columns = ["SampleID", "FLEN"]
    output["FLEN"] = flen_scoredf
    
    #####---------------------------------------------------------------------------#####
    ##### EM
    #####---------------------------------------------------------------------------#####
    input_feature = "EM"
    inputdf = input_featuredf[input_feature].set_index("SampleID").T
    
    em_shannondf = pd.DataFrame(data = inputdf.columns, columns = ["SampleID"])
    def calculate_em_shannon(x, inputdf):
        tmpdf = inputdf[x].values
        shannon = -np.sum([item * np.log2(item) for item in tmpdf])/256
        return(shannon)
    em_shannondf["em_shannon"] = em_shannondf["SampleID"].apply(lambda x: calculate_em_shannon(x, inputdf))
    em_shannondf = em_shannondf[["SampleID", "em_shannon"]]
    em_shannondf.columns = ["SampleID", "EM"]
    output["EM"] = em_shannondf
    
    #####---------------------------------------------------------------------------#####
    ##### NUCLEOSOME SCORE
    #####---------------------------------------------------------------------------#####
    inputdf = input_featuredf["NUCLEOSOME"].set_index("SampleID").T.reset_index()
    inputdf.columns = ["size"] + list(inputdf.columns)[1:]
    healthy_nucdf = pd.read_csv(os.path.join(path_to_model_files, "healthy_nucdf.csv"))
    healthy_nucdf.columns = ["size", "Healthy"]
    
    healthy_nucdf["size"] = healthy_nucdf["size"].astype(str)
    inputdf["size"] = inputdf["size"].astype(str)
    inputdf = inputdf.merge(healthy_nucdf, right_on = "size", left_on = "size")
    for sampleid in [item for item in inputdf.columns if item not in ["size", "Healthy"]]:
        inputdf[sampleid] = abs(inputdf[sampleid] - inputdf["Healthy"])
    
    nuc_scoredf = inputdf.drop(["Healthy", "size"], axis = 1).sum().reset_index()
    nuc_scoredf.columns = ["SampleID", "NUCLEOSOME"]
    output["NUCLEOSOME"] = nuc_scoredf
    
    #####---------------------------------------------------------------------------#####
    ##### OPTIMAL TRANSPORT FOR FLEN
    #####---------------------------------------------------------------------------#####
    def calculate_ot_distance_to_healthy_flen(sample1, bary_l2, inputdf):
        n = len(range(50, 351))
        x = np.arange(n, dtype=np.float64)
        M = ot.dist(x.reshape((n, 1)), x.reshape((n, 1)), 'euclidean')
        M /= M.max() * 0.1
        a = inputdf[sample1].values
        a = np.array(a)
        b = bary_l2
        d_emd = ot.emd2(a, b, M)  # direct computation of OT loss
        return d_emd
    
    input_feature = "FLEN"
    inputdf = input_featuredf[input_feature].set_index("SampleID").T
    with open(os.path.join(path_to_model_files, "bary_l2_flen.npy"), 'rb') as f:
        bary_l2 = np.load(f)
    
    ot_flendf = pd.DataFrame(data = inputdf.columns, columns = ["SampleID"])
    ot_flendf["OT_FLEN"] = ot_flendf["SampleID"].apply(lambda x: calculate_ot_distance_to_healthy_flen(x, bary_l2, inputdf))
    ot_flendf = ot_flendf[["SampleID", "OT_FLEN"]]
    output["OT_FLEN"] = ot_flendf
    
    #####---------------------------------------------------------------------------#####
    ##### OPTIMAL TRANSPORT FOR NUC
    #####---------------------------------------------------------------------------#####
    def calculate_ot_distance_to_healthy_nuc(sample1, bary_l2, inputdf):
        n = len(range(-300, 301))
        x = np.arange(n, dtype=np.float64)
        M = ot.dist(x.reshape((n, 1)), x.reshape((n, 1)), 'euclidean')
        M /= M.max() * 0.1
        a = inputdf[sample1].values
        a = np.array(a)
        b = bary_l2
        d_emd = ot.emd2(a, b, M)  # direct computation of OT loss
        return d_emd
        
    input_feature = "NUCLEOSOME"
    inputdf = input_featuredf[input_feature].set_index("SampleID").T
    with open(os.path.join(path_to_model_files, "bary_l2_nuc.npy"), 'rb') as f:
        bary_l2 = np.load(f)
    
    ot_nucdf = pd.DataFrame(data = inputdf.columns, columns = ["SampleID"])
    ot_nucdf["OT_NUC"] = ot_nucdf["SampleID"].apply(lambda x: calculate_ot_distance_to_healthy_nuc(x, bary_l2, inputdf))
    ot_nucdf.columns = ["SampleID", "OT_NUCLEOSOME"]
    output["OT_NUC"] = ot_nucdf
    
    #####---------------------------------------------------------------------------#####
    ##### NMF FLEN
    #####---------------------------------------------------------------------------#####
    X = input_featuredf["FLEN"].set_index("SampleID")
    filename = os.path.join(path_to_model_files, 'NMF_flen.sav')
    model = pickle.load(open(filename, 'rb'))
    W = model.transform(X)
    nmfdf = pd.DataFrame(data = W, columns = ["V1", "V2"])
    nmfdf["SampleID"] = list(X.index)
    nmfdf["V1_scale"] = nmfdf[["V1", "V2"]].apply(lambda x: x[0]/sum(x), axis = 1)
    nmfdf["V2_scale"] = nmfdf[["V1", "V2"]].apply(lambda x: x[1]/sum(x), axis = 1)
    nmfdf = nmfdf[["SampleID", "V{}_scale".format(nmf_flen_cancer_signal)]]
    nmfdf.columns = ["SampleID", "NMF_FLEN"]
    output["NMF_FLEN"] = nmfdf

    
    #####---------------------------------------------------------------------------#####
    ##### NMF NUCLEOSOME
    #####---------------------------------------------------------------------------#####
    X = input_featuredf["NUCLEOSOME"].set_index("SampleID")
    filename = os.path.join(path_to_model_files, 'NMF_nuc.sav')
    model = pickle.load(open(filename, 'rb'))
    W = model.transform(X)
    nmfdf = pd.DataFrame(data = W, columns = ["V1", "V2"])
    nmfdf["SampleID"] = list(X.index)
    nmfdf["V1_scale"] = nmfdf[["V1", "V2"]].apply(lambda x: x[0]/sum(x), axis = 1)
    nmfdf["V2_scale"] = nmfdf[["V1", "V2"]].apply(lambda x: x[1]/sum(x), axis = 1)
    nmfdf = nmfdf[["SampleID", "V1_scale"]]
    nmfdf.columns = ["SampleID", "NMF_NUCLEOSOME"]
    output["NMF_NUCLEOSOME"] = nmfdf
    
    #####---------------------------------------------------------------------------#####
    ##### OUTPUT
    #####---------------------------------------------------------------------------#####
    outputdf = pd.DataFrame(data = all_samples, columns = ["SampleID"])
    for k in output.keys():
        outputdf = outputdf.merge(output[k], right_on = "SampleID", left_on = "SampleID")
    outputdf = outputdf.merge(input_featuredf["IchorCNA"], right_on = "SampleID", left_on = "SampleID")
    feature_orders = ['SampleID', 'ichorCNA', 'FLEN', 'EM', 'NUCLEOSOME', 'OT_FLEN', 'OT_NUCLEOSOME', 'NMF_FLEN', 'NMF_NUCLEOSOME']
    outputdf = outputdf[feature_orders]
    
    for feature in feature_orders[1:]:
        outputdf[feature] = round(outputdf[feature], 6)

    outputdf.to_excel(os.path.join(path_to_save_output, "{}.xlsx".format(output_filename)))
    return(outputdf)

nmf_flen_cancer_signal = 2
predict(featuredf, path_to_save_output, path_to_model_files, nmf_flen_cancer_signal)