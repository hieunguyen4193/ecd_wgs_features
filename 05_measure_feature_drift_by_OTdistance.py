import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pathlib
from tqdm import tqdm
from feature_class import *
from helper_functions import *
import ot
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ot_ref_version', type=str, required=True, help='Version of the reference barycenter, use in calculating optimal transport distance drifts.')
    parser.add_argument('--input', type=str, required=False, help='Path to the feature directory of the batch')
    parser.add_argument('--metadata', type=str, required=False, help='Path to the metadata of the batch')
    parser.add_argument('--output', type=str, required=True, help='Path to the motif order file')
    args = parser.parse_args()
    path_to_feature_dir = args.input
    path_to_metadata = args.metadata
    reference_version = args.ot_ref_version
    outputdir = args.output
    os.system(f"mkdir -p {outputdir}")  
        
    obj = WGS_GW_features(path_to_feature_dir=path_to_feature_dir,
                        path_to_metadata=path_to_metadata)

    if path_to_metadata is not None:
        batch_metadata = obj.match_metadata.copy()
        control_samples = batch_metadata[batch_metadata["Label"] == "Control"]["SampleID"].unique()

        ##### for calculating feature drift, we use healthy control samples only. 
        flendf = obj.generate_flen_matrix()[control_samples]
        emdf = obj.generate_em_matrix()[control_samples]
        nucdf = obj.generate_nuc_matrix()[control_samples]
    else:
        # if no metadata is provided, calculate distance for all samples
        flendf = obj.generate_flen_matrix()
        emdf = obj.generate_em_matrix()
        nucdf = obj.generate_nuc_matrix()

    ##### keep this path default, the feature_drift_reference always goes with the repo
    flen_barycenter = pd.read_csv(f"feature_drift_reference/OT/{reference_version}/flen_barycenter.csv")
    em_barycenter = pd.read_csv(f"feature_drift_reference/OT/{reference_version}/em_barycenter.csv")
    nuc_barycenter = pd.read_csv(f"feature_drift_reference/OT/{reference_version}/nuc_barycenter.csv")

    flen_distdf = pd.DataFrame(data = flendf.columns, columns=["SampleID"])
    flen_distdf["dist_to_ref"] = flen_distdf["SampleID"].apply(lambda x: calculate_ot_distance_to_ref(x, 
                                                                                            flen_barycenter["flen_barycenter"].to_numpy(), 
                                                                                            flendf))

    em_distdf = pd.DataFrame(data = emdf.columns, columns=["SampleID"])
    em_distdf["dist_to_ref"] = em_distdf["SampleID"].apply(lambda x: calculate_ot_distance_to_ref(x, 
                                                                                            em_barycenter["em_barycenter"].to_numpy(), 
                                                                                            emdf, n = 256))

    nuc_distdf = pd.DataFrame(data = nucdf.columns, columns=["SampleID"])
    nuc_distdf["dist_to_ref"] = nuc_distdf["SampleID"].apply(lambda x: calculate_ot_distance_to_ref(x, 
                                                                                            nuc_barycenter["nuc_barycenter"].to_numpy(), 
                                                                                            nucdf, n = 601))

    flen_distdf.to_csv(f"{outputdir}/flen_dist_to_ref.csv", index=False)
    em_distdf.to_csv(f"{outputdir}/em_dist_to_ref.csv", index=False)
    nuc_distdf.to_csv(f"{outputdir}/nuc_dist_to_ref.csv", index=False)

if __name__ == '__main__':
    main()

