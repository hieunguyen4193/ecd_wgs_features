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
        # if no metadata is provided, generate drift measures for all samples in the batch
        flendf = obj.generate_flen_matrix()
        emdf = obj.generate_em_matrix()
        nucdf = obj.generate_nuc_matrix()
        
    median_flendf = flendf.median(axis=1)
    median_emdf = emdf.median(axis=1)
    median_nucdf = nucdf.median(axis=1)

    ##### keep this path default, the feature_drift_reference always goes with the repo
    flen_barycenter = pd.read_csv(f"feature_drift_reference/OT/{reference_version}/flen_barycenter.csv")
    em_barycenter = pd.read_csv(f"feature_drift_reference/OT/{reference_version}/em_barycenter.csv")
    nuc_barycenter = pd.read_csv(f"feature_drift_reference/OT/{reference_version}/nuc_barycenter.csv")

    median_ref_flendf = pd.read_csv(f"feature_drift_reference/APE/{reference_version}/median_flendf.csv")
    median_ref_emdf = pd.read_csv(f"feature_drift_reference/APE/{reference_version}/median_emdf.csv")
    median_ref_nucdf = pd.read_csv(f"feature_drift_reference/APE/{reference_version}/median_nucdf.csv")

    ##### OT dist
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

    ##### APE
    # flen has too many median 0 features, remove them before calculating APE
    ape_flendf = pd.DataFrame(data = median_flendf, columns=["median_flen"])
    ape_flendf["ref_median_flen"] = median_ref_flendf["0"].to_numpy()
    ape_flendf = ape_flendf[(ape_flendf["median_flen"] != 0) & (ape_flendf["ref_median_flen"] != 0)].reset_index()

    ape_flen = abs(ape_flendf["median_flen"].to_numpy() - ape_flendf["ref_median_flen"].to_numpy()) / ape_flendf["ref_median_flen"].to_numpy()
    ape_em = abs(median_emdf - median_ref_emdf["0"].to_numpy()) / median_ref_emdf["0"].to_numpy()
    ape_nuc = abs(median_nucdf - median_ref_nucdf["0"].to_numpy()) / median_ref_nucdf["0"].to_numpy()

    ape_flen = pd.DataFrame(
        {"feat": ape_flendf.feat.unique(),
        "ape": ape_flen
        }
    ).reset_index().drop("index", axis = 1)

    ape_em = pd.DataFrame(
        {"feat": ape_em.index,
        "ape": ape_em
        }    
    ).reset_index().drop("index", axis = 1)

    ape_nuc = pd.DataFrame(
        {"feat": ape_nuc.index,
        "ape": ape_nuc
        }
    ).reset_index().drop("index", axis = 1)

    flen_distdf.to_csv(f"{outputdir}/flen_OTdist_to_ref.csv", index=False)
    em_distdf.to_csv(f"{outputdir}/em_OTdist_to_ref.csv", index=False)
    nuc_distdf.to_csv(f"{outputdir}/nuc_OTdist_to_ref.csv", index=False)
    
    ape_flen.to_csv(f"{outputdir}/flen_APE_to_ref.csv", index=False)
    ape_em.to_csv(f"{outputdir}/em_APE_to_ref.csv", index=False)
    ape_nuc.to_csv(f"{outputdir}/nuc_APE_to_ref.csv", index=False)

if __name__ == '__main__':
    main()

