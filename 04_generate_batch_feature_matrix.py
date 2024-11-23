import argparse
import pandas as pd 
import os
from feature_class import *
##### Feature class definition
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, required=True, help='Path to batch output folder, a folder containing many sample, 1 sample per folder')
    parser.add_argument('--output', type=str, required=False, help='Path to save output feature csv files')
    
    args = parser.parse_args()
    input_batch = args.input
    outputdir = args.output
    
    batch_name = os.path.basename(input_batch)
    all_features = WGS_GW_features(path_to_feature_dir=input_batch,
                               path_to_metadata=None)

    emdf = all_features.generate_em_matrix().reset_index()
    nucdf = all_features.generate_nuc_matrix().reset_index()
    flendf = all_features.generate_flen_matrix().reset_index()
    ndrdf = all_features.generate_ndr_matrix(binary_or_TOO = "TOO").reset_index()
    ndrbdf = all_features.generate_ndr_matrix(binary_or_TOO = "binary").reset_index()
    os.system(f"mkdir -p {outputdir}/{batch_name}")
    
    emdf.to_csv(os.path.join(outputdir, batch_name, 'EM.csv'), index = None)
    flendf.to_csv(os.path.join(outputdir, batch_name, 'FLEN.csv'), index = None)
    nucdf.to_csv(os.path.join(outputdir, batch_name, 'NUCLEOSOME.csv'), index = None)
    ndrdf.to_csv(os.path.join(outputdir, batch_name, 'NDR.csv'), index = None)
    ndrbdf.to_csv(os.path.join(outputdir, batch_name, 'NDRb.csv'), index = None)

if __name__ == '__main__':
    main()