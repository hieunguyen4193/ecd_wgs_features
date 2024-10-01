import argparse
import pandas as pd 
import os
from feature_class import *
##### Feature class definition
def main():
    parser = argparse.ArgumentParser(description='Generate an image matrix.')
    parser.add_argument('--input', type=str, required=True, help='Path to the pre-processed frag.tsv files from 01 and 02')
    parser.add_argument('--output', type=str, required=False, help='Path to save output feature csv files')
    parser.add_argument('--motif_order_path', type=str, required=True, help='Path to the motif order file')
    parser.add_argument('--feature_version', type=str, required=True, help='Name of the feature version')
    
    args = parser.parse_args()
    input_tsv = args.input
    motif_order_path = args.motif_order_path 
    outputdir = args.output
    feature_version = args.feature_version
    
    output_obj = WGS_GW_Image_features(input_tsv = input_tsv,
                             motif_order_path = motif_order_path,
                             outputdir = outputdir,
                             feature_version = feature_version)
    
    ##### generate GW features and save features to output dir
    output_obj.generate_flen_feature()
    output_obj.generate_em_feature()
    # output_obj.generate_nuc_feature()
    output_obj.generate_nuc_feature_1()

    ##### generaet IMAGES feature and save features to output dir
    output_obj.generate_EM_flen_feature()
    output_obj.generate_forwardNUC_flen_feature()
    output_obj.generate_reverseNUC_flen_feature()
    output_obj.generate_EM_pairs_all_flen()
    output_obj.generate_EM_pairs_short_flen()
    output_obj.generate_EM_pairs_long_flen()
    output_obj.generate_EM_forwardNUC()
    output_obj.generate_EM_reverseNUC()

if __name__ == '__main__':
    main()