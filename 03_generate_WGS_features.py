import argparse
import pandas as pd 
import os
from tqdm import tqdm
from feature_class import *

##### Feature class definition
def main():
    parser = argparse.ArgumentParser(description='Generate an image matrix.')
    parser.add_argument('--input', type=str, required=True, help='Path to the pre-processed frag.tsv files from 01 and 02')
    parser.add_argument('--output', type=str, required=False, help='Path to save output feature csv files')
    parser.add_argument('--motif_order_path', type=str, required=True, help='Path to the motif order file')
    parser.add_argument('--feature_version', type=str, required=True, help='Name of the feature version')
    parser.add_argument('--old_nuc', type=str, required=True, help='Path to the old nucleosome file')
    parser.add_argument('--generate_feature', type=str, required=True, help='Path to the old nucleosome file')
    
    args = parser.parse_args()
    input_tsv = args.input
    motif_order_path = args.motif_order_path 
    outputdir = args.output
    feature_version = args.feature_version
    path_to_old_nuc = args.old_nuc
    generate_feature = args.generate_feature
    
    output_obj = WGS_GW_Image_features(input_tsv = input_tsv,
                             motif_order_path = motif_order_path,
                             outputdir = outputdir,
                             path_to_old_nuc = path_to_old_nuc,
                             feature_version = feature_version)
    
    ##### generate GW features and save features to output dir
    if generate_feature == "all":
        output_obj.generate_flen_feature()
        output_obj.generate_em_feature()
        if path_to_old_nuc == "none":
            output_obj.generate_nuc_feature()
        else:
            output_obj.generate_nuc_feature_1()
        # output_obj.generate_ndr_feature(binary_or_TOO="TOO")
        # output_obj.generate_ndr_feature(binary_or_TOO="binary")

        ##### generate IMAGES feature and save features to output dir
        output_obj.generate_EM_flen_feature()
        output_obj.generate_forwardNUC_flen_feature()
        output_obj.generate_reverseNUC_flen_feature()
        output_obj.generate_EM_pairs_all_flen()
        output_obj.generate_EM_pairs_short_flen()
        output_obj.generate_EM_pairs_long_flen()
        output_obj.generate_allEM_forwardNUC()
        output_obj.generate_allEM_reverseNUC()
        output_obj.generate_reverseEM_reverseNUC()
        output_obj.generate_forwardEM_forwardNUC()

        # for binary_or_TOO in ["TOO", "binary"]:
            # output_obj.generate_forwardNDR_flen_feature(binary_or_TOO = binary_or_TOO)
            # output_obj.generate_reverseNDR_flen_feature(binary_or_TOO = binary_or_TOO)
            # output_obj.generate_allEM_forwardNDR(binary_or_TOO = binary_or_TOO)
            # output_obj.generate_allEM_reverseNDR(binary_or_TOO = binary_or_TOO)
            # output_obj.generate_reverseEM_reverseNDR(binary_or_TOO = binary_or_TOO)
            # output_obj.generate_forwardEM_forwardNDR(binary_or_TOO = binary_or_TOO)
    
    elif generate_feature == "image_only":
        output_obj.generate_EM_flen_feature()
        output_obj.generate_forwardNUC_flen_feature()
        output_obj.generate_reverseNUC_flen_feature()
        output_obj.generate_EM_pairs_all_flen()
        output_obj.generate_EM_pairs_short_flen()
        output_obj.generate_EM_pairs_long_flen()
        output_obj.generate_allEM_forwardNUC()
        output_obj.generate_allEM_reverseNUC()
        output_obj.generate_reverseEM_reverseNUC()
        output_obj.generate_forwardEM_forwardNUC()

        # for binary_or_TOO in ["TOO", "binary"]:
        #     output_obj.generate_forwardNDR_flen_feature(binary_or_TOO = binary_or_TOO)
        #     output_obj.generate_reverseNDR_flen_feature(binary_or_TOO = binary_or_TOO)
        #     output_obj.generate_allEM_forwardNDR(binary_or_TOO = binary_or_TOO)
        #     output_obj.generate_allEM_reverseNDR(binary_or_TOO = binary_or_TOO)
        #     output_obj.generate_reverseEM_reverseNDR(binary_or_TOO = binary_or_TOO)
        #     output_obj.generate_forwardEM_forwardNDR(binary_or_TOO = binary_or_TOO)    

    elif generate_feature == "GW_only":
        output_obj.generate_flen_feature()
        output_obj.generate_em_feature()
        if path_to_old_nuc == "none":
            output_obj.generate_nuc_feature()
        else:
            output_obj.generate_nuc_feature_1()
        # output_obj.generate_ndr_feature(binary_or_TOO="TOO")
        # output_obj.generate_ndr_feature(binary_or_TOO="binary")
    # elif generate_feature == "new_feature_v0.1.9_only":
    #     output_obj.generate_ndr_feature(binary_or_TOO="TOO")
    #     output_obj.generate_ndr_feature(binary_or_TOO="binary")
        
        output_obj.generate_allEM_forwardNUC()
        output_obj.generate_allEM_reverseNUC()
        output_obj.generate_reverseEM_reverseNUC()
        output_obj.generate_forwardEM_forwardNUC()
        # for binary_or_TOO in ["TOO", "binary"]:
        #     output_obj.generate_forwardNDR_flen_feature(binary_or_TOO = binary_or_TOO)
        #     output_obj.generate_reverseNDR_flen_feature(binary_or_TOO = binary_or_TOO)
        #     output_obj.generate_allEM_forwardNDR(binary_or_TOO = binary_or_TOO)
        #     output_obj.generate_allEM_reverseNDR(binary_or_TOO = binary_or_TOO)
        #     output_obj.generate_reverseEM_reverseNDR(binary_or_TOO = binary_or_TOO)
        #     output_obj.generate_forwardEM_forwardNDR(binary_or_TOO = binary_or_TOO)    
    else:
        raise ValueError("Please specify the correct feature type to generate")

if __name__ == '__main__':
    main()