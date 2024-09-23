import os, sys, argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

class FeatureTransform:
    def __init__(self, 
                 motif_order_path: str = "/mnt/NAS_PROJECT/vol_ECDteam/DATA_TANPHAM/WGS_image_feature/motif_order.csv"):
        # self.MASTER_FOLDER = master_folder_path

        self.motif_order_path = motif_order_path
        self.motif_order = pd.read_csv(motif_order_path)["motif_order"].values
        self.all_4bp_motifs = [
            "{}{}{}{}".format(i,j,k,l) 
            for i in ["A", "T", "G", "C"] 
            for j in ["A", "T", "G", "C"] 
            for k in ["A", "T", "G", "C"] 
            for l in ["A", "T", "G", "C"]
        ]
    
    def EM_flen_feature(self,
                        feature_df: pd.DataFrame,
                        save_path: str):
        # Forward EM
        feature_df["forward_EM"] = feature_df["forward_EM"].apply(lambda x: x.upper())
        # Reverse EM
        feature_df["reverse_EM"] = feature_df["reverse_EM"].apply(lambda x: x.upper())

        # Flen
        feature_df = feature_df[(feature_df["flen"] >= 70) & (feature_df["flen"] <= 280)]
        
        ##### Generate EM - FLEN dataframe
        # Forward EM - flen
        forward_em_flen = feature_df[["forward_EM", "flen"]].copy()
        forward_em_flen.columns = ["EM", "flen"]

        # Reverse EM - flen
        reverse_em_flen = feature_df[["reverse_EM", "flen"]].copy()
        reverse_em_flen.columns = ["EM", "flen"]

        # EM - flen df
        em_flen_df = pd.concat([forward_em_flen, reverse_em_flen], axis = 0)
        em_flen_df = em_flen_df[~em_flen_df["EM"].str.contains("N")]
        countdf = em_flen_df.reset_index() \
                            .groupby(["EM", "flen"])["index"] \
                            .count() \
                            .reset_index() \
                            .pivot_table(index='flen', 
                                         columns='EM', 
                                         values='index', 
                                         fill_value=0)
        
        ##### fill values so that the output matrix always 50:350 x 256
        flen_range_df = pd.DataFrame(
            {
                'flen': range(70, 281)
            }
        )
        countdf = pd.merge(flen_range_df, 
                           countdf, 
                           on='flen', 
                           how='outer')
        countdf.fillna(0, inplace=True)
        
        countdf = countdf.set_index("flen")
        missing_motifs = [item for item in countdf.columns if item not in self.all_4bp_motifs]

        if len(missing_motifs) != 0:
            for motif in missing_motifs:
                countdf[motif] = 0
        
        scaled_countdf = countdf/countdf.sum().sum()
        scaled_countdf = scaled_countdf[self.motif_order]
        scaled_countdf.to_csv(save_path, index=False)
        
    def nucleosome_forward_flen(self,
                                feature_df: pd.DataFrame,
                                save_path: str):
        nucdf_forward = feature_df[["flen", "forward_NUC"]].copy()
        nucdf_forward.columns = ["flen", "nuc_dist"]
        nucdf_forward = nucdf_forward[
            (nucdf_forward["nuc_dist"] <= 300) & (nucdf_forward["nuc_dist"] >= -300)
        ]

        nuc_countdf = nucdf_forward.reset_index() \
                           .groupby(["nuc_dist", "flen"])["index"] \
                           .count() \
                           .reset_index() \
                           .pivot_table(index='flen', 
                                     columns='nuc_dist', 
                                     values='index', 
                                     fill_value=0)
        flen_range_df = pd.DataFrame(
            {
                'flen': range(70, 281)
            }
        )
        
        nuc_countdf = pd.merge(flen_range_df,
                               nuc_countdf,
                               on='flen',
                               how='outer')
        nuc_countdf.fillna(0, 
                           inplace=True)

        nuc_countdf = nuc_countdf[(nuc_countdf['flen'] >= 70) & (nuc_countdf['flen'] <= 280)]

        nuc_countdf = nuc_countdf.set_index("flen")
        nuc_countdf = nuc_countdf/nuc_countdf.sum().sum()
        
        assert nuc_countdf.shape[0] == 211, f"[NUC forward - flen] flen failed!"
        assert nuc_countdf.shape[1] == 601, f"[NUC forward - flen] NUC failed"
        
        nuc_countdf.to_csv(save_path, index=False)
        
    def nucleosome_reverse_flen(self,
                                feature_df: pd.DataFrame,
                                save_path: str):
        nucdf_reverse = feature_df[["flen", "reverse_NUC"]].copy()
        nucdf_reverse.columns = ["flen", "nuc_dist"]
        nucdf_reverse = nucdf_reverse[
            (nucdf_reverse["nuc_dist"] <= 300) & (nucdf_reverse["nuc_dist"] >= -300)
        ]

        nuc_countdf = nucdf_reverse.reset_index() \
                           .groupby(["nuc_dist", "flen"])["index"] \
                           .count() \
                           .reset_index() \
                           .pivot_table(index='flen', 
                                     columns='nuc_dist', 
                                     values='index', 
                                     fill_value=0)
        flen_range_df = pd.DataFrame(
            {
                'flen': range(70, 281)
            }
        )
        nuc_countdf = pd.merge(flen_range_df,
                               nuc_countdf,
                               on='flen',
                               how='outer')
        nuc_countdf.fillna(0, 
                           inplace=True)

        nuc_countdf = nuc_countdf[(nuc_countdf['flen'] >= 70) & (nuc_countdf['flen'] <= 280)]

        nuc_countdf = nuc_countdf.set_index("flen")
        nuc_countdf = nuc_countdf/nuc_countdf.sum().sum()
        
        assert nuc_countdf.shape[0] == 211, f"[NUC reverse - flen] flen failed!"
        assert nuc_countdf.shape[1] == 601, f"[NUC reverse - flen] NUC failed"
        
        nuc_countdf.to_csv(save_path, index=False)
    
    def nucleosome_flen_feature(self,
                                feature_df: pd.DataFrame,
                                save_path: str):
        nucdf_reverse = feature_df[["flen", "reverse_NUC"]].copy()
        nucdf_reverse.columns = ["flen", "nuc_dist"]

        nucdf_forward = feature_df[["flen", "forward_NUC"]].copy()
        nucdf_forward.columns = ["flen", "nuc_dist"]
        nucdf = pd.concat([nucdf_forward, nucdf_reverse], axis = 0)
        nucdf = nucdf[
            (nucdf["nuc_dist"] <= 300) & (nucdf["nuc_dist"] >= -300)
        ]

        nuc_countdf = nucdf.reset_index() \
                           .groupby(["nuc_dist", "flen"])["index"] \
                           .count() \
                           .reset_index() \
                           .pivot_table(index='flen', 
                                     columns='nuc_dist', 
                                     values='index', 
                                     fill_value=0)
        
        flen_range_df = pd.DataFrame(
            {
                'flen': range(70, 281)
            }
        )
        
        nuc_countdf = pd.merge(flen_range_df,
                               nuc_countdf,
                               on='flen',
                               how='outer')

        nuc_countdf.fillna(0, 
                           inplace=True)

        nuc_countdf = nuc_countdf.set_index("flen")
        nuc_countdf = nuc_countdf/nuc_countdf.sum().sum()
        nuc_countdf = nuc_countdf[self.motif_order]
        nuc_countdf.to_csv(save_path, index=False)

    def motif_pairs_all_flen(self,
                             feature_df: pd.DataFrame,
                             save_path: str):

        count_pair_EM = feature_df[(~ feature_df["reverse_EM"].str.contains("N")) & (~feature_df["forward_EM"].str.contains("N"))] \
                        .groupby(["forward_EM", "reverse_EM"])["readID"] \
                        .count() \
                        .reset_index() \
                        .reset_index() \
                        .pivot_table(index='forward_EM', 
                                     columns='reverse_EM', 
                                     values='readID', 
                                     fill_value=0)

        count_pair_EM = count_pair_EM/count_pair_EM.sum().sum()
        
        for motif in self.motif_order:
            if motif not in count_pair_EM.keys():
                count_pair_EM[motif] = [0 for _ in range(len(count_pair_EM))]
        
        for motif in self.motif_order:
            if motif not in count_pair_EM.index:
                count_pair_EM.loc[len(count_pair_EM.index)] = [0 for _ in range(len(count_pair_EM.columns))]
                count_pair_EM = count_pair_EM.rename(index={len(count_pair_EM) - 1: motif})
        
        count_pair_EM = count_pair_EM.loc[self.motif_order][self.motif_order]
        
        assert count_pair_EM.shape[0] == 256, f"Motif pairs all flen failed!"
        assert count_pair_EM.shape[1] == 256, f"Motif pairs all flen failed!"
        
        count_pair_EM.to_csv(save_path, index=False)
        
    def motif_pairs_short_flen(self,
                               feature_df: pd.DataFrame,
                               save_path: str):
        
        feature_df_short = feature_df[feature_df["flen"] <= 150]

        count_pair_EM = feature_df_short[(~ feature_df_short["reverse_EM"].str.contains("N")) & (~feature_df_short["forward_EM"].str.contains("N"))] \
                        .groupby(["forward_EM", "reverse_EM"])["readID"] \
                        .count() \
                        .reset_index() \
                        .reset_index() \
                        .pivot_table(index='forward_EM', 
                                    columns='reverse_EM', 
                                    values='readID', 
                                    fill_value=0)

        count_pair_EM = count_pair_EM/count_pair_EM.sum().sum()

        for motif in self.motif_order:
            if motif not in count_pair_EM.keys():
                count_pair_EM[motif] = [0 for _ in range(len(count_pair_EM))]
        
        for motif in self.motif_order:
            if motif not in count_pair_EM.index:
                count_pair_EM.loc[len(count_pair_EM.index)] = [0 for _ in range(len(count_pair_EM.columns))]
                count_pair_EM = count_pair_EM.rename(index={len(count_pair_EM) - 1: motif})
        
        count_pair_EM = count_pair_EM.loc[self.motif_order][self.motif_order]

        assert count_pair_EM.shape[0] == 256, f"Motif pairs short flen failed!"
        assert count_pair_EM.shape[1] == 256, f"Motif pairs short flen failed!"

        count_pair_EM.to_csv(save_path, index=False)
        
    def motif_pairs_long_flen(self,
                              feature_df: pd.DataFrame,
                              save_path: str):

        feature_df_long = feature_df[feature_df["flen"] > 150]
        count_pair_EM = feature_df_long[(~ feature_df_long["reverse_EM"].str.contains("N")) & (~feature_df_long["forward_EM"].str.contains("N"))] \
                    .groupby(["forward_EM", "reverse_EM"])["readID"] \
                    .count() \
                    .reset_index() \
                    .reset_index() \
                        .pivot_table(index='forward_EM', 
                                     columns='reverse_EM', 
                                     values='readID', 
                                     fill_value=0)

        count_pair_EM = count_pair_EM/count_pair_EM.sum().sum()

        for motif in self.motif_order:
            if motif not in count_pair_EM.keys():
                count_pair_EM[motif] = [0 for _ in range(len(count_pair_EM))]
        
        for motif in self.motif_order:
            if motif not in count_pair_EM.index:
                count_pair_EM.loc[len(count_pair_EM.index)] = [0 for _ in range(len(count_pair_EM.columns))]
                count_pair_EM = count_pair_EM.rename(index={len(count_pair_EM) - 1: motif})
        
        count_pair_EM = count_pair_EM.loc[self.motif_order][self.motif_order]

        assert count_pair_EM.shape[0] == 256, f"Motif pairs long flen failed!"
        assert count_pair_EM.shape[1] == 256, f"Motif pairs long flen failed!"
        
        count_pair_EM.to_csv(save_path, index=False)
    
    def EM_forward_nucleosome(self,
                              feature_df: pd.DataFrame,
                              save_path: str):

        ##### generate EM - forward_NUC dataframe
        # Forward EM - Forward nucleosome distance
        forward_em_forward_NUC = feature_df[["forward_EM", "forward_NUC"]].copy()
        forward_em_forward_NUC.columns = ["EM", "forward_NUC"]
        
        # Reverse EM - forward nucleosome distance
        reverse_em_forward_NUC = feature_df[["reverse_EM", "forward_NUC"]].copy()
        reverse_em_forward_NUC.columns = ["EM", "forward_NUC"]

        em_forward_NUC_df = pd.concat([forward_em_forward_NUC, reverse_em_forward_NUC], axis = 0)
        em_forward_NUC_df = em_forward_NUC_df[~em_forward_NUC_df["EM"].str.contains("N")]
        em_forward_NUC_df = em_forward_NUC_df[(em_forward_NUC_df["forward_NUC"] >= -300) & (em_forward_NUC_df["forward_NUC"] <= 300)]
        countdf = em_forward_NUC_df.reset_index() \
                                   .groupby(["EM", "forward_NUC"])["index"] \
                                   .count() \
                                   .reset_index() \
                                   .pivot_table(index='forward_NUC', 
                                                columns='EM', 
                                                values='index', 
                                                fill_value=0)

        ##### fill values so that the output matrix always 50:350 x 256
        forward_NUC_range_df = pd.DataFrame(
            {
                'forward_NUC': range(70, 281)
            }
        )

        countdf = pd.merge(forward_NUC_range_df, countdf, on='forward_NUC', how='outer')
        countdf.fillna(0, inplace=True)

        countdf = countdf.set_index("forward_NUC")
        missing_motifs = [item for item in countdf.columns if item not in self.all_4bp_motifs]

        if len(missing_motifs) != 0:
            for motif in missing_motifs:
                countdf[motif] = 0
                
        countdf = countdf/countdf.sum().sum()
        countdf = countdf[self.motif_order]
        countdf.to_csv(save_path, index=False)

    def EM_reverse_nucleosome(self,
                              feature_df: pd.DataFrame,
                              save_path: str):
        
        ##### generate EM - reverse_NUC dataframe
        reverse_EM_reverse_NUC = feature_df[["reverse_EM", "reverse_NUC"]].copy()
        reverse_EM_reverse_NUC.columns = ["EM", "reverse_NUC"]
        reverse_em_reverse_NUC = feature_df[["reverse_EM", "reverse_NUC"]].copy()
        reverse_em_reverse_NUC.columns = ["EM", "reverse_NUC"]
        em_reverse_NUC_df = pd.concat([reverse_EM_reverse_NUC, reverse_em_reverse_NUC], axis = 0)
        em_reverse_NUC_df = em_reverse_NUC_df[~em_reverse_NUC_df["EM"].str.contains("N")]
        
        em_reverse_NUC_df = em_reverse_NUC_df[
            (em_reverse_NUC_df["reverse_NUC"] >= -300) 
            & (em_reverse_NUC_df["reverse_NUC"] <= 300)
        ]
        
        countdf = em_reverse_NUC_df.reset_index() \
                                    .groupby(["EM", "reverse_NUC"])["index"] \
                                    .count() \
                                    .reset_index() \
                                    .pivot_table(index='reverse_NUC', 
                                                 columns='EM', 
                                                 values='index', 
                                                 fill_value=0)

        ##### fill values so that the output matrix always 50:350 x 256
        reverse_NUC_range_df = pd.DataFrame(
            {
                'reverse_NUC': range(70, 281)
            }
        )
        countdf = pd.merge(reverse_NUC_range_df, 
                           countdf, 
                           on='reverse_NUC', 
                           how='outer')
        countdf.fillna(0, 
                       inplace=True)

        
        countdf = countdf.set_index("reverse_NUC")
        missing_motifs = [item for item in countdf.columns if item not in self.all_4bp_motifs]

        if len(missing_motifs) != 0:
            for motif in missing_motifs:
                countdf[motif] = 0
        countdf = countdf/countdf.sum().sum()
        countdf = countdf[self.motif_order]
        countdf.to_csv(save_path, index=False)
    
    def transform_feature_to_image(self, 
                                    tsv_file_path: str = None,
                                    feature_type: str = None,
                                    save_path: str = ""):
        """
        feature_type (str) -- Type of features 
        Currently support for ['EM_flen', 'nucleosome_flen', 'motif_pairs_all_flen', 'motif_pairs_long_flen', 'motif_pairs_short_flen', 'EM_forward_nucleosome', 'EM_reverse_nucleosome']
        """
        assert tsv_file_path != None
        # input_data = os.path.join(self.MASTER_FOLDER, tsv_file_path)
        input_data = tsv_file_path
        
        print(input_data)
        
        assert feature_type != None
        feature_df = pd.read_csv(input_data, sep = "\t", header = None)
        feature_df = feature_df[[0, 1, 2, 3, 4, 8, 9, 10, 11, 12]]
        feature_df.columns = ["readID", "chrom", "start", "cigar", "flen", "readID_extra", "forward_NUC", "reverse_NUC", "forward_EM", "reverse_EM"]

        # DropNA
        feature_df = feature_df.dropna()
        
        if feature_type == 'EM_flen':
            self.EM_flen_feature(feature_df,
                                 save_path)
        elif feature_type == "nucleosome_flen":
            self.nucleosome_flen_feature(feature_df,
                                         save_path)
        elif feature_type == "motif_pairs_all_flen":
            self.motif_pairs_all_flen(feature_df,
                                      save_path)
        elif feature_type == "motif_pairs_long_flen":
            self.motif_pairs_long_flen(feature_df,
                                       save_path)
        elif feature_type == "motif_pairs_short_flen":
            self.motif_pairs_short_flen(feature_df,
                                        save_path)
        elif feature_type == "EM_forward_nucleosome":
            self.EM_forward_nucleosome(feature_df,
                                       save_path)
        elif feature_type == "EM_reverse_nucleosome":
            self.EM_reverse_nucleosome(feature_df,
                                       save_path)
        elif feature_type == "nucleosome_forward_flen":
            self.nucleosome_forward_flen(feature_df,
                                         save_path)
        elif feature_type == "nucleosome_reverse_flen":
            self.nucleosome_reverse_flen(feature_df,
                                         save_path)
        else:
            raise NotImplementedError(f"Feature {feature_type} is not yet supported. Currently support for ['EM_flen', 'nucleosome_flen', 'motif_pairs_all_flen', 'motif_pairs_long_flen', 'motif_pairs_short_flen', 'EM_forward_nucleosome', 'EM_reverse_nucleosome']")

def transform(args):
    feature_to_image = FeatureTransform(motif_order_path=args.motif_order_path)
    feature_to_image.transform_feature_to_image(tsv_file_path=args.tsv_file_path,
                                            feature_type=args.feature_type,
                                            save_path=args.save_path)

def parse_arguments(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_id", type=str, action="store")
    parser.add_argument("--motif_order_path", type=str, action="store")
    parser.add_argument("--tsv_file_path", type=str, action="store")
    parser.add_argument("--feature_type", type=str, action="store")
    parser.add_argument("--save_path", type=str, action="store")
    return parser.parse_args(argv)


if __name__ == "__main__":
    transform(parse_arguments(sys.argv[1:]))
