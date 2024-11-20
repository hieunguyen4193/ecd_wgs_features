import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pathlib
from tqdm import tqdm

class WGS_GW_features:
    def __init__(self,
                 path_to_feature_dir,
                 path_to_metadata):
        self.path_to_feature_dir = path_to_feature_dir
        self.all_flen_features = [
            item for item in pathlib.Path(self.path_to_feature_dir).glob("*/*_GWfeature_flen.csv")
        ]
        self.all_nuc_features = [
            item for item in pathlib.Path(self.path_to_feature_dir).glob("*/*_GWfeature_Nucleosome.csv")
        ]
        self.all_em_features = [
            item for item in pathlib.Path(self.path_to_feature_dir).glob("*/*_GWfeature_EM.csv")
        ]
        assert len(self.all_flen_features) == len(self.all_em_features)
        assert len(self.all_flen_features) == len(self.all_nuc_features)
        if path_to_metadata is not None:        
            self.path_to_metadata = path_to_metadata
            if ".xlsx" in self.path_to_metadata:
                self.metadata = pd.read_excel(self.path_to_metadata)
            elif ".csv" in self.path_to_metadata:
                self.metadata = pd.read_csv(self.path_to_metadata)
            elif ".tsv" in self.path_to_metadata:
                self.metadata = pd.read_csv(self.path_to_metadata, sep="\t")  
                
            tmp = pd.DataFrame(
                data = [file.name for file in self.all_flen_features],
                columns = ["Filename"]
            )
            tmp["SampleID"] = tmp["Filename"].apply(lambda x: x.split("_")[0] if "-" not in x else x.split("_")[0].split("-")[1])
            self.real_metadata = tmp.copy()
            self.match_metadata = tmp.merge(self.metadata, right_on = "SampleID", left_on = "SampleID")
        
    def generate_flen_matrix(self):
        maindf = pd.DataFrame(data = range(50, 351), columns = ["feat"])
        for file in tqdm(self.all_flen_features):
            tmpdf = pd.read_csv(file, index_col = [0])[["freq"]].reset_index()
            sampleid = file.name.split("_")[0]
            if "-" in sampleid:
                sampleid = sampleid.split("-")[1]
            tmpdf.columns = ["feat", sampleid]
            maindf = maindf.merge(tmpdf, right_on = "feat", left_on="feat")
        maindf = maindf.set_index("feat")
        return maindf

    def generate_nuc_matrix(self):
        maindf = pd.DataFrame(data = range(-300, 301), columns = ["feat"])
        for file in tqdm(self.all_nuc_features):
            sampleid = file.name.split("_")[0]
            if "-" in sampleid:
                sampleid = sampleid.split("-")[1]
            tmpdf = pd.read_csv(file)
            tmpdf.columns = ["feat", sampleid]
            maindf = maindf.merge(tmpdf, right_on = "feat", left_on="feat")
            
        maindf = maindf.set_index("feat")
        return maindf
    
    def generate_em_matrix(self):
        maindf = pd.DataFrame(data = ["{}{}{}{}".format(i,j,k,l) 
                                      for i in ["A", "T", "G", "C"] 
                                      for j in ["A", "T", "G", "C"] 
                                      for k in ["A", "T", "G", "C"] 
                                      for l in ["A", "T", "G", "C"]], 
                              columns = ["feat"])
        for file in tqdm(self.all_em_features):
            sampleid = file.name.split("_")[0]
            if "-" in sampleid:
                sampleid = sampleid.split("-")[1]
            tmpdf = pd.read_csv(file)[["motif", "freq"]]
            tmpdf.columns = ["feat", sampleid]
            maindf = maindf.merge(tmpdf, right_on = "feat", left_on="feat")
        maindf = maindf.set_index("feat")
        return maindf
    
################################################################################
class WGS_GW_Image_features:
    def __init__(self,
                 input_tsv,
                 motif_order_path,
                 outputdir,
                 path_to_old_nuc = "none",
                 feature_version = "20241001",
                 use_softmask = False):
        self.input_tsv = input_tsv
        self.sampleid = input_tsv.split("/")[-1].split(".")[0]
        print("reading in the input frag.tsv data")
        if feature_version == "20241001": 
            self.maindf = pd.read_csv(input_tsv, sep = "\t", header = None)
            self.maindf.columns = ["chr", "start", "end", "flen", "readID", "QC", "forward_NUC", "reverse_NUC", "forward_EM", "reverse_EM", "forward_NDR", "reverse_NDR"]
        else:
            # use only for reading the first version of *.final_output.tsv file for GW-Image features.
            self.maindf = pd.read_csv(input_tsv, sep = "\t", header = None)
            self.maindf = self.maindf[[0, 1, 2, 3, 4, 8, 9, 10, 11, 12]]
            self.maindf.columns = ["readID", "chr", "start", "cigar", "flen", "readID_extra", "forward_NUC", "reverse_NUC", "forward_EM", "reverse_EM", "forward_NDR", "reverse_NDR"]
        
        if use_softmask:
            self.maindf["forward_EM"] = self.maindf["forward_EM"].apply(lambda x: x.upper())
            self.maindf["reverse_EM"] = self.maindf["reverse_EM"].apply(lambda x: x.upper())

        self.motif_order_path = motif_order_path
        self.motif_order = pd.read_csv(motif_order_path)["motif_order"].values
        self.all_4bp_motifs = [
            "{}{}{}{}".format(i,j,k,l) 
            for i in ["A", "T", "G", "C"] 
            for j in ["A", "T", "G", "C"] 
            for k in ["A", "T", "G", "C"] 
            for l in ["A", "T", "G", "C"]
        ]
        self.maindf_filter_chr = self.maindf[(self.maindf["chr"].isin([f"chr{i}" for i in range(1, 23)])) & (self.maindf["flen"] > 0)]
        # self.outputdir = os.path.join(outputdir, self.sampleid)
        self.outputdir = outputdir
        os.system(f"mkdir -p {self.outputdir}")
        self.feature_version = feature_version
        self.path_to_old_nuc = path_to_old_nuc
        
    #####-------------------------------------------------------------#####
    ##### Distribution of fragment lengths
    #####-------------------------------------------------------------#####
    def generate_flen_feature(self, 
                              save_feature = True):
        flendf = self.maindf[["flen"]].copy()
        flendf["abs_flen"] = flendf["flen"].abs()
        if not flendf.empty:
            flen_count = flendf["abs_flen"].value_counts().reset_index()
            flen_count.columns = ["size", "count"]
            ##### keep only fragments that are between 50 and 350 bp
            flen_count = flen_count[(flen_count["size"] >= 50) & (flen_count["size"] <= 350)]
            flen_count["freq"] = flen_count["count"] / flen_count["count"].sum()
            flen_count = flen_count.sort_values("size")
            output_flendf = pd.DataFrame({"size": range(50, 351)})
            output_flendf = output_flendf.merge(flen_count, on="size", how="left").fillna(0)
            output_flendf = output_flendf[["size", "freq", "count"]]
            if save_feature:
                output_flendf.to_csv(os.path.join(self.outputdir, f"{self.sampleid}_GWfeature_flen.csv"), index=False)
            return output_flendf
        
    #####-------------------------------------------------------------#####
    ##### 4bp end motif
    #####-------------------------------------------------------------#####    
    def generate_em_feature(self, 
                            save_feature = True):
        if self.feature_version == "20241001":
            emdf1 = self.maindf[["reverse_EM", "QC", "flen"]].copy() 
            emdf2 = self.maindf[["forward_EM", "QC", "flen"]].copy()
            emdf1 = emdf1[(emdf1["QC"] >= 30) & (emdf1["flen"] < 0)].drop(["QC", "flen"], axis = 1)
            emdf2 = emdf2[(emdf2["QC"] >= 30) & (emdf2["flen"] > 0)].drop(["QC", "flen"], axis = 1)
        else:
            emdf1 = self.maindf[["reverse_EM", "flen"]].copy() 
            emdf2 = self.maindf[["forward_EM", "flen"]].copy()
            emdf1 = emdf1[(emdf1["flen"] < 0)].drop(["flen"], axis = 1)
            emdf2 = emdf2[(emdf2["flen"] > 0)].drop(["flen"], axis = 1)

        emdf1.columns = ["motif"]
        emdf2.columns = ["motif"]
        
        emdf = pd.concat([emdf1, emdf2], axis = 0)
        emdf.columns = ["motif"]
        emdf = emdf[emdf["motif"].isna() == False]
        emdf["motif"] = emdf["motif"].str.upper()
        output_emdf = emdf["motif"].value_counts().reset_index()
        if not output_emdf.empty:
            output_emdf.columns = ["motif", "count"]
            output_emdf = output_emdf[~output_emdf["motif"].str.contains("N")]
            output_emdf["freq"] = output_emdf["count"] / output_emdf["count"].sum()
            output_emdf = output_emdf[["motif", "freq"]]
            if save_feature:
                output_emdf.to_csv(os.path.join(self.outputdir, f"{self.sampleid}_GWfeature_EM.csv"), index=False)
            return output_emdf

    #####-------------------------------------------------------------#####    
    ##### distribution of distance read-to-nearest nucleosome
    #####-------------------------------------------------------------#####    
    def generate_nuc_feature(self, 
                            save_feature = True):
        nucdf1 = pd.DataFrame(data = self.maindf["reverse_NUC"].values,
                     columns = ["feat"])
        nucdf2 = pd.DataFrame(data = self.maindf["forward_NUC"].values,
                     columns = ["feat"])
        nucdf = pd.concat([nucdf1, nucdf2], axis = 0)
        nucdf = nucdf[(nucdf["feat"] >= -300) & (nucdf["feat"] <= 300)]
        output_nucdf = nucdf.reset_index().groupby("feat")["index"].count().reset_index()
        output_nucdf["index"] = output_nucdf["index"].apply(lambda x: x/output_nucdf["index"].sum())
        output_nucdf.columns = ["dist", "freq"]
        if save_feature:
            output_nucdf.to_csv(os.path.join(self.outputdir, f"{self.sampleid}_GWfeature_Nucleosome.2.csv"), index=False)
        return nucdf
    
    def generate_nuc_feature_1(self, 
                               save_feature = True):
        if self.path_to_old_nuc != "none":
            print("Generate features Nucleosome from old data, bedtools closest -t all, not -t first ...")
            output_nucdf = pd.read_csv(self.path_to_old_nuc, sep = "\t", header = None)
            output_nucdf = output_nucdf[(output_nucdf[3] >= -300) & (output_nucdf[3] <= 300)]
            output_nucdf = output_nucdf.groupby(3)[0].count().reset_index()
            output_nucdf.columns = ["dist", "freq"]
            sum_nuc = output_nucdf["freq"].sum()
            output_nucdf["freq"] = output_nucdf["freq"].apply(lambda x: x/sum_nuc)
            if save_feature:
                output_nucdf.to_csv(os.path.join(self.outputdir, f"{self.sampleid}_GWfeature_Nucleosome.csv"), index=False)
        else:
            error = "Please provide the path to the old nucleosome feature file"
            raise ValueError(error)
    
    #####-------------------------------------------------------------#####
    ##### generate NDR features
    #####-------------------------------------------------------------#####
    def generate_ndr_feature(self, 
                            save_feature = True):
        NDRdf1 = pd.DataFrame(data = self.maindf_filter_chr["reverse_NDR"].values,
                     columns = ["feat"])
        NDRdf2 = pd.DataFrame(data = self.maindf_filter_chr["forward_NDR"].values,
                     columns = ["feat"])
        NDRdf = pd.concat([NDRdf1, NDRdf2], axis = 0)
        NDRdf = NDRdf[(NDRdf["feat"] >= -1000) & (NDRdf["feat"] <= 1000)]
        output_NDRdf = NDRdf.reset_index().groupby("feat")["index"].count().reset_index()
        output_NDRdf["index"] = output_NDRdf["index"].apply(lambda x: x/output_NDRdf["index"].sum())
        output_NDRdf.columns = ["dist", "freq"]
        if save_feature:
            output_NDRdf.to_csv(os.path.join(self.outputdir, f"{self.sampleid}_GWfeature_NDR.csv"), index=False)
        return NDRdf

    #####-------------------------------------------------------------#####    
    ##### EM - FLEN features
    #####-------------------------------------------------------------#####    
    def generate_EM_flen_feature(self,
                                 save_feature = True):
            # Forward EM
            feature_df = self.maindf_filter_chr.copy()
            feature_df = feature_df[(feature_df["forward_EM"].isna() == False) & 
                                    (feature_df["reverse_EM"].isna() == False)]
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
            if save_feature:
                scaled_countdf.to_csv(os.path.join(self.outputdir, f"{self.sampleid}_EM_FLEN.csv"), index=False)

    #####-------------------------------------------------------------#####    
    ##### FORWARD NUCLEOSOME - FLEN
    #####-------------------------------------------------------------#####    
    def generate_forwardNUC_flen_feature(self,
                                        save_feature = True):
        feature_df = self.maindf_filter_chr.copy()
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
        print(nuc_countdf.shape)
        assert nuc_countdf.shape[0] == 211, f"[NUC forward - flen] flen failed!"
        assert nuc_countdf.shape[1] == 601, f"[NUC forward - flen] NUC failed"
        if save_feature:
            nuc_countdf.to_csv(os.path.join(self.outputdir, f"{self.sampleid}_forwardNUC_FLEN.csv"), index=False)

    #####-------------------------------------------------------------#####    
    ##### FORWARD NDR - FLEN
    #####-------------------------------------------------------------#####    
    def generate_forwardNDR_flen_feature(self,
                                        save_feature = True):
        feature_df = self.maindf_filter_chr.copy()
        NDRdf_forward = feature_df[["flen", "forward_NDR"]].copy()
        NDRdf_forward.columns = ["flen", "NDR_dist"]
        NDRdf_forward = NDRdf_forward[
            (NDRdf_forward["NDR_dist"] <= 1000) & (NDRdf_forward["NDR_dist"] >= -1000)
        ]

        NDR_countdf = NDRdf_forward.reset_index() \
                           .groupby(["NDR_dist", "flen"])["index"] \
                           .count() \
                           .reset_index() \
                           .pivot_table(index='flen', 
                                     columns='NDR_dist', 
                                     values='index', 
                                     fill_value=0)
        flen_range_df = pd.DataFrame(
            {
                'flen': range(70, 281)
            }
        )
        
        NDR_countdf = pd.merge(flen_range_df,
                               NDR_countdf,
                               on='flen',
                               how='outer')
        NDR_countdf.fillna(0, 
                           inplace=True)

        NDR_countdf = NDR_countdf[(NDR_countdf['flen'] >= 70) & (NDR_countdf['flen'] <= 280)]

        NDR_countdf = NDR_countdf.set_index("flen")
        NDR_countdf = NDR_countdf/NDR_countdf.sum().sum()
        if save_feature:
            NDR_countdf.to_csv(os.path.join(self.outputdir, f"{self.sampleid}_forwardNDR_FLEN.csv"), index=False)

    #####-------------------------------------------------------------#####    
    ##### REVERSE NUCLEOSOME - FLEN
    #####-------------------------------------------------------------#####  
    def generate_reverseNUC_flen_feature(self,
                                         save_feature = True):
        feature_df = self.maindf_filter_chr.copy()
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
        
        if save_feature:
            nuc_countdf.to_csv(os.path.join(self.outputdir, f"{self.sampleid}_reverseNUC_FLEN.csv"), index=False)

    #####-------------------------------------------------------------#####    
    ##### REVERSE NDR - FLEN
    #####-------------------------------------------------------------#####  
    def generate_reverseNDR_flen_feature(self,
                                         save_feature = True):
        feature_df = self.maindf_filter_chr.copy()
        NDRdf_reverse = feature_df[["flen", "reverse_NDR"]].copy()
        NDRdf_reverse.columns = ["flen", "NDR_dist"]
        NDRdf_reverse = NDRdf_reverse[
            (NDRdf_reverse["NDR_dist"] <= 1000) & (NDRdf_reverse["NDR_dist"] >= -1000)
        ]

        NDR_countdf = NDRdf_reverse.reset_index() \
                           .groupby(["NDR_dist", "flen"])["index"] \
                           .count() \
                           .reset_index() \
                           .pivot_table(index='flen', 
                                     columns='NDR_dist', 
                                     values='index', 
                                     fill_value=0)
        flen_range_df = pd.DataFrame(
            {
                'flen': range(70, 281)
            }
        )
        NDR_countdf = pd.merge(flen_range_df,
                               NDR_countdf,
                               on='flen',
                               how='outer')
        NDR_countdf.fillna(0, 
                           inplace=True)

        NDR_countdf = NDR_countdf[(NDR_countdf['flen'] >= 70) & (NDR_countdf['flen'] <= 280)]

        NDR_countdf = NDR_countdf.set_index("flen")
        NDR_countdf = NDR_countdf/NDR_countdf.sum().sum()
        
        assert NDR_countdf.shape[0] == 211, f"[NDR reverse - flen] flen failed!"
        assert NDR_countdf.shape[1] == 601, f"[NDR reverse - flen] NDR failed"
        
        if save_feature:
            NDR_countdf.to_csv(os.path.join(self.outputdir, f"{self.sampleid}_reverseNDR_FLEN.csv"), index=False)

    #####-------------------------------------------------------------#####
    ##### EM - EM, all flen
    #####-------------------------------------------------------------#####
    def generate_EM_pairs_all_flen(self,
                              save_feature = True):
        feature_df = self.maindf_filter_chr.copy()
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
        
        if save_feature:
            count_pair_EM.to_csv(os.path.join(self.outputdir, f"{self.sampleid}_EM_all_fragments.csv"), index=False)
    
    #####-------------------------------------------------------------#####
    ##### EM -EM, short fragments only
    #####-------------------------------------------------------------#####
    def generate_EM_pairs_short_flen(self,
                              save_feature = True):
        feature_df = self.maindf_filter_chr.copy()
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

        if save_feature:
            count_pair_EM.to_csv(os.path.join(self.outputdir, f"{self.sampleid}_EM_short_fragments.csv"), index=False)
    
    #####-------------------------------------------------------------#####
    ##### EM - EM, long fragments only
    #####-------------------------------------------------------------#####
    def generate_EM_pairs_long_flen(self,
                              save_feature = True):
        feature_df = self.maindf_filter_chr.copy()
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
        
        if save_feature:
            count_pair_EM.to_csv(os.path.join(self.outputdir, f"{self.sampleid}_EM_long_fragments.csv"), index=False)
    
    #####-------------------------------------------------------------#####
    ##### all EM - forward NUC
    #####-------------------------------------------------------------#####
    def generate_allEM_forwardNUC(self,
                              save_feature = True):
        # IMPORTANT NOTE: TAKE BOTH REVERSE AND FORWARD EM + FORWARD NUC
        feature_df = self.maindf_filter_chr.copy()
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
                'forward_NUC': range(-300, 301)
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
        if save_feature:
            # countdf.to_csv(os.path.join(self.outputdir, f"{self.sampleid}_EM_forwardNUC.csv"), index=False) # <<<<< OLD NAME!!!!!
            countdf.to_csv(os.path.join(self.outputdir, f"{self.sampleid}_allEM_forwardNUC.csv"), index=False)

    #####-------------------------------------------------------------#####
    ##### all EM - reverse NUC
    #####-------------------------------------------------------------#####
    def generate_allEM_reverseNUC(self,
                              save_feature = True):
        # IMPORTANT NOTE: TAKE BOTH REVERSE AND FORWARD EM + REVERSE NUC
        feature_df = self.maindf_filter_chr.copy()
        ##### generate EM - REVERSE_NUC dataframe
        # Forward EM - Forward nucleosome distance
        forward_em_reverse_NUC = feature_df[["forward_EM", "reverse_NUC"]].copy()
        forward_em_reverse_NUC.columns = ["EM", "reverse_NUC"]
        
        # Reverse EM - forward nucleosome distance
        reverse_em_reverse_NUC = feature_df[["reverse_EM", "reverse_NUC"]].copy()
        reverse_em_reverse_NUC.columns = ["EM", "reverse_NUC"]

        em_reverse_NUC_df = pd.concat([forward_em_reverse_NUC, reverse_em_reverse_NUC], axis = 0)
        em_reverse_NUC_df = em_reverse_NUC_df[~em_reverse_NUC_df["EM"].str.contains("N")]
        em_reverse_NUC_df = em_reverse_NUC_df[(em_reverse_NUC_df["reverse_NUC"] >= -300) & (em_reverse_NUC_df["reverse_NUC"] <= 300)]
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
                'reverse_NUC': range(-300, 301)
            }
        )

        countdf = pd.merge(reverse_NUC_range_df, countdf, on='reverse_NUC', how='outer')
        countdf.fillna(0, inplace=True)

        countdf = countdf.set_index("reverse_NUC")
        missing_motifs = [item for item in countdf.columns if item not in self.all_4bp_motifs]

        if len(missing_motifs) != 0:
            for motif in missing_motifs:
                countdf[motif] = 0
                
        countdf = countdf/countdf.sum().sum()
        countdf = countdf[self.motif_order]
        if save_feature:
            countdf.to_csv(os.path.join(self.outputdir, f"{self.sampleid}_allEM_reverseNUC.csv"), index=False)

    #####-------------------------------------------------------------#####
    ##### reverse EN .reverse NUC
    #####-------------------------------------------------------------#####
    def generate_reverseEM_reverseNUC(self,
                              save_feature = True):
        # IMPORTANT NOTE: JUST TAKE REVERSE EM + REVERSE NUC
        feature_df = self.maindf_filter_chr.copy()
        reverse_em_reverse_NUC = feature_df[["reverse_EM", "reverse_NUC"]].copy()
        reverse_em_reverse_NUC.columns = ["EM", "reverse_NUC"]
        em_reverse_NUC_df = reverse_em_reverse_NUC.copy()
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
                'reverse_NUC': range(-300, 301)
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
        if save_feature:
            # countdf.to_csv(os.path.join(self.outputdir, f"{self.sampleid}_EM_reverseNUC.csv"), index=False) #  <<<<< OLD NAME!!!!!
            countdf.to_csv(os.path.join(self.outputdir, f"{self.sampleid}_reverseEM_reverseNUC.csv"), index=False)
    
    #####-------------------------------------------------------------#####
    ##### forward EN .forward NUC
    #####-------------------------------------------------------------#####
    def generate_forwardEM_forwardNUC(self,
                              save_feature = True):
        # IMPORTANT NOTE: JUST TAKE forward EM + forward NUC
        feature_df = self.maindf_filter_chr.copy()
        forward_em_forward_NUC = feature_df[["forward_EM", "forward_NUC"]].copy()
        forward_em_forward_NUC.columns = ["EM", "forward_NUC"]
        em_forward_NUC_df = forward_em_forward_NUC.copy()
        em_forward_NUC_df = em_forward_NUC_df[~em_forward_NUC_df["EM"].str.contains("N")]
        
        em_forward_NUC_df = em_forward_NUC_df[
            (em_forward_NUC_df["forward_NUC"] >= -300) 
            & (em_forward_NUC_df["forward_NUC"] <= 300)
        ]
        
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
                'forward_NUC': range(-300, 301)
            }
        )
        countdf = pd.merge(forward_NUC_range_df, 
                           countdf, 
                           on='forward_NUC', 
                           how='outer')
        countdf.fillna(0, 
                       inplace=True)

        
        countdf = countdf.set_index("forward_NUC")
        missing_motifs = [item for item in countdf.columns if item not in self.all_4bp_motifs]

        if len(missing_motifs) != 0:
            for motif in missing_motifs:
                countdf[motif] = 0
        countdf = countdf/countdf.sum().sum()
        countdf = countdf[self.motif_order]
        if save_feature:
            countdf.to_csv(os.path.join(self.outputdir, f"{self.sampleid}_forwardEM_forwardNUC.csv"), index=False)

    #####-------------------------------------------------------------#####
    ##### EM - forward NDR
    #####-------------------------------------------------------------#####
    def generate_allEM_forwardNDR(self,
                              save_feature = True):
        feature_df = self.maindf_filter_chr.copy()
        ##### generate EM - forward_NDR dataframe
        # Forward EM - Forward NDRleosome distance
        forward_em_forward_NDR = feature_df[["forward_EM", "forward_NDR"]].copy()
        forward_em_forward_NDR.columns = ["EM", "forward_NDR"]
        
        # Reverse EM - forward NDRleosome distance
        reverse_em_forward_NDR = feature_df[["reverse_EM", "forward_NDR"]].copy()
        reverse_em_forward_NDR.columns = ["EM", "forward_NDR"]

        em_forward_NDR_df = pd.concat([forward_em_forward_NDR, reverse_em_forward_NDR], axis = 0)
        em_forward_NDR_df = em_forward_NDR_df[~em_forward_NDR_df["EM"].str.contains("N")]
        em_forward_NDR_df = em_forward_NDR_df[
            (em_forward_NDR_df["forward_NDR"] >= -1000) & 
            (em_forward_NDR_df["forward_NDR"] <= 1000)]
        countdf = em_forward_NDR_df.reset_index() \
                                   .groupby(["EM", "forward_NDR"])["index"] \
                                   .count() \
                                   .reset_index() \
                                   .pivot_table(index='forward_NDR', 
                                                columns='EM', 
                                                values='index', 
                                                fill_value=0)

        ##### fill values so that the output matrix always 50:350 x 256
        forward_NDR_range_df = pd.DataFrame(
            {
                'forward_NDR': range(-1000, 1001)
            }
        )

        countdf = pd.merge(forward_NDR_range_df, countdf, on='forward_NDR', how='outer')
        countdf.fillna(0, inplace=True)

        countdf = countdf.set_index("forward_NDR")
        missing_motifs = [item for item in countdf.columns if item not in self.all_4bp_motifs]

        if len(missing_motifs) != 0:
            for motif in missing_motifs:
                countdf[motif] = 0
                
        countdf = countdf/countdf.sum().sum()
        countdf = countdf[self.motif_order]
        if save_feature:
            countdf.to_csv(os.path.join(self.outputdir, f"{self.sampleid}_allEM_forwardNDR.csv"), index=False)

    #####-------------------------------------------------------------#####
    ##### EM - reverse NDR
    #####-------------------------------------------------------------#####
    def generate_allEM_reverseNDR(self,
                              save_feature = True):
        feature_df = self.maindf_filter_chr.copy()
        ##### generate EM - reverse_NDR dataframe
        # Forward EM - Forward NDRleosome distance
        forward_em_reverse_NDR = feature_df[["forward_EM", "reverse_NDR"]].copy()
        forward_em_reverse_NDR.columns = ["EM", "reverse_NDR"]
        
        # Reverse EM - forward NDRleosome distance
        reverse_em_reverse_NDR = feature_df[["reverse_EM", "reverse_NDR"]].copy()
        reverse_em_reverse_NDR.columns = ["EM", "reverse_NDR"]

        em_reverse_NDR_df = pd.concat([forward_em_reverse_NDR, reverse_em_reverse_NDR], axis = 0)
        em_reverse_NDR_df = em_reverse_NDR_df[~em_reverse_NDR_df["EM"].str.contains("N")]
        em_reverse_NDR_df = em_reverse_NDR_df[
            (em_reverse_NDR_df["reverse_NDR"] >= -1000) & 
            (em_reverse_NDR_df["reverse_NDR"] <= 1000)]
        countdf = em_reverse_NDR_df.reset_index() \
                                   .groupby(["EM", "reverse_NDR"])["index"] \
                                   .count() \
                                   .reset_index() \
                                   .pivot_table(index='reverse_NDR', 
                                                columns='EM', 
                                                values='index', 
                                                fill_value=0)

        ##### fill values so that the output matrix always 50:350 x 256
        reverse_NDR_range_df = pd.DataFrame(
            {
                'reverse_NDR': range(-1000, 1001)
            }
        )

        countdf = pd.merge(reverse_NDR_range_df, countdf, on='reverse_NDR', how='outer')
        countdf.fillna(0, inplace=True)

        countdf = countdf.set_index("reverse_NDR")
        missing_motifs = [item for item in countdf.columns if item not in self.all_4bp_motifs]

        if len(missing_motifs) != 0:
            for motif in missing_motifs:
                countdf[motif] = 0
                
        countdf = countdf/countdf.sum().sum()
        countdf = countdf[self.motif_order]
        if save_feature:
            countdf.to_csv(os.path.join(self.outputdir, f"{self.sampleid}_allEM_reverseNDR.csv"), index=False)

    
    #####-------------------------------------------------------------#####
    ##### reverse EN .reverse NDR
    #####-------------------------------------------------------------#####
    def generate_reverseEM_reverseNDR(self,
                              save_feature = True):
        # IMPORTANT NOTE: JUST TAKE REVERSE EM + REVERSE NDR
        feature_df = self.maindf_filter_chr.copy()
        reverse_em_reverse_NDR = feature_df[["reverse_EM", "reverse_NDR"]].copy()
        reverse_em_reverse_NDR.columns = ["EM", "reverse_NDR"]
        em_reverse_NDR_df = reverse_em_reverse_NDR.copy()
        em_reverse_NDR_df = em_reverse_NDR_df[~em_reverse_NDR_df["EM"].str.contains("N")]
        
        em_reverse_NDR_df = em_reverse_NDR_df[
            (em_reverse_NDR_df["reverse_NDR"] >= -1000) 
            & (em_reverse_NDR_df["reverse_NDR"] <= 1000)
        ]
        
        countdf = em_reverse_NDR_df.reset_index() \
                                    .groupby(["EM", "reverse_NDR"])["index"] \
                                    .count() \
                                    .reset_index() \
                                    .pivot_table(index='reverse_NDR', 
                                                 columns='EM', 
                                                 values='index', 
                                                 fill_value=0)

        ##### fill values so that the output matrix always 50:350 x 256
        reverse_NDR_range_df = pd.DataFrame(
            {
                'reverse_NDR': range(-1000, 1001)
            }
        )
        countdf = pd.merge(reverse_NDR_range_df, 
                           countdf, 
                           on='reverse_NDR', 
                           how='outer')
        countdf.fillna(0, 
                       inplace=True)

        
        countdf = countdf.set_index("reverse_NDR")
        missing_motifs = [item for item in countdf.columns if item not in self.all_4bp_motifs]

        if len(missing_motifs) != 0:
            for motif in missing_motifs:
                countdf[motif] = 0
        countdf = countdf/countdf.sum().sum()
        countdf = countdf[self.motif_order]
        if save_feature:
            countdf.to_csv(os.path.join(self.outputdir, f"{self.sampleid}_reverseEM_reverseNDR.csv"), index=False)
    
    #####-------------------------------------------------------------#####
    ##### forward EN .forward NDR
    #####-------------------------------------------------------------#####
    def generate_forwardEM_forwardNDR(self,
                              save_feature = True):
        # IMPORTANT NOTE: JUST TAKE forward EM + forward NDR
        feature_df = self.maindf_filter_chr.copy()
        forward_em_forward_NDR = feature_df[["forward_EM", "forward_NDR"]].copy()
        forward_em_forward_NDR.columns = ["EM", "forward_NDR"]
        em_forward_NDR_df = forward_em_forward_NDR.copy()
        em_forward_NDR_df = em_forward_NDR_df[~em_forward_NDR_df["EM"].str.contains("N")]
        
        em_forward_NDR_df = em_forward_NDR_df[
            (em_forward_NDR_df["forward_NDR"] >= -1000) 
            & (em_forward_NDR_df["forward_NDR"] <= 1000)
        ]
        
        countdf = em_forward_NDR_df.reset_index() \
                                    .groupby(["EM", "forward_NDR"])["index"] \
                                    .count() \
                                    .reset_index() \
                                    .pivot_table(index='forward_NDR', 
                                                 columns='EM', 
                                                 values='index', 
                                                 fill_value=0)

        ##### fill values so that the output matrix always 50:350 x 256
        forward_NDR_range_df = pd.DataFrame(
            {
                'forward_NDR': range(-1000, 1001)
            }
        )
        countdf = pd.merge(forward_NDR_range_df, 
                           countdf, 
                           on='forward_NDR', 
                           how='outer')
        countdf.fillna(0, 
                       inplace=True)

        
        countdf = countdf.set_index("forward_NDR")
        missing_motifs = [item for item in countdf.columns if item not in self.all_4bp_motifs]

        if len(missing_motifs) != 0:
            for motif in missing_motifs:
                countdf[motif] = 0
        countdf = countdf/countdf.sum().sum()
        countdf = countdf[self.motif_order]
        if save_feature:
            countdf.to_csv(os.path.join(self.outputdir, f"{self.sampleid}_forwardEM_forwardNDR.csv"), index=False)