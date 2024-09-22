import argparse
import pandas as pd 
import os

##### Feature class definition

class WGS_GW_features:
    def __init__(self,
                 input_tsv,
                 motif_order_path,
                 outputdir):
        self.input_tsv = input_tsv
        self.sampleid = input_tsv.split("/")[-1].split(".")[0]
        print("reading in the input frag.tsv data")
        self.maindf = pd.read_csv(input_tsv, sep = "\t", header = None)
        self.maindf.columns = ["chr", "start", "end", "flen", "readID", "forward_NUC", "reverse_NUC", "forward_EM", "reverse_EM"]
        self.motif_order_path = motif_order_path
        self.motif_order = pd.read_csv(motif_order_path)["motif_order"].values
        self.all_4bp_motifs = [
            "{}{}{}{}".format(i,j,k,l) 
            for i in ["A", "T", "G", "C"] 
            for j in ["A", "T", "G", "C"] 
            for k in ["A", "T", "G", "C"] 
            for l in ["A", "T", "G", "C"]
        ]
        self.maindf_filter_chr = self.maindf[self.maindf["chr"].isin([f"chr{i}" for i in range(1, 22)])]
        self.outputdir = outputdir
        
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
                output_flendf.to_csv(os.path.join(self.outputdir, f"{self.sampleid}.flen.csv"), index=False)
            return output_flendf
        
    #####-------------------------------------------------------------#####
    ##### 4bp end motif
    #####-------------------------------------------------------------#####    
    def generate_em_feature(self, 
                            save_feature = True):
        emdf = pd.DataFrame(data = [item for item in self.maindf["reverse_EM"].values + self.maindf["forward_EM"].values if item != "NA"],
                            columns = ["EM"])
        emdf.columns = ["motif"]
        emdf["motif"] = emdf["motif"].str.upper()
        output_emdf = emdf["motif"].value_counts().reset_index()
        if not output_emdf.empty:
            output_emdf.columns = ["motif", "count"]
            output_emdf = output_emdf[~output_emdf["motif"].str.contains("N")]
            output_emdf["freq"] = output_emdf["count"] / output_emdf["count"].sum()
            output_emdf = output_emdf[["motif", "freq"]]
            if save_feature:
                output_emdf.to_csv(os.path.join(self.outputdir, f"{self.sampleid}.EM.csv"), index=False)
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
        output_nucdf["index"].plot()
        if save_feature:
            output_nucdf.to_csv(os.path.join(self.outputdir, f"{self.sampleid}.NUC.csv"), index=False)
        return nucdf

def main():
    parser = argparse.ArgumentParser(description='Generate an image matrix.')
    parser.add_argument('--input', type=str, required=True, help='Path to the pre-processed frag.tsv files from 01 and 02')
    parser.add_argument('--output', type=str, required=False, help='Path to save output feature csv files')
    parser.add_argument('--motif_order_path', type=str, required=True, help='Path to the motif order file')
    
    args = parser.parse_args()
    input_tsv = args.input
    motif_order_path = args.motif_order_path 
    outputdir = args.output
    
    output_obj = WGS_GW_features(input_tsv = input_tsv,
                             motif_order_path = motif_order_path,
                             outputdir = outputdir)
    
    ##### generate and save features to output dir
    output_obj.generate_flen_feature()
    output_obj.generate_em_feature()
    output_obj.generate_nuc_feature()
    
if __name__ == '__main__':
    main()