import numpy as np 
import pandas as pd 

import matplotlib.pyplot as plt
from helper_functions_from_TSMA import *
import pysam

path_to_bam_file = "/media/hieunguyen/GSHD_HN01/raw_data/bam_files/WGShg19.bam"
path_to_all_fa = "/media/hieunguyen/GSHD_HN01/storage/resources/hg19"
path_to_bed_file = "./methyl_regions/TSMA.bed"

# bedfile = pd.read_csv(path_to_bed_file, sep="\t", header=None)
# bedfile.columns = ["chrom", "start", "end", "region_name"]
# bedfile["region"] = bedfile[["chrom", "start", "end"]].apply(lambda x: "{}:{}-{}".format(x[0], x[1], x[2]), axis=1)

cna_bin = pd.read_csv("./methyl_regions/CNA_bins.bed", header = None)[0].values

outputdf = pd.DataFrame({"CGN_motif": ["CGN", "NCG", "NA"]})
for bin in tqdm(cna_bin):
    readdf = fetch_reads(path_to_bam_file, "chr{}".format(bin))

    # cnadf = pd.read_csv("/home/hieunguyen/Downloads/1-0ACKD53A11_S7501-S7701.CNA.csv")
    # cnadf[(cnadf["Unnamed: 0"].str.contains("X") == False) & (cnadf["Unnamed: 0"].str.contains("Y") == False)][["Unnamed: 0"]].to_csv("./methyl_regions/CNA_bins.bed", 
    #                                                                                                                                   sep = "\t", 
    #                                                                                                                                   header = False, 
    #                                                                                                                                   index = False)

    sampleid = str(path_to_bam_file).split("/")[-1].split(".")[0]

    region_chrom = bin.split(":")[0]
    region_start = int(bin.split(":")[1].split("-")[0])
    region_end = int(bin.split(":")[1].split("-")[1])

    refseq_at_cluster = get_refseq(path_to_all_fa = path_to_all_fa, 
                                    chrom = region_chrom, 
                                    start = region_start, 
                                    end = region_end + 1)
    all_cpg_in_cluster = [m.start(0) for m in re.finditer("CG", refseq_at_cluster)]
    cpg_coords = [item + region_start for item in all_cpg_in_cluster]
    cpg_coords_minus1 = [item -1 for item in cpg_coords]

    cpg_coords = [str(item) for item in cpg_coords]
    cpg_coords_minus1 = [str(item) for item in cpg_coords_minus1]

    def assign_status_read_start(x, cpg_coords = cpg_coords, cpg_coords_minus1 = cpg_coords_minus1):
        if x in cpg_coords:
            return "CGN"
        elif x in cpg_coords_minus1:
            return "NCG"
        else:
            return "NA"
    readdf["check_CG_motif"] = readdf["start"].apply(lambda x: assign_status_read_start(x))

    tmp_outputdf = readdf.groupby("check_CG_motif")["start"].count().reset_index()
    tmp_outputdf.columns = ["CGN_motif", bin]
    outputdf = outputdf.merge(tmp_outputdf,right_on = "CGN_motif", left_on = "CGN_motif")