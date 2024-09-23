# Generate features from low-depth WGS data

## Pre-processing `BAM` file

If we use `BAM` file as input, we perform some pre-processing steps. Run `01_BAM_to_FRAGS.sh`. Example command: 

```sh
bash 01_BAM_to_FRAGS.sh -i /path/to/sorted_and_indexed_bam.bam  -o /path/to/outputdir -n ${num_threads_use_in_samtools}
```

This function takes a sorted and indexed `BAM` file as input and generate a `.frag.tsv` files. The first 3 columns of this `.frag.tsv` file are: `chromosome`, `fragment start` 
and `fragment end`. 

If the input file is already in a `.frag.tsv` format, skip `01_BAM_to_FRAGS.sh` and just run `02_calculate_EM_FLEN_NUC_from_FRAG.sh`.

Example command for `02_calculate_EM_FLEN_NUC_from_FRAG.sh`
```sh
bash 02_calculate_EM_FLEN_NUC_from_FRAG.sh \
    -i /path/to/file.frag.tsv  \
    -o /path/to/output \
    -f /path/to/hg19.fa \
    -r /path/to/rpr_map_Budhraja_STM2023.bed
```

The file `hg19.fa` can be downloaded from UCSC, `rpr_map_Budhraja_STM2023.bed` can be download from  https://zenodo.org/records/7402091.

## Generate ready-to-use feature `.csv` files

To generate ready-to-use feature matrix, we take the file `SampleID.final_output.tsv` output from `02_calculate_EM_FLEN_NUC_from_FRAG.sh` as input. 

The script `03_generate_WGS_features.py` generates `EM`, `FLEN`, `NUCLEOSOME` (distance of reads to nearest nucleosome) and all 11 paired-feature in *matrix/image* format. 

Example command:
```sh
python 03_generate_WGS_features.py --input /path/to/SampleID.final_output.tsv --output /path/to/output --motif_order_path /path/to/motif_order.csv
```

### Bin-wise (1M) features: CNA and Flen ratio
To generate CNA and fragment length ratio in 1M bin across the whole genome, split the input `BAM` file into short, long fragment `BAM` files and run the script `04_generate_BinWise_features.R`.

Example command:
```sh
split_BAM_to_short_long_BAMs.sh -i /path/to/file.bam -o /path/to/output;

Rscript 04_generate_BinWise_features.R --input_short_bam /path/to/file.short.bam --input_long_bam /path/to/file.long.bam --input_full_bam /path/to/file.bam --output /path/to/output 
```