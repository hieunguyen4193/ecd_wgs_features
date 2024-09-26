
# ECD_WGS_features_HG38

Description: nonBS samples with HG38


## How to use

Step 1: Clone the repository.

    git clone https://gitlab.genesolutions.vn/ecd-data/ecd_wgs_features_hg38.git

Step 2:  Change to the project directory.

    cd ecd_wgs_features_hg38

Step 3: Fill in the `run/run_template.sh` file with the following variables:
    
- `log`: Path to the log file.

- `input_dir`: Path to the input directory (fastq).

- `output_dir`: Path to the output directory.

- `working_dir`: Path to the project directory (ecd_wgs_features_hg38).

- `submit_type`: How to submit job. Options: `'slurm'`, `'local'`.

- `BWA_REPO`: Path to the HG38 directory (HPC-R:/datassd/hieunho/hg38_selected)

- `resource`: Path to the resource directory (HPC-R:/datassd/hieunho/ECD_WGS_resource).


Step 4: Submit the job using SLURM.

    sbatch run/run_template.sh

## Docker image

`biocontainers/fastqc:v0.11.9_cv8`: FastQC Version v0.11.9

`staphb/trimmomatic:latest`: Trimmomatic Version 0.39

`quay.io/biocontainers/bwa:0.7.17--h5bf99c6_8`: Bwa Version 0.7.17

`quay.io/biocontainers/samtools:1.13--h8c37831_0`: Samtools Version 1.13
`gene110/gw_pipeline:v2`: Python Version 3.10.6

`tronghieunguyen/ecd_features:latest`: R version 4.3.3

`gene110/samtools_python:v5`: Python Version 3.8.10 and Samtools Version 1.16

`dqpham/ichorcna:latest`: IchorCNA
