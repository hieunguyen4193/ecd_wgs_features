#!/bin/bash 
#SBATCH --job-name=MRD
#SBATCH --partition=research_hi
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH --output=/mnt/DATASM14/hieunho/hieu_pipeline/MRD_GW_UMI_Phuc/log.log

source activate /home/hieutran/.conda/envs/hieunho

############################ NEED TO INPUT ##################################
input_fastq=/mnt/NAS_PROJECT/vol_ECDteam/hieunho/link_fastq/Research-mrdgw-all-060624 # path_to_input_dir
output_dir=/datassd/hieunho/MRD_GW_UMI_Phuc/Research-mrdgw-all-060624 # path_to_output_dir
source_dir=/mnt/DATASM14/hieunho/hieu_pipeline/MRD_GW_UMI_Phuc # path_to_project_dir
##############################################################################

mkdir -p $output_dir
work=$output_dir/work

# if run slurm
nextflow run $source_dir \
    -profile slurm,docker \
    --fqdir $input_fastq --outdir $output_dir -w $work