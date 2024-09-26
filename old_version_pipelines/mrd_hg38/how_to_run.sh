############################ NEED TO INPUT ##################################
input_fastq=".." # path_to_input_dir
output_dir=".." # path_to_output_dir
source_dir="../MRD_GW_UMI_Phuc_HG38" # path_to_project_dir
##############################################################################

mkdir -p $output_dir
cd $output_dir
work=$output_dir/work

# if run slurm
nextflow run $source_dir \
    -profile slurm,docker \
    --fqdir $input_fastq --outdir $output_dir -w $work

# # if run local
# nextflow run $source_dir \
#     -profile docker \
#     --fqdir $input_fastq --outdir $output_dir -w $work