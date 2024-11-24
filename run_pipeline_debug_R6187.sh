export PATH=/home/hieunguyen/bedtools2/bin/:$PATH
export PATH=/home/hieunguyen/samtools/bin/:$PATH

# inputbam="/media/hieunguyen/GSHD_HN01/raw_data/bam_files/WGShg19.bam";
inputbam=$1;

outputdir="./R6187/";
num_threads=20;
path_to_fa="/media/hieunguyen/GSHD_HN01/storage/resources/hg19.fa";
path_to_ref="/media/hieunguyen/GSHD_HN01/storage/resources/rpr_map_EXP0779.sorted.bed";
motif_order_path="./motif_order.csv";
final_Feature_dir="./final_Feature_dir";
ndr_ref="/media/hieunguyen/HNSD01/src/ecd_wgs_features/NDR_cancer_specific_location_5common_CENTER.sorted.bed"
ndrb_ref="/media/hieunguyen/HNSD01/src/ecd_wgs_features/NNDR_cancer_specific_location_binary.sorted.bed"

mkdir -p ${outputdir};
mkdir -p ${final_Feature_dir};

sampleid=$(echo ${inputbam} | xargs -n 1 basename)
sampleid=${sampleid%.bam*}

if [ ! -f "${outputdir}/${sampleid}/${sampleid}.final_output.tsv" ]; then
    # bash 01_BAM_to_FRAGS.sh \
    bash debug_R6187.sh \
    -i  ${inputbam} \
    -o ${outputdir} \
    -n ${num_threads} \
    -f ${path_to_fa};

    bash 02_calculate_EM_FLEN_NUC_from_FRAG.sh \
        -i ${outputdir}/${sampleid}/${sampleid}.frag.tsv  \
        -o ${outputdir} \
        -f ${path_to_fa} \
        -r ${path_to_ref} \
        -n ${ndr_ref} \
        -b ${ndrb_ref} \
        -c false;
        
    echo -e ${outputdir}/${sampleid}/${sampleid}.final_output.tsv "exists";
    echo -e "Running script 03 to generate features ..."
    python 03_generate_WGS_features.py \
        --input ${outputdir}/${sampleid}/${sampleid}.final_output.tsv \
        --output ${outputdir}/${sampleid} \
        --motif_order_path ${motif_order_path} \
        --feature_version "20241001" \
        --old_nuc ${outputdir}/${sampleid}/${sampleid}.full_Nucleosome.dist.final.bed \
        --generate_feature "all" 

    echo -e "Finished generating features, saving csv files"

    echo -e "Running script 04 to generate feature matrix ..."
    python 04_generate_batch_feature_matrix.py \
        --input ${outputdir} \
        --output ${final_Feature_dir}
    echo -e "Finished generating feature matrix"

else 
    # echo -e ${outputdir}/${sampleid}/${sampleid}.final_output.tsv "exists";
    echo -e "Running script 03 to generate features ..."
    python 03_generate_WGS_features.py \
        --input ${outputdir}/${sampleid}/${sampleid}.final_output.tsv \
        --output ${outputdir}/${sampleid} \
        --motif_order_path ${motif_order_path} \
        --feature_version "20241001" \
        --old_nuc ${outputdir}/${sampleid}/${sampleid}.full_Nucleosome.dist.final.bed \
        --generate_feature "all" 

    echo -e "Finished generating features, saving csv files"

    echo -e "Running script 04 to generate feature matrix ..."
    python 04_generate_batch_feature_matrix.py \
        --input ${outputdir} \
        --output ${final_Feature_dir}
    echo -e "Finished generating feature matrix"
    
fi



