input_FinaleDB_frag="/Users/hieunguyen/data/FinaleDB/EE87519.hg19.frag.tsv";
outputdir="./batch_test_FinaleDB/";
num_threads=10;
path_to_fa="/Users/hieunguyen/data/resources/hg19.fa";
path_to_ref="/Users/hieunguyen/data/resources/rpr_map_EXP0779.sorted.bed";
motif_order_path="./motif_order.csv";
final_Feature_dir="./final_Feature_dir";

mkdir -p ${outputdir};
mkdir -p ${final_Feature_dir};

sampleid=$(echo ${inputbam} | xargs -n 1 basename)
sampleid=${sampleid%.bam*}

if [ ! -f "${outputdir}/${sampleid}/${sampleid}.final_output.tsv" ]; then
    echo -e "preprocessing the input file from FinaleDB ..."
    sampleid=$(echo ${input_FinaleDB_frag} | xargs -n 1 basename)
    sampleid=${sampleid%.frag.tsv*}
    bash preprocess_FinaleDB_input.sh -i ${input_FinaleDB_frag} -o ${outputdir}/${sampleid};
    echo -e "old input frag: " ${input_FinaleDB_frag}
    input_FinaleDB_frag=${outputdir}/${sampleid}/${sampleid}.frag.tsv;
    echo -e "new input frag: " ${input_FinaleDB_frag}
    echo -e "finished pre-processing input finaleDB" 

    echo -e "Running script 02 to calculate EM and FLEN from FRAG file ..."
    bash 02_calculate_EM_FLEN_NUC_from_FRAG.sh \
        -i ${input_FinaleDB_frag}  \
        -o ${outputdir} \
        -f ${path_to_fa} \
        -r ${path_to_ref} \
        -c false;
else 
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
    python 04_generate_batch_feature_matrix.py \
        --input ${outputdir} \
        --output ${final_Feature_dir}
fi



