inputbam="/media/hieunguyen/GSHD_HN01/raw_data/bam_files/WGShg19.bam";
outputdir="./batch_test/";
num_threads=20;
path_to_fa="/media/hieunguyen/GSHD_HN01/storage/resources/hg19.fa";
path_to_ref="/media/hieunguyen/GSHD_HN01/storage/resources/rpr_map_EXP0779.sorted.bed";
motif_order_path="./motif_order.csv";
final_Feature_dir="./final_Feature_dir";

mkdir -p ${outputdir};
mkdir -p ${final_Feature_dir};

sampleid=$(echo ${inputbam} | xargs -n 1 basename)
sampleid=${sampleid%.bam*}

bash 06_preprocess_bam_for_CNA_ratioFLEN_features.sh \
    -i ${inputbam} \
    -o ${outputdir} \
    -n ${num_threads} \
    -t "raw";