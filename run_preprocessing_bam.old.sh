export PATH=/Users/hieunguyen/bedtools2/bin/:$PATH;
export PATH=/Users/hieunguyen/samtools/bin:$PATH;

path_to_bam_file="/Users/hieunguyen/data/bam_files/WGShg19.bam";
path_to_output="./output/test_compare_pipeline_old_new";

mkdir -p $path_to_output;

path_to_fa="/Users/hieunguyen/data/resources/hg19.fa";
sampleid="test_old_pipeline"
path_to_SRC="/Users/hieunguyen/src/ecd_wgs_features/old_version_pipelines/ecd_wgs_features_hg19/src";
path_to_REF="/Users/hieunguyen/data/resources/rpr_map_EXP0779.sorted.bed";

bash preprocessing_bam.old.sh $path_to_bam_file $path_to_output $path_to_fa $sampleid $path_to_SRC $path_to_REF;