
path_to_final_output="/Users/hieunguyen/src/ecd_wgs_features/batch_test_FinaleDB/EE86537.hg19_prep/EE86537.hg19_prep.final_output.tsv";
path_to_old_nuc="/Users/hieunguyen/src/ecd_wgs_features/batch_test_FinaleDB/EE86537.hg19_prep/EE86537.hg19_prep.full_Nucleosome.dist.final.bed";
outputdir="/Users/hieunguyen/src/ecd_wgs_features/batch_test_FinaleDB";
motif_order_path="./motif_order.csv";
final_Feature_dir=${outputdir}
echo -e "Running script 03 to generate features ..."
python 03_generate_WGS_features.py \
    --input ${path_to_final_output} \
    --output ${outputdir}/EE86537 \
    --motif_order_path ${motif_order_path} \
    --feature_version "20241001" \
    --old_nuc ${path_to_old_nuc} \
    --generate_feature "all" 

echo -e "Finished generating features, saving csv files"
python 04_generate_batch_feature_matrix.py \
    --input ${outputdir} \
    --output ${final_Feature_dir}