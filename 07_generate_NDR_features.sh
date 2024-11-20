# while getopts "i:o:t:" opt; do
#   case ${opt} in
#     i )
#       finalOutputFile=$OPTARG
#       ;;
#     o )
#       outputdir=$OPTARG
#       ;;
#     t )
#       tmpdir=$OPTARG
#       ;;
#     \? )
#       echo "Usage: cmd [-i] inputbam [-o] outputdir [-n] samtools_num_threads [-t] input_type"
#       exit 1
#       ;;
#   esac
# done

# finalOutputFile="/media/hieunguyen/HNSD01/src/ecd_wgs_features/batch_test/WGShg19/WGShg19.final_output.tsv";
finalOutputFile="/media/hieunguyen/HNSD01/src/ecd_wgs_features/old_final_output.tsv";
outputdir="debug_NDR";
tmpdir="tmp_debug_NDR";
mkdir -p ${outputdir};
mkdir -p ${tmpdir};

file_name=$(basename ${finalOutputFile});
echo -e "Working on sample " $file_name;

echo -e "generating file *_reordered.bed"
cat ${finalOutputFile} | awk -v OFS='\t' '{print $2,$7,$8,$1,$3,$4,$5,$6,$9,$10,$11,$12,$13}' > ${tmpdir}/${file_name}_reordered.bed

# awk '{print NF}' NDR_cancer_specific_location_5common.HNmodified.csv | sort | uniq -c
# awk '{print NF}' ${tmpdir}/${file_name}_reordered.bed | sort | uniq -c

echo -e "generating file *_filtered.bed with betdtools intersect"
bedtools intersect -b NDR_cancer_specific_location_5common.HNmodified.csv -a ${tmpdir}/${file_name}_reordered.bed -wa > ${tmpdir}/${file_name}_filtered.bed

echo -e "generating file NDR_<filename> ..."
cat ${tmpdir}/${file_name}_filtered.bed | awk -v OFS='\t' '{print $4,$1,$5,$6,$7,$8,$2,$3,$9,$10,$11,$12,$13}' > ${tmpdir}/NDR_${file_name}
# rm ${tmpdir}/${file_name}_filtered.bed
# rm ${tmpdir}/${file_name}_reordered.bed

echo -e "Sort the NDR bed file ..."
if [ ! -f NDR_cancer_specific_location_5common_CENTER.sorted.bed ]; then
    sort -k 1V,1 -k 2n,2 NDR_cancer_specific_location_5common_CENTER.csv -o NDR_cancer_specific_location_5common_CENTER.sorted.bed
fi

  cat ${tmpdir}/NDR_${file_name} | awk -v OFS='\t' '{print $2,$3,$7,$9,$5,$10,$11,$12,$13}' > ${outputdir}/${file_name}_forward.bed
  cat ${tmpdir}/NDR_${file_name} | awk -v OFS='\t' '{print $2,$6,$8,$9,$5,$10,$11,$12,$13}' > ${outputdir}/${file_name}_reverse.bed
  sort -k 1V,1 -k 2n,2 ${outputdir}/${file_name}_forward.bed -o ${outputdir}/${file_name}_forward.sorted.bed
  sort -k 1V,1 -k 2n,2 ${outputdir}/${file_name}_reverse.bed -o ${outputdir}/${file_name}_reverse.sorted.bed

  bedtools closest -a ${outputdir}/${file_name}_forward.sorted.bed -b NDR_cancer_specific_location_5common_CENTER.sorted.bed -t first | \
  awk -v OFS='\t' '{$13=$2 - $11; print $1,$2,$4,$5,$6,$7,$8,$9,$13}' > ${outputdir}/${file_name}_forward.tsv

  bedtools closest -a ${outputdir}/${file_name}_reverse.sorted.bed -b NDR_cancer_specific_location_5common_CENTER.sorted.bed -t first | \
  awk -v OFS='\t' '{$13=$2 - $11; print $1,$2,$4,$5,$6,$7,$8,$9,$13}' > ${outputdir}/${file_name}_reverse.tsv

#   rm ${outputdir}/${file_name}_forward.bed
#   rm ${outputdir}/${file_name}_forward.sorted.bed
#   rm ${outputdir}/${file_name}_reverse.bed
#   rm ${outputdir}/${file_name}_reverse.sorted.bed

