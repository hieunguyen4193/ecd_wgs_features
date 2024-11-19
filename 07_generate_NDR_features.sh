while getopts "i:o:t:" opt; do
  case ${opt} in
    i )
      finalOutputFile=$OPTARG
      ;;
    o )
      outputdir=$OPTARG
      ;;
    t )
      tmpdir=$OPTARG
      ;;
    \? )
      echo "Usage: cmd [-i] inputbam [-o] outputdir [-n] samtools_num_threads [-t] input_type"
      exit 1
      ;;
  esac
done
mkdir -p ${outputdir};
mkdir -p ${tmpdir};

file_name=$(basename ${finalOutputFile});
echo "$file_name"

cat ${finalOutputFile} | awk -v OFS='\t' '{print $2,$7,$8,$1,$3,$4,$5,$6,$9,$10,$11,$12,$13}' > ${tmpdir}/${file_name}_reordered.bed
bedtools intersect -b NDR_cancer_specific_location_5common.csv -a ${tmpdir}/${file_name}_reordered.bed -wa > ${tmpdir}/${file_name}_filtered.bed
cat ${tmpdir}/${file_name}_filtered.bed | awk -v OFS='\t' '{print $4,$1,$5,$6,$7,$8,$2,$3,$9,$10,$11,$12,$13}' > ${tmpdir}/NDR_${file_name}
rm ${tmpdir}/${file_name}_filtered.bed
rm ${tmpdir}/${file_name}_reordered.bed

if [ ! -f NDR_cancer_specific_location_5common_CENTER.sorted.bed ]; then
    sort -k 1V,1 -k 2n,2 NDR_cancer_specific_location_5common_CENTER.csv -o NDR_cancer_specific_location_5common_CENTER.sorted.bed
fi

  cat ${tmpdir}/NDR_${file_name} | awk -v OFS='\t' '{print $2, $3,$7,$9,$5,$10,$11,$12,$13}' > ${outputdir}/${file_name}_forward.bed
  cat ${tmpdir}/NDR_${file_name} | awk -v OFS='\t' '{print $2, $6,$8,$9,$5,$10,$11,$12,$13}' > ${outputdir}/${file_name}_reverse.bed
  sort -k 1V,1 -k 2n,2 ${outputdir}/${file_name}_forward.bed -o ${outputdir}/${file_name}_forward.sorted.bed
  sort -k 1V,1 -k 2n,2 ${outputdir}/${file_name}_reverse.bed -o ${outputdir}/${file_name}_reverse.sorted.bed

  bedtools closest -a ${outputdir}/${file_name}_forward.sorted.bed -b NDR_cancer_specific_location_5common_CENTER.sorted.bed -t first | \
  awk -v OFS='\t' '{$13=$2 - $11; print $1,$2,$4,$5,$6,$7,$8,$9,$13}' > ${outputdir}/${file_name}_forward.tsv

  bedtools closest -a ${outputdir}/${file_name}_reverse.sorted.bed -b NDR_cancer_specific_location_5common_CENTER.sorted.bed -t first | \
  awk -v OFS='\t' '{$13=$2 - $11; print $1,$2,$4,$5,$6,$7,$8,$9,$13}' > ${outputdir}/${file_name}_reverse.tsv

  rm ${outputdir}/${file_name}_forward.bed
  rm ${outputdir}/${file_name}_forward.sorted.bed
  rm ${outputdir}/${file_name}_reverse.bed
  rm ${outputdir}/${file_name}_reverse.sorted.bed