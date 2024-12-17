export PATH=/Users/hieunguyen/samtools/bin:$PATH 
# bash 01_BAM_to_FRAGS.sh -i /Users/hieunguyen/data/tmp/debug_ecd_wgs_features/data/1-ZLAAO90NB_S7509-S7709.sorted.bam -o ./output_debug/ -n 10 -f /Users/hieunguyen/data/resources/hg19.fa 
#####----------------------------------------------------------------------#####
##### input args
#####----------------------------------------------------------------------#####
while getopts "i:o:n:d:r:" opt; do
  case ${opt} in
    i )
      finalOutputFile=$OPTARG
      ;;
    o )
      outputdir=$OPTARG
      ;;
    n )
      ndr_ref1=$OPTARG
      ;;
    d )
      ndr_ref2=$OPTARG
      ;;
    r )
      ndr_ref3=$OPTARG
      ;;
      

    \? )
      echo "Usage: cmd [-i] finalOutputFile [-o] outputdir [-n] NDR ref 1 [-d] NDR ref 2 [-r] NDR ref 3"
      exit 1
      ;;
  esac
done


sampleid=$(echo ${finalOutputFile} | xargs -n 1 basename)
sampleid=${sampleid%.final_output.tsv*}
outputdir=${outputdir}/${sampleid}

mkdir -p ${outputdir}

cat ${finalOutputFile} | cut -f1,2,4,5 | \
    awk -v OFS='\t' '{$5=$2 + 1; print $1 "\t" $2 "\t" $5 "\t" $4 "\t" $3}' \
    > ${outputdir}/${sampleid}.forward_NDR.bed
cat ${finalOutputFile} | cut -f1,3,4,5 | \
awk -v OFS='\t' '{$5=$2 + 1; print $1 "\t" $2 "\t" $5 "\t" $4 "\t" $3}' \
> ${outputdir}/${sampleid}.reverse_NDR.bed

# Sort your generated BED files
sort -k 1V,1 -k 2n,2 ${outputdir}/${sampleid}.forward_NDR.bed -o ${outputdir}/${sampleid}.sortedNuc.forward_NDR.bed
sort -k 1V,1 -k 2n,2 ${outputdir}/${sampleid}.reverse_NDR.bed -o ${outputdir}/${sampleid}.sortedNuc.reverse_NDR.bed

echo -e "generating NDR REF 1 ..."
bedtools closest -a ${outputdir}/${sampleid}.sortedNuc.forward_NDR.bed -b ${ndr_ref1} -t first | awk -v OFS='\t' '{$9=$7 - $2;print $0}' > ${outputdir}/${sampleid}.NDR1_forward.dist.bed
bedtools closest -a ${outputdir}/${sampleid}.sortedNuc.reverse_NDR.bed -b ${ndr_ref1} -t first | awk -v OFS='\t' '{$9=$7 - $2;print $0}' > ${outputdir}/${sampleid}.NDR1_reverse.dist.bed

echo -e "sorting forward NDR REF 1 file"
sort -k4,4 ${outputdir}/${sampleid}.NDR1_forward.dist.bed > ${outputdir}/${sampleid}.NDR1_forward.dist.sorted.bed
echo -e "sorting reverse NDR REF 1 file"
sort -k4,4 ${outputdir}/${sampleid}.NDR1_reverse.dist.bed > ${outputdir}/${sampleid}.NDR1_reverse.dist.sorted.bed

echo -e "generating NRD REF 2 ..."
bedtools closest -a ${outputdir}/${sampleid}.sortedNuc.forward_NDR.bed -b ${ndr_ref2} -t first | awk -v OFS='\t' '{$9=$7 - $2;print $0}' > ${outputdir}/${sampleid}.NDR2_forward.dist.bed
bedtools closest -a ${outputdir}/${sampleid}.sortedNuc.reverse_NDR.bed -b ${ndr_ref2} -t first | awk -v OFS='\t' '{$9=$7 - $2;print $0}' > ${outputdir}/${sampleid}.NDR2_reverse.dist.bed

echo -e "sorting forward NDR REF 2 file"
sort -k4,4 ${outputdir}/${sampleid}.NDR2_forward.dist.bed > ${outputdir}/${sampleid}.NDR2_forward.dist.sorted.bed
echo -e "sorting reverse NDR REF 2 file"
sort -k4,4 ${outputdir}/${sampleid}.NDR2_reverse.dist.bed > ${outputdir}/${sampleid}.NDR2_reverse.dist.sorted.bed

echo -e "generating NRD REF 3 ..."
bedtools closest -a ${outputdir}/${sampleid}.sortedNuc.forward_NDR.bed -b ${ndr_ref2} -t first | awk -v OFS='\t' '{$9=$7 - $2;print $0}' > ${outputdir}/${sampleid}.NDR3_forward.dist.bed
bedtools closest -a ${outputdir}/${sampleid}.sortedNuc.reverse_NDR.bed -b ${ndr_ref2} -t first | awk -v OFS='\t' '{$9=$7 - $2;print $0}' > ${outputdir}/${sampleid}.NDR3_reverse.dist.bed

echo -e "sorting forward NDR REF 3 file"
sort -k4,4 ${outputdir}/${sampleid}.NDR3_forward.dist.bed > ${outputdir}/${sampleid}.NDR3_forward.dist.sorted.bed
echo -e "sorting reverse NDR REF 3 file"
sort -k4,4 ${outputdir}/${sampleid}.NDR3_reverse.dist.bed > ${outputdir}/${sampleid}.NDR3_reverse.dist.sorted.bed

##### column $10: distance of forward read to the nearest NDR
cat ${outputdir}/${sampleid}.NDR1_forward.dist.sorted.bed | cut -f9 > ${outputdir}/forward_ndr1.tmp.txt
paste ${finalOutputFile} ${outputdir}/forward_ndr1.tmp.txt  > ${outputdir}/${sampleid}.modified1.tsv

##### column $11: distance of reverse read to the nearest NDR
cat ${outputdir}/${sampleid}.NDR1_reverse.dist.sorted.bed | cut -f9 > ${outputdir}/reverse_ndr1.tmp.txt
paste ${outputdir}/${sampleid}.modified1.tsv ${outputdir}/reverse_ndr1.tmp.txt  > ${outputdir}/${sampleid}.modified2.tsv

##### column $10: distance of forward read to the nearest NDR
cat ${outputdir}/${sampleid}.NDR2_forward.dist.sorted.bed | cut -f9 > ${outputdir}/forward_ndr2.tmp.txt
paste ${outputdir}/${sampleid}.modified2.tsv ${outputdir}/forward_ndr2.tmp.txt  > ${outputdir}/${sampleid}.modified3.tsv

##### column $11: distance of reverse read to the nearest NDR
cat ${outputdir}/${sampleid}.NDR2_reverse.dist.sorted.bed | cut -f9 > ${outputdir}/reverse_ndr2.tmp.txt
paste ${outputdir}/${sampleid}.modified3.tsv ${outputdir}/reverse_ndr2.tmp.txt  > ${outputdir}/${sampleid}.modified4.tsv

##### column $10: distance of forward read to the nearest NDR
cat ${outputdir}/${sampleid}.NDR2_forward.dist.sorted.bed | cut -f9 > ${outputdir}/forward_ndr3.tmp.txt
paste ${outputdir}/${sampleid}.modified4.tsv ${outputdir}/forward_ndr3.tmp.txt  > ${outputdir}/${sampleid}.modified5.tsv

##### column $11: distance of reverse read to the nearest NDR
cat ${outputdir}/${sampleid}.NDR2_reverse.dist.sorted.bed | cut -f9 > ${outputdir}/reverse_ndr3.tmp.txt
paste ${outputdir}/${sampleid}.modified5.tsv ${outputdir}/reverse_ndr3.tmp.txt  > ${outputdir}/${sampleid}.modified6.tsv


mv ${outputdir}/${sampleid}.modified6.tsv ${outputdir}/${sampleid}.final_output.tsv
rm -rf ${outputdir}/${sampleid}.modified{1,2,3,4,5,6}.tsv
rm -rf ${outputdir}/*ndr*.tmp.txt

##### example cmd
# finalOutputFile="/media/hieunguyen/HNSD01/src/ecd_wgs_features/testRUN/WGShg19/WGShg19.final_output.tsv";
# bedfile="/media/hieunguyen/HNSD01/src/ecd_wgs_features/bin1M.sorted.bed";
# outputdir="/media/hieunguyen/HNSD01/src/ecd_wgs_features/testRUN/WGShg19/NDR_output";
# ndr_ref1="/media/hieunguyen/HNSD01/src/ecd_wgs_features/NDR_bed_files/NDR_cancer_specific_location_10cancer-types.csv";
# ndr_ref1="/media/hieunguyen/HNSD01/src/ecd_wgs_features/NDR_bed_files/NDR_cancer_specific_location_SingleCancer-CRC-1.csv";
# ndr_ref1="/media/hieunguyen/HNSD01/src/ecd_wgs_features/NDR_bed_files/NDR_cancer_specific_location_SingleCancer-LUNG.csv";
# bash 07_add_NDR_feature_to_current_final_output.sh -i ${finalOutputFile} -o ${outputdir} -n ${ndr_ref1} -d ${ndr_ref2} -r ${ndr_ref3}
