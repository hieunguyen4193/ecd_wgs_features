export PATH=/Users/hieunguyen/samtools/bin:$PATH 
# bash 01_BAM_to_FRAGS.sh -i /Users/hieunguyen/data/tmp/debug_ecd_wgs_features/data/1-ZLAAO90NB_S7509-S7709.sorted.bam -o ./output_debug/ -n 10 -f /Users/hieunguyen/data/resources/hg19.fa 
#####----------------------------------------------------------------------#####
##### input args
#####----------------------------------------------------------------------#####
while getopts "i:o:n:b:" opt; do
  case ${opt} in
    i )
      finalOutputFile=$OPTARG
      ;;
    o )
      outputdir=$OPTARG
      ;;
    n )
      ndr_ref=$OPTARG
      ;;
    b )
      ndrb_ref=$OPTARG
      ;;

    \? )
      echo "Usage: cmd [-i] finalOutputFile [-o] outputdir [-n] NDR ref [-b] NDR binary ref"
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

echo -e "generating NRD ..."
bedtools closest -a ${outputdir}/${sampleid}.sortedNuc.forward_NDR.bed -b ${ndr_ref} -t first | awk -v OFS='\t' '{$9=$7 - $2;print $0}' > ${outputdir}/${sampleid}.NDR_forward.dist.bed
bedtools closest -a ${outputdir}/${sampleid}.sortedNuc.reverse_NDR.bed -b ${ndr_ref} -t first | awk -v OFS='\t' '{$9=$7 - $2;print $0}' > ${outputdir}/${sampleid}.NDR_reverse.dist.bed

echo -e "sorting forward NDR file"
sort -k4,4 ${outputdir}/${sampleid}.NDR_forward.dist.bed > ${outputdir}/${sampleid}.NDR_forward.dist.sorted.bed
echo -e "sorting reverse NDR file"
sort -k4,4 ${outputdir}/${sampleid}.NDR_reverse.dist.bed > ${outputdir}/${sampleid}.NDR_reverse.dist.sorted.bed

echo -e "generating NRD BINARY ..."
bedtools closest -a ${outputdir}/${sampleid}.sortedNuc.forward_NDR.bed -b ${ndrb_ref} -t first | awk -v OFS='\t' '{$9=$7 - $2;print $0}' > ${outputdir}/${sampleid}.NDRb_forward.dist.bed
bedtools closest -a ${outputdir}/${sampleid}.sortedNuc.reverse_NDR.bed -b ${ndrb_ref} -t first | awk -v OFS='\t' '{$9=$7 - $2;print $0}' > ${outputdir}/${sampleid}.NDRb_reverse.dist.bed

echo -e "sorting forward NDR file"
sort -k4,4 ${outputdir}/${sampleid}.NDRb_forward.dist.bed > ${outputdir}/${sampleid}.NDRb_forward.dist.sorted.bed
echo -e "sorting reverse NDR file"
sort -k4,4 ${outputdir}/${sampleid}.NDRb_reverse.dist.bed > ${outputdir}/${sampleid}.NDRb_reverse.dist.sorted.bed

##### column $10: distance of forward read to the nearest NDR
cat ${outputdir}/${sampleid}.NDR_forward.dist.sorted.bed | cut -f9 > ${outputdir}/forward_ndr.tmp.txt
paste ${finalOutputFile} ${outputdir}/forward_ndr.tmp.txt  > ${outputdir}/${sampleid}.modified1.tsv

##### column $11: distance of reverse read to the nearest NDR
cat ${outputdir}/${sampleid}.NDR_reverse.dist.sorted.bed | cut -f9 > ${outputdir}/reverse_ndr.tmp.txt
paste ${outputdir}/${sampleid}.modified1.tsv ${outputdir}/reverse_ndr.tmp.txt  > ${outputdir}/${sampleid}.modified2.tsv

##### column $10: distance of forward read to the nearest NDR
cat ${outputdir}/${sampleid}.NDRb_forward.dist.sorted.bed | cut -f9 > ${outputdir}/forward_ndrb.tmp.txt
paste ${outputdir}/${sampleid}.modified2.tsv ${outputdir}/forward_ndrb.tmp.txt  > ${outputdir}/${sampleid}.modified3.tsv

##### column $11: distance of reverse read to the nearest NDR
cat ${outputdir}/${sampleid}.NDRb_reverse.dist.sorted.bed | cut -f9 > ${outputdir}/reverse_ndrb.tmp.txt
paste ${outputdir}/${sampleid}.modified3.tsv ${outputdir}/reverse_ndrb.tmp.txt  > ${outputdir}/${sampleid}.modified4.tsv

mv ${outputdir}/${sampleid}.modified4.tsv ${outputdir}/${sampleid}.final_output.tsv
rm -rf ${outputdir}/${sampleid}.modified{1,2,3,4}.tsv