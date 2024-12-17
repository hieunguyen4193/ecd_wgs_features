echo "START time: $(date)";
export PATH=/home/hieunguyen/bedtools2/bin/:$PATH
export PATH=/home/hieunguyen/samtools/bin/:$PATH

while getopts "i:o:b:" opt; do
  case ${opt} in
    i )
      finalOutputFile=$OPTARG
      ;;
    o )
      outputdir=$OPTARG
      ;;
    b )
      bedfile=$OPTARG
      ;;
    \? )
      echo "Usage: cmd [-i] finalOutputFile [-o] outputdir [-b] bed file CNA bins"
      exit 1
      ;;
  esac
done
# if bed file is not sorted. 
# sort -k 1V,1 -k 2n,2 bin1M.bed -o bin1M.sorted.bed

sampleid=$(echo ${finalOutputFile} | xargs -n 1 basename)
sampleid=${sampleid%.final_output.tsv*}
outputdir=${outputdir}/${sampleid}
mkdir -p ${outputdir};

cat ${finalOutputFile} | cut -f1,2,3 | awk '{print $0"\t"NR}' > ${outputdir}/coordCNA.final_output.bed
bedtools intersect -a ${outputdir}/coordCNA.final_output.bed -b ${bedfile} -wa -c | cut -f5 > ${outputdir}/coordCNA.final_output.count.bed

paste ${outputdir}/coordCNA.final_output.bed ${outputdir}/coordCNA.final_output.count.bed > ${outputdir}/coordCNA.final_output.addCount.bed

cat ${outputdir}/coordCNA.final_output.addCount.bed | \
    awk -v OFS="\t" '{if($5 == 1) print $0}' > ${outputdir}/coordCNA.final_output.UniqueCount.bed
    
bedtools intersect -a ${outputdir}/coordCNA.final_output.UniqueCount.bed -b ${bedfile} -wa -wb -loj | cut -f4,9 > ${outputdir}/coordCNA.final_output.intersect.bed

echo -e "start joining table...";
echo "START JOINNING TABLE time: $(date)";
join -1 4 -2 1 -t $'\t' -a1 ${outputdir}/coordCNA.final_output.bed ${outputdir}/coordCNA.final_output.intersect.bed | cut -f5 > ${outputdir}/coordCNA.final.tsv;

paste ${finalOutputFile} ${outputdir}/coordCNA.final.tsv > ${outputdir}/${sampleid}.final_output.addCNA.tsv ;
rm -rf ${outputdir}/coordCNA*
echo "FINAL END time: $(date)";

##### example cmd
# finalOutputFile="/media/hieunguyen/HNSD01/src/ecd_wgs_features/testRUN/WGShg19/WGShg19.final_output.tsv";
# bedfile="/media/hieunguyen/HNSD01/src/ecd_wgs_features/bin1M.sorted.bed";
# outputdir="/media/hieunguyen/HNSD01/src/ecd_wgs_features/testRUN/WGShg19/intersect_output";
# bash 08_add_CNA_bin_to_fragments.sh -i ${finalOutputFile} -o ${outputdir} -b ${bedfile}

