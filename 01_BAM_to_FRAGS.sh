#####-----------------------------------------------#####
##### INTRODUCTION
#####-----------------------------------------------#####
# This script pre-process an input BAM file to a 
# fragment-wise data features, which can be use to calculate
# several fragmentomics features. 
export PATH=/Users/hieunguyen/samtools/bin:$PATH

#####-----------------------------------------------#####
##### input args
#####-----------------------------------------------#####
while getopts "i:o:t:r:" opt; do
  case ${opt} in
    i )
      inputbam=$OPTARG
      ;;
    o )
      outputdir=$OPTARG
      ;;
    \? )
      echo "Usage: cmd [-i] inputbam [-o] outputdir"
      exit 1
      ;;
  esac
done


echo -e "input bam file: " ${inputbam}
# Check if the input BAM file exists
if [ ! -f "${inputbam}" ]; then
    echo "Input BAM file does not exist: ${inputbam}"
    exit 1
fi

sampleid=$(echo ${inputbam} | xargs -n 1 basename)
sampleid=${sampleid%.bam*}

#####-----------------------------------------------#####
##### pre-processing steps for an input bam file:
#####-----------------------------------------------#####
samtools view -f 3 ${inputbam} -@ ${num_threads} -b -o ${outputdir}/${sampleid}.prep.bam
samtools index -@ ${num_threads} ${outputdir}/${sampleid}.prep.bam
samtools view ${outputdir}/${sampleid}.prep.bam | cut -f1,3,4,6,9 > ${outputdir}/${sampleid}.prep.tsv

##### get true fragment end
echo -e "modify fragment end, fragment end = foward read start + fragment length"

##### explain new columns in the data:
# $3: start of read 1 (fragment start)
# $6: start of read 2 (fragment end)
# $7 = $3 + 1: fragment start + 1, to use in nucleosome footprint feature
# $8 = $6 + 1: fragment end + 1, to use in nucleosome footprint feature

##### explain the sign of TLEN field $9
# The TLEN field is positive for the leftmost segment of the template, 
# negative for the rightmost, and the sign for any middle segment is undefined. 
# If segments cover the same coordinates then the choice of which is leftmost 
# and rightmost is arbitrary, but the two ends must still have differing signs
# https://samtools.github.io/hts-specs/SAMv1.pdf
awk -v OFS='\t' '{if ($5 > 0){$6=$3+$5; $7=$3+1; $8=$6+1; $9=$1"_"$2"_"$3; print $0}}' \
  ${outputdir}/${sampleid}.prep.tsv \
  | sort -k9,9 > ${outputdir}/${sampleid}.modified.tsv


