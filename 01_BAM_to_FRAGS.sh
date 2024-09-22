#####----------------------------------------------------------------------#####
##### INTRODUCTION
#####----------------------------------------------------------------------#####
# This script pre-process an input BAM file to a 
# fragment-wise data features, which can be use to calculate
# several fragmentomics features. 
export PATH=/Users/hieunguyen/samtools/bin:$PATH

# bash 01_BAM_to_FRAGS.sh -i /Users/hieunguyen/data/bam_files/WGShg19.bam  -o ./output/ -n 10

#####----------------------------------------------------------------------#####
##### input args
#####----------------------------------------------------------------------#####
while getopts "i:o:n:" opt; do
  case ${opt} in
    i )
      inputbam=$OPTARG
      ;;
    o )
      outputdir=$OPTARG
      ;;
    n )
      samtools_num_threads=$OPTARG
      ;;
    \? )
      echo "Usage: cmd [-i] inputbam [-o] outputdir [-n] samtools_num_threads"
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

mkdir -p ${outputdir}
sampleid=$(echo ${inputbam} | xargs -n 1 basename)
sampleid=${sampleid%.bam*}

#####----------------------------------------------------------------------#####
##### pre-processing steps for an input bam file:
#####----------------------------------------------------------------------#####
# remove unpaired and unmapped reads in BAM files.
# input bam file must be sorted and index
# samtools sort input.bam > input.sorted.bam
# samtools index input.sorted.bam

if [ ! -f "${outputdir}/${sampleid}.frag.tsv" ]; then
  echo -e "Pre-processing BAM files and convert them to fragmentomics format ..."

  if [ ! -f "${outputdir}/${sampleid}.prep.tsv" ]; then
  echo -e "remove unpaired and unmapped reads in BAM files, generate prep.tsv file";
  samtools view -h -f 3 ${inputbam} | samtools sort -n -@ 5 -o tmp.bam;
  samtools view -h tmp.bam | awk -f preprocessing_script.awk - > tmp.sam;
  samtools sort -@ 5 -O BAM -o ${outputdir}/${sampleid}.sorted.bam tmp.sam;
  samtools index -@ 5 ${outputdir}/${sampleid}.sorted.bam

  samtools view ${outputdir}/${sampleid}.sorted.bam | cut -f1,3,4,6,9 > ${outputdir}/${sampleid}.prep.tsv
  fi

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
  # this remove half of the reads, keep only fragment-wise information. 
  # 1 line = 1 fragment, chrom - start - end and a read ID, use this read ID to sort and
  # match the information later when calculating End motif and nucleosome footprint. 
  awk -v OFS='\t' '{if ($5 > 0){$6=$3+$5; $7=$1"_"$2"_"$3; print $0}}' \
    ${outputdir}/${sampleid}.prep.tsv \
    | sort -k9,9 \
    | awk '{ print $2 "\t" $3 "\t" $6 "\t" $5 "\t" $7}' > ${outputdir}/${sampleid}.frag.tsv
  rm -rf ${outputdir}/${sampleid}.prep.tsv
fi


