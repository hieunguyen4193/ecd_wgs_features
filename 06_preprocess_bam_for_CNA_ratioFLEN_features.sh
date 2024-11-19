#####----------------------------------------------------------------------#####
##### INTRODUCTION
#####----------------------------------------------------------------------#####
# This script pre-process an input BAM file to a 
# fragment-wise data features, which can be use to calculate
# several fragmentomics features. 
# export PATH=/Users/hieunguyen/samtools/bin:$PATH
# bash 01_BAM_to_FRAGS.sh -i /Users/hieunguyen/data/bam_files/WGShg19.bam  -o ./output/ -n 10

# export PATH=/home/hieunguyen/samtools/bin:$PATH
# bash 01_BAM_to_FRAGS.sh -i /media/hieunguyen/GSHD_HN01/raw_data/bam_files/WGShg19.bam  -o ./output/ -n 40

export PATH=/Users/hieunguyen/samtools/bin:$PATH 
# bash 01_BAM_to_FRAGS.sh -i /Users/hieunguyen/data/tmp/debug_ecd_wgs_features/data/1-ZLAAO90NB_S7509-S7709.sorted.bam -o ./output_debug/ -n 10 -f /Users/hieunguyen/data/resources/hg19.fa 
#####----------------------------------------------------------------------#####
##### input args
#####----------------------------------------------------------------------#####
while getopts "i:o:n:t" opt; do
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
    t )
      input_type=$OPTARG
      ;;
    
    \? )
      echo "Usage: cmd [-i] inputbam [-o] outputdir [-n] samtools_num_threads [-t] input_type"
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
outputdir=${outputdir}/${sampleid}

mkdir -p ${outputdir}

if [ "${input_type}" = "markdup" ]; then
# if the input BAM file is already PREPROCESSED and MARKDUP
  echo "Processing as markdup BAM file"
  bash split_bam_short_long.sh -i ${inputbam} -o ${outputdir} -n ${samtools_num_threads}
else
  echo "Processing as RAW BAM file"
  echo -e "remove unpaired and unmapped reads in BAM files, generate prep.tsv file";
  samtools view -h -f 3 ${inputbam} | samtools sort -n -@ ${samtools_num_threads} -o ${outputdir}/tmp.bam;
  samtools view -h ${outputdir}/tmp.bam | awk -f preprocessing_script.awk - > ${outputdir}/tmp.sam;
  samtools sort -@ ${samtools_num_threads} -O BAM -o ${outputdir}/${sampleid}.tmp.sorted.bam ${outputdir}/tmp.sam;

  ##### mark duplicates
  java -Xms512m -Xmx4g -jar ./picard.jar MarkDuplicates \
      I=${outputdir}/${sampleid}.tmp.sorted.bam \
      O=${outputdir}/${sampleid}.sorted.markdup.bam \
      M=${outputdir}/${sampleid}.marked_dup_metrics.txt

  rm -rf ${outputdir}/${sampleid}.tmp.sorted.bam
  samtools index -@ ${samtools_num_threads} ${outputdir}/${sampleid}.sorted.markdup.bam
  
  bash split_bam_short_long.sh -i ${outputdir}/${sampleid}.sorted.markdup.bam -o ${outputdir} -n ${samtools_num_threads}
fi