export PATH=/Users/hieunguyen/samtools/bin:$PATH 
# bash 01_BAM_to_FRAGS.sh -i /Users/hieunguyen/data/tmp/debug_ecd_wgs_features/data/1-ZLAAO90NB_S7509-S7709.sorted.bam -o ./output_debug/ -n 10 -f /Users/hieunguyen/data/resources/hg19.fa 
#####----------------------------------------------------------------------#####
##### input args
#####----------------------------------------------------------------------#####
while getopts "i:o:n:t:" opt; do
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

sampleid=$(echo ${inputbam} | xargs -n 1 basename)
sampleid=${sampleid%.bam*}
outputdir=${outputdir}/${sampleid}

##### Split bam into short and long bam
echo -e "splitting BAM file into short and long BAM files ..."

samtools view -f 2 -h ${inputbam} | \
  awk 'substr($0,1,1)=="@" || ($9 >= 50 && $9 <= 150) || ($9 <= -50 && $9 >= -150)' | \
  samtools view -b > ${outputdir}/${sampleid}_50_150.short.bam;
samtools index ${outputdir}/${sampleid}_50_150.short.bam
samtools view -f 2 -h ${inputbam} | \
  awk 'substr($0,1,1)=="@" || ($9 > 150 && $9 <= 250) || ($9 < -150 && $9 >= -250)' | \
  samtools view -b > ${outputdir}/${sampleid}_151_250.long.bam;
samtools index ${outputdir}/${sampleid}_151_250.long.bam;

##### Split bam into short and long bam with new cutoff
samtools view -f 2 -h ${inputbam} | \
  awk 'substr($0,1,1)=="@" || ($9 >= 100 && $9 <= 150) || ($9 <= -100 && $9 >= -150)' | \
  samtools view -b > ${outputdir}/${sampleid}.100_150.short.bam;
samtools index ${outputdir}/${sampleid}.100_150.short.bam;
samtools view -f 2 -h ${inputbam} | \
  awk 'substr($0,1,1)=="@" || ($9 > 151 && $9 <= 220) || ($9 < -151 && $9 >= -220)' | \
  samtools view -b > ${outputdir}/${sampleid}.151_220.long.bam;
samtools index ${outputdir}/${sampleid}.151_220.long.bam;

