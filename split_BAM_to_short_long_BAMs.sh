
export PATH=/home/hieunguyen/samtools/bin:$PATH
# bash split_BAM_to_short_long_BAMs.sh -i /media/hieunguyen/GSHD_HN01/raw_data/bam_files/WGShg19.bam  -o ./output/

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

export PATH=/home/hieunguyen/samtools/bin:$PATH

samtools view -f 2 -h ${inputbam} | \
  awk 'substr($0,1,1)=="@" || ($9 >= 100 && $9 <= 150) || ($9 <= -100 && $9 >= -150)' | \
  samtools view -b > ${outputdir}/${sampleid}.short.bam;
samtools index ${outputdir}/${sampleid}.short.bam
samtools view -f 2 -h ${inputbam} | \
  awk 'substr($0,1,1)=="@" || ($9 > 151 && $9 <= 220) || ($9 < -151 && $9 >= -220)' | \
  samtools view -b > ${outputdir}/${sampleid}.long.bam;
samtools index ${outputdir}/${sampleid}.long.bam;