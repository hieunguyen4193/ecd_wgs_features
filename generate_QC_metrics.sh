# export PATH=/home/hieunguyen/samtools/bin:$PATH
export PATH=/Users/hieunguyen/samtools/bin:$PATH

# note for HieuNguyen: on your macM3 do not use the conda env pytorch, use the base env. It has java 23.
# bash generate_QC_metrics.sh -i /media/hieunguyen/GSHD_HN01/raw_data/bam_files/WGShg19.bam  -o ./output/ -n 40 -f /media/hieunguyen/HNSD01/resources/hg19.fa
# bash generate_QC_metrics.sh -i /Users/hieunguyen/src/ecd_wgs_features/output/WGShg19.sorted.bam -o ./output -n 10 -f /Users/hieunguyen/data/resources/hg19.fa

#####----------------------------------------------------------------------#####
##### input args
#####----------------------------------------------------------------------#####
while getopts "i:o:n:f:" opt; do
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
    f )
      ref=$OPTARG
      ;;
    \? )
      echo "Usage: cmd [-i] inputbam [-o] outputdir [-n] samtools_num_threads [-f] reference genome"
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

#####----------------------------------------------------#####
##### generate QC metrics for the input BAM file. 
#####----------------------------------------------------#####

# samtools flagstat
mkdir -p ${outputdir}/QC;
echo -e "runnning samtools flagstats"
samtools flagstat ${inputbam} -@ ${samtools_num_threads} --output-fmt 'tsv' > ${outputdir}/QC/${sampleid}.flagstat.txt
echo -e "finished"

# samtools 
echo -e "running samtools stats"
samtools stat ${inputbam} -@ ${samtools_num_threads} > ${outputdir}/QC/${sampleid}.sortedbamfilename.stats
plot-bamstats -p ${outputdir}/QC/samtools_stat_plots/${sampleid} ${outputdir}/QC/${sampleid}.sortedbamfilename.stats
rm -rf ${outputdir}/QC/samtools_stat_plot/*.gp
echo -e "finished"

# Some picard QC metrics
echo -e "running picard collect alignment summary metrics"
java -Xms512m -Xmx4g -jar ./picard.jar CollectAlignmentSummaryMetrics \
    R=${ref} \
    I=${inputbam} \
    O=${outputdir}/QC/${sampleid}.alignment_summary_metrics.txt
echo -e "finished"

echo -e "running collect insert size metrics";
java -Xms512m -Xmx4g -jar ./picard.jar CollectInsertSizeMetrics \
    I=${inputbam} \
    O=${outputdir}/QC/${sampleid}.insert_size_metrics.txt \
    H=${outputdir}/QC/${sampleid}.insert_size_histogram.pdf
echo -e "finished"

echo -e "running collect gc bias metrics";″
java -Xms512m -Xmx4g -jar ./picard.jar CollectGcBiasMetrics \
    R=${ref} \
    I=${inputbam} \
    O=${outputdir}/QC/${sampleid}.gc_bias_metrics.txt \
    CHART=${outputdir}/QC/${sampleid}.gc_bias_metrics.pdf \
    S=${outputdir}/QC/${sampleid}.summary_metrics.txt
echo -e "finished"

echo -e "running collect gc bias metrics";″
java -Xms512m -Xmx4g -jar ./picard.jar CollectWgsMetrics \
    R=${ref} \
    I=${inputbam} \
    O=${outputdir}/QC/${sampleid}.gc_bias_metrics.txt
echo -e "finished"

