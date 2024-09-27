# export PATH=/home/hieunguyen/samtools/bin:$PATH
export PATH=/Users/hieunguyen/samtools/bin:$PATH

# note for HieuNguyen: on your macM3 do not use the conda env pytorch, use the base env. It has java 23.
# bash generate_QC_metrics.sh -i /Users/hieunguyen/src/ecd_wgs_features/output/WGShg19.sorted.bam -o ./output -n 10 -f /Users/hieunguyen/data/resources/hg19.fa -m true

#####----------------------------------------------------------------------#####
##### input args
#####----------------------------------------------------------------------#####
while getopts "i:o:n:f:m:" opt; do
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
    m )
      rerun_markdup=$OPTARG
      ;;
    \? )
      echo "Usage: cmd [-i] inputbam [-o] outputdir [-n] samtools_num_threads [-f] reference genome [-m] rerun_markdup"
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
mkdir -p ${outputdir}/QC/${sampleid};

#####----------------------------------------------------#####
##### Perform mark duplication
#####----------------------------------------------------#####
if [ "${rerun_markdup}" = "true" ]; then
    echo "Run PICARD MARK DUPLICATE for the input bam file" $inputbam
    java -Xms512m -Xmx4g -jar ./picard.jar MarkDuplicates \
    I=${inputbam} \
    O=${outputdir}/${sampleid}.sorted.markdup.bam \
    M=${outputdir}/QC/${sampleid}/${sampleid}.marked_dup_metrics.txt

    # replace inputbam by the new markdup bam
    inputbam=${outputdir}/${sampleid}.sorted.markdup.bam
    echo -e "New markdup bam file: " ${inputbam}
    echo -e "Finished mark duplication"
elif [ "${rerun_markdup}" = "false" ]; then
    echo "Assume the input file is alread MARK DUPLICATED ..."
    # Add your code here for when the flag is false
else
    echo "Invalid value for boolean flag: ${rerun_markdup}. Use 'true' or 'false' only."
    exit 1
fi

#####----------------------------------------------------#####
##### generate QC metrics for the input BAM file. 
#####----------------------------------------------------#####

# samtools flagstat
echo -e "runnning samtools flagstats"
samtools flagstat ${inputbam} -@ ${samtools_num_threads} --output-fmt 'tsv' > ${outputdir}/QC/${sampleid}/${sampleid}.flagstat.txt
echo -e "finished"

# samtools 
echo -e "running samtools stats"
samtools stat ${inputbam} -@ ${samtools_num_threads} > ${outputdir}/QC/${sampleid}/${sampleid}.sortedbamfilename.stats
plot-bamstats -p ${outputdir}/QC/${sampleid}/samtools_stat_plots/${sampleid} ${outputdir}/QC/${sampleid}/${sampleid}.sortedbamfilename.stats
rm -rf ${outputdir}/QC/${sampleid}/samtools_stat_plot/*.gp
echo -e "finished"

# Some picard QC metrics
echo -e "running picard collect alignment summary metrics"
java -Xms512m -Xmx4g -jar ./picard.jar CollectAlignmentSummaryMetrics \
    R=${ref} \
    I=${inputbam} \
    O=${outputdir}/QC/${sampleid}/${sampleid}.alignment_summary_metrics.txt
sed '/^#/d' ${outputdir}/QC/${sampleid}/${sampleid}.alignment_summary_metrics.txt > ${outputdir}/QC/${sampleid}/${sampleid}.alignment_summary_metrics.final.txt
echo -e "finished"

echo -e "running collect insert size metrics";
java -Xms512m -Xmx4g -jar ./picard.jar CollectInsertSizeMetrics \
    I=${inputbam} \
    O=${outputdir}/QC/${sampleid}/${sampleid}.insert_size_metrics.txt \
    H=${outputdir}/QC/${sampleid}/${sampleid}.insert_size_histogram.pdf

sed '/^#/d' ${outputdir}/QC/${sampleid}/${sampleid}.insert_size_metrics.txt > ${outputdir}/QC/${sampleid}/${sampleid}.insert_size_metrics.final.txt
echo -e "finished"

echo -e "running collect gc bias metrics";″
java -Xms512m -Xmx4g -jar ./picard.jar CollectGcBiasMetrics \
    R=${ref} \
    I=${inputbam} \
    O=${outputdir}/QC/${sampleid}/${sampleid}.gc_bias_metrics.txt \
    CHART=${outputdir}/QC/${sampleid}/${sampleid}.gc_bias_metrics.pdf \
    S=${outputdir}/QC/${sampleid}/${sampleid}.summary_metrics.txt

sed '/^#/d' ${outputdir}/QC/${sampleid}/${sampleid}.summary_metrics.txt > ${outputdir}/QC/${sampleid}/${sampleid}.summary_metrics.final.txt
echo -e "finished"

echo -e "running collect gc bias metrics";″
java -Xms512m -Xmx4g -jar ./picard.jar CollectWgsMetrics \
    R=${ref} \
    I=${inputbam} \
    O=${outputdir}/QC/${sampleid}/${sampleid}.wgs_metrics.txt

sed '/^#/d' ${outputdir}/QC/${sampleid}/${sampleid}.wgs_metrics.txt > ${outputdir}/QC/${sampleid}/${sampleid}.wgs_metrics.final.txt
echo -e "finished"

