#!/bin/bash
# inputbam=$1
# outputdir=$2
# num_threads=$3

# get samtools to work
export PATH=/Users/hieunguyen/samtools/bin:$PATH 
# samtools --version

# default input values
inputbam="/Volumes/HNSD02/data/WGS_bam/9-ZMC014NB_S95025-S97025.sorted.bam";
outputdir="/Volumes/HNSD02/outdir/ecd_wgs_features";

mkdir -p ${outputdir};

# Parse command line arguments
while getopts "i:o:t:" opt; do
  case ${opt} in
    i )
      inputbam=$OPTARG
      ;;
    o )
      outputdir=$OPTARG
      ;;
    \? )
      echo "Usage: cmd [-i] inputbam [-o] outputdir [-t] num_threads"
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

filename=$(echo ${inputbam} | xargs -n 1 basename)
filename=${filename%.bam*}

samtools view ${inputbam}| cut -f1,3,4,9 \
| awk -v OFS='\t' '{if ($4 > 0){$5=$3+$4; print $0}}' > ${outputdir}/${filename}.coord.bed
    
