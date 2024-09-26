#!/bin/bash
#SBATCH --job-name="ECD_WGS_features_HG38"
#SBATCH --partition=research_hi
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH --output=/mnt/DATASM14/hieunho/hieu_gitlab/ecd_wgs_features_hg38/log/log.log

##########################################################################################
INDIR=/mnt/DATASM14/hieunho/hieu_pipeline/ECD_WGS_features_HG38/test_data
OUTDIR=/datassd/hieunho/ECD_WGS_features_HG38/test_data
submit_type={{submit_type}}
working_dir={{working_dir}}
BWA_REPO={{BWA_REPO}}
RESOURCE={{resource}}
##########################################################################################
mkdir -p $OUTDIR
cd $OUTDIR

nextflow_dir=${working_dir}/nextflow
SRC=${working_dir}/src
work_dir=$OUTDIR/work
report_file=$OUTDIR/report.html
work=$OUTDIR/work

nextflow run $nextflow_dir/main_bwa.nf \
    -c $nextflow_dir/main_${submit_type}.config \
    --OUTDIR $OUTDIR \
    --FQDIR $INDIR \
    --PP_REPO $working_dir \
    --BWA_REPO $BWA_REPO \
    -w $work \
    --SRC $SRC \
    --RESOURCE $RESOURCE -with-report $report_file -resume