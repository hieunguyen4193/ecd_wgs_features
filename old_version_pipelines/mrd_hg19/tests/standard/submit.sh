#!/usr/bin/bash

# activate nextflow environment
source activate /mnt/DATASSD/PIPELINESRC/OncoGS_pipeline/nextflow22

# the nextflow command
nextflow run ../.. -profile docker -resume \
    --fqdir data 