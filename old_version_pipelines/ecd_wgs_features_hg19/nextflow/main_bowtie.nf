nextflow.enable.dsl = 1

ADAPTER_FILE    = file("$params.PP_REPO/resource/adapter/TruSeq3-PE-2.fa")
INDEX_FILE      = file("$params.PP_REPO/resource/bt2")

params.SRC=""
SRC=file(params.SRC)

params.preprocessing_script="$SRC/preprocessing_script.awk"
preprocessing_script=file(params.preprocessing_script)

params.RESOURCE=""
RESOURCE=file(params.RESOURCE)

params.OUTDIR=""
OUTDIR=file(params.OUTDIR)

params.hg_type="19"
params.bin_size="1000"
params.window_size="1000000"
params.ref_name="rpr_map_Budhraja_STM2023"

//--------------------------------------- PROCESS ------------------------------------------------//
Channel
    .fromFilePairs("$params.FQDIR/*_R{1,2}*")
    .ifEmpty { error "Cannot find any reads matching: ${params.FQDIR}/*_R{1,2}*" }
    .into {Fastqc_ch; Dedupfastq_ch}

process FastQC {
    cache "deep"; tag "$sample_id"
    storeDir "$params.OUTDIR/fastqc"
    
    input:
        tuple sample_id, file(READS) from Fastqc_ch
    output:
        tuple sample_id, "${sample_id}_R1*fastqc.html", "${sample_id}_R2*fastqc.html"
    script:
        "fastqc -t $params.FASTQC_THREAD --quiet --outdir \$(pwd) ${READS[0]} ${READS[1]}"
}

process Dedupfastq {
    cache "deep"; tag "$sample_id"
    storeDir "$params.OUTDIR/dedup"

    input:
        tuple sample_id, file(READS) from Dedupfastq_ch
    output:
        tuple sample_id, "${sample_id}_R*.dedup.fastq.gz" into Trim_ch, Trim_Wcx_ch
    script:
        """
        clumpify.sh -Xmx5g in=${READS[0]} in2=${READS[1]} \
        out=${sample_id}_R1.dedup.fastq.gz out2=${sample_id}_R2.dedup.fastq.gz dedupe subs=0
        """
}

process Trim {
    cache "deep"; tag "$sample_id"
    storeDir "$params.OUTDIR/trim"
    
    input:
        tuple sample_id, file(READS) from Trim_ch
        file ADAPTER_FILE
    output:
        tuple sample_id, "${sample_id}_R{1,2}_trim50.fastq.gz" into Align_ch   
    script:
       """
       trimmomatic PE -phred33 -threads $params.TRIM_THREAD ${READS[0]} ${READS[1]} \
       ${sample_id}_R1_trim50.fastq.gz ${sample_id}_R1_UP.fastq.gz \
       ${sample_id}_R2_trim50.fastq.gz ${sample_id}_R2_UP.fastq.gz \
       ILLUMINACLIP:${ADAPTER_FILE}:2:30:10 LEADING:3 TRAILING:3 \
       SLIDINGWINDOW:4:15 CROP:$params.READLENGTH MINLEN:36

       rm -rf ${sample_id}_R*_UP.fastq.gz
       """
}

process Align {
    cache "deep"; tag "$sample_id"
    storeDir "$params.OUTDIR/align_bowtie"
    
    input:
        tuple sample_id, file ("${sample_id}_R*_trim50.fastq.gz") from Align_ch
        file INDEX_FILE
    output:
        tuple sample_id, "${sample_id}.sam" into SortBam_ch
    
    script:
        """
        bowtie2 --threads $params.ALIGN_THREAD --mm \
        --very-fast \
        --no-mixed \
        --no-discordant \
        --no-unal \
        --rg-id ${sample_id} \
        --rg SM:${sample_id} \
        -x ${INDEX_FILE}/hg19 \
        -1 ${sample_id}_R1_trim50.fastq.gz -2 ${sample_id}_R2_trim50.fastq.gz \
        -S ${sample_id}.sam
        """
}

process SortBam {
    cache "deep"; tag "$sample_id"
    storeDir "$params.OUTDIR/align"
    
    input:
        tuple sample_id, file(BAM) from SortBam_ch
    output:
        tuple sample_id, file("${sample_id}.sorted.bam") into GW_feature_ch
        tuple sample_id, file("${sample_id}.sorted.bam"), file("${sample_id}.sorted.bam.bai") into IchorCNA_ch
    script:
        """
        samtools sort --threads $params.SAMTOOLSORT_THREAD -m 2G ${BAM[0]} -o ${sample_id}.sorted.bam
        samtools index ${sample_id}.sorted.bam 
        """
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process Preprocessing {
    cache "deep"; tag "$sample_id"
    // publishDir "${params.OUTDIR}/Preprocessing", mode: "copy"
    storeDir "${params.OUTDIR}/Preprocessing"

    input:
        tuple sample_id, file(Bam) from GW_feature_ch
        file SRC
        file RESOURCE
        file preprocessing_script
    output:
        tuple sample_id, file("$sample_id") into Calculate_feature_ch
    script:
    """
    bash $SRC/preprocess_bam_file.GW.sh ${Bam} . hg${params.hg_type} $RESOURCE/hg19.fa $RESOURCE/hg38.fa $sample_id $SRC $RESOURCE/${params.ref_name}.bed
    """
}

process Calculate_feature {
    cache "deep"; tag "$sample_id"
    // publishDir "${params.OUTDIR}/Feature", mode: "copy"
    storeDir "${params.OUTDIR}/Feature"

    input:
        tuple sample_id, file(Input_dir) from Calculate_feature_ch
        file SRC
    output:
        tuple sample_id, file("${sample_id}_features")
    script:
    """
    Rscript $SRC/features_from_GW_data.R --input $Input_dir --output .
    """
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process Read_count {
    cache "deep"; tag "$SampleID"
    storeDir "$params.OUTDIR/Read_count"

    input:
        set SampleID, file(bam), file(bai) from IchorCNA_ch
    output:
        set SampleID, file("${SampleID}.wig") into Run_IchorCNA_ch

    script:
    """
    readCounter --window ${params.window_size} --quality 20 \
                --chromosome "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY" $bam > "${SampleID}.wig"
    """
}


process Run_ichorCNA {
    cache "deep"; tag "$SampleID"
    storeDir "$params.OUTDIR/Final_result"

    input:
        set SampleID, file(wig) from Run_IchorCNA_ch
        file SRC
    output:
        file("${SampleID}{,.cna.seg,.correctedDepth.txt,.params.txt,.RData,.seg,.seg.txt}") into Summary_IchorCNA_ch
    script:
    """
    centromere_path=$SRC/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt
    normalPanel_path=$SRC/inst/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds

    Rscript $SRC/scripts/runIchorCNA.R --id $SampleID \
            --WIG $wig --ploidy "c(2)" --normal "c(0.95, 0.99, 0.995, 0.999)" --maxCN 3 \
            --gcWig $SRC/inst/extdata/gc_hg${params.hg_type}_${params.bin_size}kb.wig \
            --mapWig $SRC/inst/extdata/map_hg${params.hg_type}_${params.bin_size}kb.wig \
            --centromere \$centromere_path \
            --normalPanel \$normalPanel_path \
            --includeHOMD False --chrs "c(1:22)" --chrTrain "c(1:22)" \
            --estimateNormal True --estimatePloidy True --estimateScPrevalence False \
            --scStates "c()" --txnE 0.9999 --txnStrength 10000 --outDir ./ --seginfofilepath $SRC/inst/extdata/seqinfo_hg${params.hg_type}_ucsc.rds

    
    """
}

process Summary_tumor_fraction {
    cache "deep"; tag "$SampleID"
    storeDir "$params.OUTDIR/IchorCNA_tumor_fraction"

    input:
        file "*" from Summary_IchorCNA_ch.collect()
        file SRC
    output:
        file("IchorCNA_tumor_fraction.csv")
    script:
    """
        python3 $SRC/get_tumor_fraction.py . .
        
    """
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////