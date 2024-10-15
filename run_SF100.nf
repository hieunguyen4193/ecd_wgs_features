// ----- ----- ----- CHANNEL ----- ----- -----
params.SampleSheet=""
Channel
    .fromPath( params.SampleSheet )
    .splitCsv( header:true )
    .map { row -> tuple(row.SampleID,  file(row.path_to_bam_file))}  
    .view()
    .set { Input_ch }


process run_SF100 {  
    cache "deep";
    publishDir "$params.output/SF100", mode: 'copy'
    // errorStrategy 'retry'
    // maxRetries 1
    maxForks 20

    input:
        tuple sample_id, file(bamfile) from Input_ch
    output:
        file("*") into output_ch
    script:
    """
    samtools view -f 2 -q 10 $bamfile | \
    awk -v sample=$sample_id '
    BEGIN {
        OFS="\t";
        print "Sample", "Total", "Short", "SF100", "SF100_status", "cpm_chrY", "zscoreMAPQ10chrY"
    }
    {
        if (\$9 > 0) {
            total++;
            if (\$3 == "chrY") {
                total_chrY++;
                if (˙$9 < 100) short_chrY++;
            };
        }
    }
    END {
        SF100 = total_chrY > 0 ? short_chrY / total_chrY : 0;
        SF100_status = SF100 < 0.15 ? "pass" : "fail";
       
        cpm_chrY = total_chrY * 1e6 / total;
        zscore_chrY = ( cpm_chrY - 119.3757 ) / 13.19449;
       
        print sample, total_chrY, short_chrY, SF100, SF100_status, cpm_chrY, zscore_chrY
    }' > ${sample_id}_SF100_zscore_10.tsv &
    """
}

// nextflow run run_SF100.nf \
// --SampleSheet /datassd/DATA_HIEUNGUYEN/2024/src/SF100/ecd_wgs_features/SampleSheet_to_run.SF100.csv \
// --output /datassd/DATA_HIEUNGUYEN/2024/outdir \
// -resume

