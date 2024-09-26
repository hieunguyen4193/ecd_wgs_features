# Quality control metric for WGS samples

## Setup 
Update your `java`: https://www.oracle.com/java/technologies/downloads/?er=221886#jdk23-mac or https://www.oracle.com/java/technologies/downloads/?er=221886#jdk23-linux
Download `picard.jar` from https://github.com/broadinstitute/picard/releases/tag/3.2.0.
Install `R` to enable plotting in `picard` functions.
Install `gnuplot` to enable `plot-bamstats`

See a sample output of the function `generate_QC_metrics.sh` in the folder `QC_demo`.
 
## Introduction 
Sử dụng `samtools` và `picard` để tạo ra một số QC metrics cho một file low-depth WGS `BAM`. Để xem QC cho file FASTQ thì xem ở phần `FASTQC`, các metric trong phần này áp dụng cho các reads sau khi được align lên reference genome hg19. 

## QC metrics description

| Tool     | QC name  | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     | Filename                          |
|----------|----------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------|
| samtools | flagstat | Chú ý các trường thông tin sau: - total: Tổng số reads, bao gồm cả QC passed và QC failed.  - mapped: tổng số reads map được lên reference genome.  - mapped %: tỉ lệ % cho mapped.  - read1, read2: số lượng read1 và read2. - properly paired: số lượng reads có đủ pair.  - properly paired %: tỉ lệ %.                                                                                                                                                                                                                                      | *.sorted.flagstat.txt             |
| samtools | stat     | Xem file *.html có interactive plot thể hiện các QC metrics.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    | folder samtools_stat_plots/*.html |
| picard   |          | 1. `Picard CollectWgsMetrics`  : https://gatk.broadinstitute.org/hc/en-us/articles/360037269351-CollectWgsMetrics-Picard 2. `Picard CollectGcBiasMetrics` : https://gatk.broadinstitute.org/hc/en-us/articles/360036801531-CollectGcBiasMetrics-Picard 3. `Picard CollectInsertSizeMetrics` : https://gatk.broadinstitute.org/hc/en-us/articles/360037055772-CollectInsertSizeMetrics-Picard 4. `Picard CollectAlignmentSummaryMetrics`  : https://gatk.broadinstitute.org/hc/en-us/articles/360040507751-CollectAlignmentSummaryMetrics-Picard |                                   |

## References
1. `Samtools flagstat`: http://www.htslib.org/doc/samtools-flagstat.html 
2. `Samtools stats`: http://www.htslib.org/doc/samtools-stats.html
3. `Picard CollectWgsMetrics`: https://gatk.broadinstitute.org/hc/en-us/articles/360037269351-CollectWgsMetrics-Picard
4. `Picard CollectGcBiasMetrics`: https://gatk.broadinstitute.org/hc/en-us/articles/360036801531-CollectGcBiasMetrics-Picard
5. `Picard CollectInsertSizeMetrics`: https://gatk.broadinstitute.org/hc/en-us/articles/360037055772-CollectInsertSizeMetrics-Picard
6. `Picard CollectAlignmentSummaryMetrics`: https://gatk.broadinstitute.org/hc/en-us/articles/360040507751-CollectAlignmentSummaryMetrics-Picard