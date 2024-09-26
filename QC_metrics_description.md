# Quality control metric for WGS samples

## Setup 
Update your `java`: https://www.oracle.com/java/technologies/downloads/?er=221886#jdk23-mac or https://www.oracle.com/java/technologies/downloads/?er=221886#jdk23-linux
Download `picard.jar` from https://github.com/broadinstitute/picard/releases/tag/3.2.0.
Install `R` to enable plotting in `picard` functions.
Install `gnuplot` to enable `plot-bamstats`

## Introduction 
Sử dụng `samtools` và `picard` để tạo ra một số QC metrics cho một file low-depth WGS `BAM`. Để xem QC cho file FASTQ thì xem ở phần `FASTQC`, các metric trong phần này áp dụng cho các reads sau khi được align lên reference genome hg19. 

## QC metrics description
Trong folder QC sẽ có các output như sau:
- 

## References
1. `Samtools flagstat`: http://www.htslib.org/doc/samtools-flagstat.html 
2. `Samtools stats`: http://www.htslib.org/doc/samtools-stats.html
3. `Picard CollectWgsMetrics`: https://gatk.broadinstitute.org/hc/en-us/articles/360037269351-CollectWgsMetrics-Picard
4. `Picard CollectGcBiasMetrics`: https://gatk.broadinstitute.org/hc/en-us/articles/360036801531-CollectGcBiasMetrics-Picard
5. `Picard CollectInsertSizeMetrics`: https://gatk.broadinstitute.org/hc/en-us/articles/360037055772-CollectInsertSizeMetrics-Picard