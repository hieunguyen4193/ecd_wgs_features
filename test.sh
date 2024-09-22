export PATH=/Users/hieunguyen/samtools/bin:$PATH
echo -e "Remove unmapped and unpaired reads from bam file ..."
path_to_bam_file="/Users/hieunguyen/data/bam_files/WGShg19.bam";
sampleid="WGShg19"
path_to_output="/Users/hieunguyen/src/ecd_wgs_features"

samtools view -h -f 3 ${path_to_bam_file} | samtools sort -n -@ 5 -o tmp.bam;
samtools view -h tmp.bam | awk -f preprocessing_script.awk - > tmp.sam;

echo -e "Sort and index bam file ..."
samtools sort -@ 5 -O BAM -o ${path_to_output}/${sampleid}.sorted.bam tmp.sam;
rm -rf tmp.sam;
rm -rf tmp.bam
samtools index -@ 5 ${path_to_output}/${sampleid}.sorted.bam

echo -e "Finished! Use the new bam file for downstream tasks"