#####----------------------------------------------------------------------------#####
##### INPUT ARGS
#####----------------------------------------------------------------------------#####
path_to_bam_file=$1;
path_to_output=$2;
path_to_fa=$3;
sampleid=$4;
path_to_SRC=$5
path_to_REF=$6

run_methyl_extract="yes";

path_to_output=${path_to_output}/${sampleid};
echo -e "Main output directory: " $path_to_output "\n"
mkdir -p ${path_to_output};
#####----------------------------------------------------------------------------#####
##### Remove unmapped and un-pair reads from bam file
#####----------------------------------------------------------------------------#####
echo -e "Remove unmapped and unpaired reads from bam file ..."
samtools view -h -f 3 ${path_to_bam_file} | samtools sort -n -@ 5 -o tmp.bam;
samtools view -h tmp.bam | awk -f preprocessing_script.awk - > tmp.sam;

echo -e "Sort and index bam file ..."
samtools sort -@ 5 -O BAM -o ${path_to_output}/${sampleid}.sorted.bam tmp.sam;
rm -rf tmp.sam;
rm -rf tmp.bam
samtools index -@ 5 ${path_to_output}/${sampleid}.sorted.bam

echo -e "Finished! Use the new bam file for downstream tasks"

#####----------------------------------------------------------------------------#####
##### re assign path_to_bam_file
path_to_bam_file=${path_to_output}/${sampleid}.sorted.bam;

echo -e "Working on sample " $sampleid "\n";
echo -e "Sample saved at " ${path_to_bam_file} "\n";


echo -e "Using the following hg file: \n";
echo $path_to_fa;

path_to_short_bam=${path_to_output}/short_bam;
mkdir -p ${path_to_short_bam};
path_to_long_bam=${path_to_output}/long_bam;
mkdir -p ${path_to_long_bam};
path_to_flen=${path_to_output}/flen;
mkdir -p ${path_to_flen};
path_to_em=${path_to_output}/EM;
mkdir -p ${path_to_em};
path_to_nucleosome=${path_to_output}/nucleosome;
mkdir -p ${path_to_nucleosome};

path_to_short_bam_modify_cutoff=${path_to_output}/short_bam_modify_cutoff;
mkdir -p ${path_to_short_bam_modify_cutoff};
path_to_long_bam_modify_cutoff=${path_to_output}/long_bam_modify_cutoff;
mkdir -p ${path_to_long_bam_modify_cutoff};





##### Nucleosome
# samtools view -@ 8 ${path_to_bam_file} | awk '{if ($9 > 0){print $0}}'  | cut -f3,4 > $path_to_nucleosome/${sampleid}.bed
# awk -v OFS='\t' '{$3=$2+1; print $0}' $path_to_nucleosome/${sampleid}.bed > $path_to_nucleosome/${sampleid}.full.bed
samtools view -@ 8 ${path_to_bam_file} | awk '{if ($9 > 0){chrom=$3;start=$4} else {chrom=$3;start=$8 - $9}; print chrom "\t" start}' > $path_to_nucleosome/${sampleid}.bed
awk -v OFS='\t' '{$3=$2+1; print $0}' $path_to_nucleosome/${sampleid}.bed > $path_to_nucleosome/${sampleid}.full.bed

python3 $path_to_SRC/convert_full_bed_nucleosome.py $path_to_nucleosome/${sampleid}.full.bed $path_to_nucleosome/${sampleid}.full.bed

bedtools closest -a $path_to_nucleosome/${sampleid}.full.bed -b $path_to_REF | cut -f1,2,10 > $path_to_nucleosome/${sampleid}.rpr_map.bed
awk -v OFS='\t' '{$4=$3-$2; print $0}' $path_to_nucleosome/${sampleid}.rpr_map.bed > $path_to_nucleosome/${sampleid}.dist.bed
python3 $path_to_SRC/count_density.py $path_to_nucleosome/${sampleid}.dist.bed $path_to_nucleosome



##### Extract fragment lengths from bam file
samtools view ${path_to_bam_file} | cut -f9 > ${path_to_flen}/${sampleid}.flen.txt;




##### Split bam into short and long bam
samtools view -f 2 -h ${path_to_bam_file} | \
  awk 'substr($0,1,1)=="@" || ($9 >= 50 && $9 <= 150) || ($9 <= -50 && $9 >= -150)' | \
  samtools view -b > ${path_to_short_bam}/${sampleid}.short.bam;
samtools index ${path_to_short_bam}/${sampleid}.short.bam
samtools view -f 2 -h ${path_to_bam_file} | \
  awk 'substr($0,1,1)=="@" || ($9 > 150 && $9 <= 250) || ($9 < -150 && $9 >= -250)' | \
  samtools view -b > ${path_to_long_bam}/${sampleid}.long.bam;
samtools index ${path_to_long_bam}/${sampleid}.long.bam;


##### Split bam into short and long bam with new cutoff
samtools view -f 2 -h ${path_to_bam_file} | \
  awk 'substr($0,1,1)=="@" || ($9 >= 100 && $9 <= 150) || ($9 <= -100 && $9 >= -150)' | \
  samtools view -b > ${path_to_short_bam_modify_cutoff}/${sampleid}.short.bam;
samtools index ${path_to_short_bam_modify_cutoff}/${sampleid}.short.bam
samtools view -f 2 -h ${path_to_bam_file} | \
  awk 'substr($0,1,1)=="@" || ($9 > 151 && $9 <= 220) || ($9 < -151 && $9 >= -220)' | \
  samtools view -b > ${path_to_long_bam_modify_cutoff}/${sampleid}.long.bam;
samtools index ${path_to_long_bam_modify_cutoff}/${sampleid}.long.bam;




##### Extract end motif from bam file
samtools view -@ 4 -q 30 -f 2 ${path_to_bam_file} | awk '{if($9 > 0) {start=$4 - 1; end= $4 - 1 + 4; name= $1"_"$9; strand = "+"} else {start=$8 - $9 - 4 - 1; end= $8 - $9 -1; name= $1"_"$9; strand = "-"}; print $3 "\t" start "\t" end "\t" name "\t" "1" "\t" strand}' > ${sampleid}.endcoord4bp.bed;

bedtools getfasta -s -name -tab -fi ${path_to_fa} -bed ${sampleid}.endcoord4bp.bed  > ${path_to_em}/${sampleid}.endmotif4bp.txt
rm ${sampleid}.endcoord4bp.bed;

#####################################################################################################################

##### Extract flagment lengths for each chromosome
path_to_flen_chr=${path_to_flen}_chr
mkdir -p $path_to_flen_chr

##### Extract em for each chromosome
path_to_em_chr=${path_to_em}_chr
mkdir -p $path_to_em_chr

##### Extract em for each chromosome
path_to_nucleosome_chr=${path_to_nucleosome}_chr
mkdir -p $path_to_nucleosome_chr



samtools idxstats $path_to_bam_file | cut -f 1 | grep -v '*' | while read chr; do 

chromosome=chr${chr##*chr}; 
if [[ $chromosome =~ ^chr([0-9]+)$ ]]; then
    chrom_number=${BASH_REMATCH[1]}
    # Check if the number is between 1 and 22
    if (( chrom_number >= 1 && chrom_number <= 22 )); then
        echo $chr " ===> PASS " $chrom_number;

        # FLEN
        samtools view -b $path_to_bam_file $chr | samtools view | cut -f9 >> ${path_to_flen_chr}/${sampleid}.chr.${chr##*chr}.flen.txt

        # EM
        samtools view -b $path_to_bam_file $chr | samtools view -@ 4 -q 30 -f 2 | awk '{if($9 > 0) {start=$4 - 1; end= $4 - 1 + 4; name= $1"_"$9; strand = "+"} else {start=$8 - $9 - 4 - 1; end= $8 - $9 -1; name= $1"_"$9; strand = "-"}; print $3 "\t" start "\t" end "\t" name "\t" "1" "\t" strand}' > ${path_to_em_chr}/${sampleid}.chr.${chr##*chr}.EM.bed
        bedtools getfasta -s -name -tab -fi ${path_to_fa} -bed ${path_to_em_chr}/${sampleid}.chr.${chr##*chr}.EM.bed  > ${path_to_em_chr}/${sampleid}.chr.${chr##*chr}.EM.txt
        rm ${path_to_em_chr}/${sampleid}.chr.${chr##*chr}.EM.bed

        ##### Nucleosome
        python3 $path_to_SRC/split_nucleosome_by_each_chr.py $path_to_nucleosome/${sampleid}.dist.bed $chr $path_to_nucleosome_chr/${sampleid}.chr.${chr##*chr}.dist.bed

        python3 $path_to_SRC/count_density.py $path_to_nucleosome_chr/${sampleid}.chr.${chr##*chr}.dist.bed $path_to_nucleosome_chr
        mv $path_to_nucleosome_chr/${sampleid}.csv $path_to_nucleosome_chr/${sampleid}.chr.${chr##*chr}.csv
        rm $path_to_nucleosome_chr/${sampleid}.chr.${chr##*chr}.dist.bed
    fi
fi

done

python3 $path_to_SRC/merge_feature_nucleosome_each_chr.py $path_to_nucleosome_chr $sampleid $path_to_nucleosome_chr

rm $path_to_nucleosome/${sampleid}{,.full,.rpr_map,.dist}.bed
rm $path_to_nucleosome_chr/${sampleid}.chr.*.csv




touch ${path_to_output}/sample_${sampleid}.finished.txt;
