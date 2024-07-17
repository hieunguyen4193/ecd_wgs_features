#!/bin/bash
# inputbam=$1
# outputdir=$2
# num_threads=$3

# get samtools to work
export PATH=/Users/hieunguyen/samtools/bin:$PATH 
# samtools --version

# default input values
# inputbam="/Volumes/HNSD02/data/WGS_bam/9-ZMC014NB_S95025-S97025.sorted.bam";
# samtools view "/Volumes/HNSD02/data/WGS_bam/9-ZMC014NB_S95025-S97025.sorted.bam" -h | head -n 20000 | samtools view -b - > input_test.bam 
inputbam="./examples/input_test.bam";
outputdir="/Volumes/HNSD02/outdir/ecd_wgs_features";
path_to_fa="/Volumes/HNSD02/resource/hg19.fa";
num_threads=10;
nucleosome_ref="/Volumes/HNSD02/resource/nucleosome_footprint_bed/rpr_map_Budhraja_STM2023.bed"

# Parse command line arguments
while getopts "i:o:t:r:" opt; do
  case ${opt} in
    i )
      inputbam=$OPTARG
      ;;
    o )
      outputdir=$OPTARG
      ;;
    t )
      num_threads=$OPTARG
      ;;
    r )
      nucleosome_ref=$OPTARG
      ;;
    \? )
      echo "Usage: cmd [-i] inputbam [-o] outputdir [-t] num_threads [-r] nucleosome reference bed file"
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

if [ ! -f "${outputdir}/${filename}.modified.tsv" ]; then
    echo -e "pre-processing BAM file to remove unmapped reads"
    echo -e "filter bam file to keep chr1-22 only ..."
    samtools view -b ${inputbam} chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 > ${outputdir}/${filename}.filterChr.bam
    samtools index -@ ${num_threads} ${outputdir}/${filename}.filterChr.bam;
    
    inputbam=${outputdir}/${filename}.filterChr.bam;
    samtools flagstat ${inputbam} > ${outputdir}/${filename}.filterChr.flagstat.txt

    echo -e "excluding unmapped reads ... "
    # samtools view --excl-flags 4 -b ${inputbam} > ${outputdir}/${filename}.tmp.bam
    # echo -e "sorting bam file according to read name ... "
    # samtools sort -n ${outputdir}/${filename}.tmp.bam -@ ${num_threads} -o ${outputdir}/${filename}.tmp.sortedN.bam;
    # echo -e "filter bam file: keep only paired mate reads..."
    # samtools view -h ${outputdir}/${filename}.tmp.sortedN.bam | awk -f preprocessing_script.awk - > ${outputdir}/${filename}.tmp.sam;
    # echo -e "re-sort bam file according to coordinate ..."
    # samtools sort -@ ${num_threads} -O BAM -o ${outputdir}/${filename}.prep.bam ${outputdir}/${filename}.tmp.sam;
    # rm ${outputdir}/${filename}.tmp.sam;

    # flag 3: read paired and read mapped in proper pair  
    samtools view -f 3 ${inputbam} -@ ${num_threads} -b -o ${outputdir}/${filename}.prep.bam

    echo -e "index final output bam ..."
    samtools index -@ ${num_threads} ${outputdir}/${filename}.prep.bam
    
    samtools view ${outputdir}/${filename}.prep.bam | cut -f1,3,4,6,9 > ${outputdir}/${filename}.prep.tsv
    
    samtools flagstat ${outputdir}/${filename}.prep.bam > ${outputdir}/${filename}.prep.flagstat.txt

    ##### get true fragment end
    echo -e "modify fragment end, fragment end = foward read start + fragment length"
    ##### explain new columns in the data:
    # $3: start of read 1 (fragment start)
    # $6: start of read 2 (fragment end)
    # $7 = $3 + 1: fragment start + 1, to use in nucleosome footprint feature
    # $8 = $6 + 1: fragment end + 1, to use in nucleosome footprint feature

    ##### explain the sign of TLEN field $9
    # The TLEN field is positive for the leftmost segment of the template, negative for the rightmost, and the sign for any middle segment is undefined. 
    # If segments cover the same coordinates then the choice of which is leftmost and rightmost is arbitrary, but the two ends must still have differing signs
    # https://samtools.github.io/hts-specs/SAMv1.pdf
    awk -v OFS='\t' '{if ($5 > 0){$6=$3+$5;$7=$3+1;$8=$6+1;$9=$1"_"$2"_"$3; print $0}}' ${outputdir}/${filename}.prep.tsv | sort -k9,9 > ${outputdir}/${filename}.modified.tsv

    count_main_file=$(cat ${outputdir}/${filename}.modified.tsv | wc -l)
fi

#####----------------------------------------------------------------------#####
##### 4bp END MOTIF
#####----------------------------------------------------------------------#####
if [ ! -f "${outputdir}/${filename}.finished_4bpEM.txt" ]; then
    echo -e "getting 4bp end motif"
    cat ${outputdir}/${filename}.modified.tsv | \
      awk '{start=$3 - 1; end= $3 - 1 + 4; name= $9; strand = "+"; print $2 "\t" start "\t" end "\t" name "\t" "1" "\t" strand}' \
    > ${outputdir}/${filename}.forward_endcoord4bp.bed;

    cat ${outputdir}/${filename}.modified.tsv | \
      awk '{start=$6 - 1 - 4; end= $6 - 1; name= $9; strand = "-"; print $2 "\t" start "\t" end "\t" name "\t" "1" "\t" strand}' \
    > ${outputdir}/${filename}.reverse_endcoord4bp.bed;

    bedtools getfasta -s -name -tab -fi ${path_to_fa} -bed ${outputdir}/${filename}.forward_endcoord4bp.bed |  awk -v OFS='\t' '{split($0, a, "::"); $1=a[1]; print $0}'  > ${outputdir}/${filename}.forward_endmotif4bp.txt
    bedtools getfasta -s -name -tab -fi ${path_to_fa} -bed ${outputdir}/${filename}.reverse_endcoord4bp.bed |  awk -v OFS='\t' '{split($0, a, "::"); $1=a[1]; print $0}'  > ${outputdir}/${filename}.reverse_endmotif4bp.txt

    rm -rf ${outputdir}/${filename}.forward_endcoord4bp.bed
    rm -rf ${outputdir}/${filename}.reverse_endcoord4bp.bed

    sort -k1,1 ${outputdir}/${filename}.forward_endmotif4bp.txt > ${outputdir}/${filename}.forward_endmotif4bp.sorted.txt
    sort -k1,1 ${outputdir}/${filename}.reverse_endmotif4bp.txt > ${outputdir}/${filename}.reverse_endmotif4bp.sorted.txt

    touch ${outputdir}/${filename}.finished_4bpEM.txt
fi

count_4bpEM_forward=$(cat ${outputdir}/${filename}.forward_endmotif4bp.sorted.txt | wc -l)
count_4bpEM_reverse=$(cat ${outputdir}/${filename}.reverse_endmotif4bp.sorted.txt | wc -l)

#####----------------------------------------------------------------------#####
##### 21bp END MOTIF
#####----------------------------------------------------------------------#####
# if [ ! -f "${outputdir}/${filename}.finished_21bpEM.txt" ]; then
#     echo -e "getting 21bp end motif"
#     cat ${outputdir}/${filename}.modified.tsv | \
#       awk '{if($3 > 11) {start=$3 - 1 - 10; end= $3 - 1 + 11; name= $9; strand = "+"}; print $2 "\t" start "\t" end "\t" name "\t" "1" "\t" strand}' \
#     > ${outputdir}/${filename}.forward_endcoord21bp.bed;

#     cat ${outputdir}/${filename}.modified.tsv | \
#       awk '{if($6 > 11) {start=$6 - 1 - 11; end= $6 - 1 + 10; name= $9; strand = "-"}; print $2 "\t" start "\t" end "\t" name "\t" "1" "\t" strand}' \
#     > ${outputdir}/${filename}.reverse_endcoord21bp.bed;

#     bedtools getfasta -s -name -tab -fi ${path_to_fa} -bed ${outputdir}/${filename}.forward_endcoord21bp.bed | awk '{split($0, a, "::"); $1=a[1]; print $0}'  > ${outputdir}/${filename}.forward_endmotif21bp.txt
#     bedtools getfasta -s -name -tab -fi ${path_to_fa} -bed ${outputdir}/${filename}.reverse_endcoord21bp.bed | awk '{split($0, a, "::"); $1=a[1]; print $0}' > ${outputdir}/${filename}.reverse_endmotif21bp.txt

#     rm -rf ${outputdir}/${filename}.forward_endcoord21bp.bed
#     rm -rf ${outputdir}/${filename}.reverse_endcoord21bp.bed
    
#     sort -k1,1 ${outputdir}/${filename}.forward_endmotif21bp.txt > ${outputdir}/${filename}.forward_endmotif21bp.sorted.txt
#     sort -k1,1 ${outputdir}/${filename}.reverse_endmotif21bp.txt > ${outputdir}/${filename}.reverse_endmotif21bp.sorted.txt
    
#     touch ${outputdir}/${filename}.finished_21bpEM.txt
# fi

# count_21bpEM_forward=$(cat ${outputdir}/${filename}.forward_endmotif21bp.sorted.txt | wc -l)
# count_21bpEM_reverse=$(cat ${outputdir}/${filename}.reverse_endmotif21bp.sorted.txt | wc -l)

#####----------------------------------------------------------------------#####
##### NUCLEOSOME FOOTPRINT
#####----------------------------------------------------------------------#####

if [ ! -f "${outputdir}/${filename}.finished_Nucleosome.txt" ]; then
  echo -e "generating nucleosome features ..."
  cat ${outputdir}/${filename}.modified.tsv | cut -f2,3,7,9 > ${outputdir}/${filename}.forward_Nucleosome.bed
  cat ${outputdir}/${filename}.modified.tsv | cut -f2,6,8,9 > ${outputdir}/${filename}.reverse_Nucleosome.bed
  # Sort your generated BED files
  sort -k1,1 -k2,2n ${outputdir}/${filename}.forward_Nucleosome.bed -o ${outputdir}/${filename}.sortedNuc.forward_Nucleosome.bed
  sort -k1,1 -k2,2n ${outputdir}/${filename}.reverse_Nucleosome.bed -o ${outputdir}/${filename}.sortedNuc.reverse_Nucleosome.bed

  # Sort the reference BED file
  sort -k1,1 -k2,2n ${nucleosome_ref} -o ${nucleosome_ref%.bed*}.sorted.bed

  bedtools closest -d -D ref -a ${outputdir}/${filename}.sortedNuc.forward_Nucleosome.bed -b ${nucleosome_ref%.bed*}.sorted.bed -t first  > ${outputdir}/${filename}.forward_Nucleosome.dist.bed
  bedtools closest -d -D ref -a ${outputdir}/${filename}.sortedNuc.reverse_Nucleosome.bed -b ${nucleosome_ref%.bed*}.sorted.bed -t first > ${outputdir}/${filename}.reverse_Nucleosome.dist.bed
# | awk -v OFS='\t' '{$14=$2-$6;print $0}'
  echo -e "sorting forward nucleosome file"
  sort -k4,4 ${outputdir}/${filename}.forward_Nucleosome.dist.bed > ${outputdir}/${filename}.forward_Nucleosome.dist.sorted.bed
  echo -e "sorting reverse nucleosome file"
  sort -k4,4 ${outputdir}/${filename}.reverse_Nucleosome.dist.bed > ${outputdir}/${filename}.reverse_Nucleosome.dist.sorted.bed

  touch ${outputdir}/${filename}.finished_Nucleosome.txt
fi

count_nuc_forward=$(cat ${outputdir}/${filename}.forward_Nucleosome.dist.sorted.bed | wc -l)
count_nuc_reverse=$(cat ${outputdir}/${filename}.reverse_Nucleosome.dist.sorted.bed | wc -l)

##### column $10: distance of forward read to the nearest nucleosome
cat ${outputdir}/${filename}.forward_Nucleosome.dist.sorted.bed | cut -f13 > ${outputdir}/forward_nuc.tmp.txt
paste ${outputdir}/${filename}.modified.tsv ${outputdir}/forward_nuc.tmp.txt  > ${outputdir}/${filename}.modified1.tsv

##### column $11: distance of forward read to the nearest nucleosome
cat ${outputdir}/${filename}.reverse_Nucleosome.dist.sorted.bed | cut -f13 > ${outputdir}/reverse_nuc.tmp.txt
paste ${outputdir}/${filename}.modified1.tsv ${outputdir}/reverse_nuc.tmp.txt  > ${outputdir}/${filename}.modified2.tsv

##### column $12: 4bp end motif from forward reads
cat ${outputdir}/${filename}.forward_endmotif4bp.sorted.txt | cut -d" " -f2 > ${outputdir}/forward_4bpEM.tmp.txt
paste ${outputdir}/${filename}.modified2.tsv ${outputdir}/forward_4bpEM.tmp.txt  > ${outputdir}/${filename}.modified3.tsv

# ##### column $13: 4bp end motif from reverse reads
cat ${outputdir}/${filename}.reverse_endmotif4bp.sorted.txt | cut -d" " -f2 > ${outputdir}/reverse_4bpEM.tmp.txt
paste ${outputdir}/${filename}.modified3.tsv ${outputdir}/reverse_4bpEM.tmp.txt  > ${outputdir}/${filename}.modified4.tsv

mv ${outputdir}/${filename}.modified4.tsv ${outputdir}/${filename}.final_output.tsv
rm -rf ${outputdir}/${filename}.modified{1,2,3,4}.tsv