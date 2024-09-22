#####----------------------------------------------------------------------#####
##### INTRODUCTION
#####----------------------------------------------------------------------#####
# This script pre-process an input BAM file to a 
# fragment-wise data features, which can be use to calculate
# several fragmentomics features. 
export PATH=/Users/hieunguyen/samtools/bin:$PATH

# bash 01_BAM_to_FRAGS.sh -i /Users/hieunguyen/src/data/bam_files/WGShg19.bam  -o ./output/

#####----------------------------------------------------------------------#####
##### input args
#####----------------------------------------------------------------------#####
while getopts "i:o:n:" opt; do
  case ${opt} in
    i )
      inputbam=$OPTARG
      ;;
    o )
      outputdir=$OPTARG
      ;;
    n )
      samtools_num_threads=$OPTARG
      ;;
    \? )
      echo "Usage: cmd [-i] inputbam [-o] outputdir [-n] samtools_num_threads"
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

mkdir -p ${outputdir}
sampleid=$(echo ${inputbam} | xargs -n 1 basename)
sampleid=${sampleid%.bam*}

#####----------------------------------------------------------------------#####
##### pre-processing steps for an input bam file:
#####----------------------------------------------------------------------#####
# remove unpaired and unmapped reads in BAM files.
if [ ! -f "${outputdir}/${sampleid}.frag.tsv" ]; then
  echo -e "Pre-processing BAM files and convert them to fragmentomics format ..."

  if [ ! -f "${outputdir}/${sampleid}.prep.tsv" ]; then
  echo -e "remove unpaired and unmapped reads in BAM files, generate prep.tsv file";
  samtools view -f 3 ${inputbam} -@ ${samtools_num_threads} -b -o ${outputdir}/${sampleid}.prep.bam
  samtools index -@ ${samtools_num_threads} ${outputdir}/${sampleid}.prep.bam
  samtools view ${outputdir}/${sampleid}.prep.bam | cut -f1,3,4,6,9 > ${outputdir}/${sampleid}.prep.tsv
  fi

  ##### get true fragment end
  echo -e "modify fragment end, fragment end = foward read start + fragment length"

  ##### explain new columns in the data:
  # $3: start of read 1 (fragment start)
  # $6: start of read 2 (fragment end)
  # $7 = $3 + 1: fragment start + 1, to use in nucleosome footprint feature
  # $8 = $6 + 1: fragment end + 1, to use in nucleosome footprint feature

  ##### explain the sign of TLEN field $9
  # The TLEN field is positive for the leftmost segment of the template, 
  # negative for the rightmost, and the sign for any middle segment is undefined. 
  # If segments cover the same coordinates then the choice of which is leftmost 
  # and rightmost is arbitrary, but the two ends must still have differing signs
  # https://samtools.github.io/hts-specs/SAMv1.pdf
  # this remove half of the reads, keep only fragment-wise information. 
  # 1 line = 1 fragment, chrom - start - end and a read ID, use this read ID to sort and
  # match the information later when calculating End motif and nucleosome footprint. 
  awk -v OFS='\t' '{if ($5 > 0){$6=$3+$5; $7=$1"_"$2"_"$3; print $0}}' \
    ${outputdir}/${sampleid}.prep.tsv \
    | sort -k9,9 \
    | awk '{ print $2 "\t" $3 "\t" $6 "\t" $5 "\t" $7}' > ${outputdir}/${sampleid}.frag.tsv
  rm -rf ${outputdir}/${sampleid}.prep.tsv
fi

#####----------------------------------------------------------------------#####
##### 4bp END MOTIF
#####----------------------------------------------------------------------#####
if [ ! -f "${outputdir}/${sampleid}.finished_4bpEM.txt" ]; then
    echo -e "getting 4bp end motif"
    cat ${outputdir}/${sampleid}.frag.tsv | \
      awk '{start=$2 - 1; end= $2 - 1 + 4; name= $9; strand = "+"; print $2 "\t" start "\t" end "\t" name "\t" "1" "\t" strand}' \
    > ${outputdir}/${sampleid}.forward_endcoord4bp.bed;

    cat ${outputdir}/${sampleid}.frag.tsv | \
      awk '{start=$3 - 1 - 4; end= $3 - 1; name= $9; strand = "-"; print $2 "\t" start "\t" end "\t" name "\t" "1" "\t" strand}' \
    > ${outputdir}/${sampleid}.reverse_endcoord4bp.bed;

    bedtools getfasta -s -name -tab -fi ${path_to_fa} -bed ${outputdir}/${sampleid}.forward_endcoord4bp.bed |  awk -v OFS='\t' '{split($0, a, "::"); $1=a[1]; print $0}'  > ${outputdir}/${sampleid}.forward_endmotif4bp.txt
    bedtools getfasta -s -name -tab -fi ${path_to_fa} -bed ${outputdir}/${sampleid}.reverse_endcoord4bp.bed |  awk -v OFS='\t' '{split($0, a, "::"); $1=a[1]; print $0}'  > ${outputdir}/${sampleid}.reverse_endmotif4bp.txt

    rm -rf ${outputdir}/${sampleid}.forward_endcoord4bp.bed
    rm -rf ${outputdir}/${sampleid}.reverse_endcoord4bp.bed

    sort -k1,1 ${outputdir}/${sampleid}.forward_endmotif4bp.txt > ${outputdir}/${sampleid}.forward_endmotif4bp.sorted.txt
    sort -k1,1 ${outputdir}/${sampleid}.reverse_endmotif4bp.txt > ${outputdir}/${sampleid}.reverse_endmotif4bp.sorted.txt

    touch ${outputdir}/${sampleid}.finished_4bpEM.txt
fi

count_4bpEM_forward=$(cat ${outputdir}/${sampleid}.forward_endmotif4bp.sorted.txt | wc -l)
count_4bpEM_reverse=$(cat ${outputdir}/${sampleid}.reverse_endmotif4bp.sorted.txt | wc -l)

#####----------------------------------------------------------------------#####
##### NUCLEOSOME FOOTPRINT
#####----------------------------------------------------------------------#####
if [ ! -f "${outputdir}/${sampleid}.finished_Nucleosome.txt" ]; then
  echo -e "generating nucleosome features ..."
  cat ${outputdir}/${sampleid}.frag.tsv | cut -f1,2,5 \
    | awk -v OFS='\t' '{$4=$2 + 1; print $1 "\t" $2 "\t" $4 "\t" $3}' \
    > ${outputdir}/${sampleid}.forward_Nucleosome.bed
  cat ${outputdir}/${sampleid}.frag.tsv | cut -f1,3,5 \
    | awk -v OFS='\t' '{$4=$2 + 1; print $1 "\t" $2 "\t" $4 "\t" $3}' \
    > ${outputdir}/${sampleid}.reverse_Nucleosome.bed

  # Sort your generated BED files
  sort -k 1V,1 -k 2n,2 ${outputdir}/${sampleid}.forward_Nucleosome.bed -o ${outputdir}/${sampleid}.sortedNuc.forward_Nucleosome.bed
  sort -k 1V,1 -k 2n,2 ${outputdir}/${sampleid}.reverse_Nucleosome.bed -o ${outputdir}/${sampleid}.sortedNuc.reverse_Nucleosome.bed

  # Sort the reference BED file
  sort -k 1V,1 -k 2n,2 ${nucleosome_ref} -o ${nucleosome_ref%.bed*}.sorted.bed

  ##### sort with -k 1V,1 to get the correct order of chromosome, add -t first to get first nucleosome only, match row numbers.
  bedtools closest -a ${outputdir}/${sampleid}.sortedNuc.forward_Nucleosome.bed -b ${nucleosome_ref%.bed*}.sorted.bed -t first | awk -v OFS='\t' '{$13=$11 - $2;print $0}' > ${outputdir}/${sampleid}.forward_Nucleosome.dist.bed
  bedtools closest -a ${outputdir}/${sampleid}.sortedNuc.reverse_Nucleosome.bed -b ${nucleosome_ref%.bed*}.sorted.bed -t first | awk -v OFS='\t' '{$13=$11 - $2;print $0}' > ${outputdir}/${sampleid}.reverse_Nucleosome.dist.bed

# | awk -v OFS='\t' '{$14=$2-$6;print $0}'
  echo -e "sorting forward nucleosome file"
  sort -k4,4 ${outputdir}/${sampleid}.forward_Nucleosome.dist.bed > ${outputdir}/${sampleid}.forward_Nucleosome.dist.sorted.bed
  echo -e "sorting reverse nucleosome file"
  sort -k4,4 ${outputdir}/${sampleid}.reverse_Nucleosome.dist.bed > ${outputdir}/${sampleid}.reverse_Nucleosome.dist.sorted.bed

  touch ${outputdir}/${sampleid}.finished_Nucleosome.txt
fi

#####----------------------------------------------------------------------#####
##### Merge all features into one single tsv output file. 
#####----------------------------------------------------------------------#####
if [ ! -f "${outputdir}/${sampleid}.final_output.tsv" ]; then
    echo -e "merge output"
    count_nuc_forward=$(cat ${outputdir}/${sampleid}.forward_Nucleosome.dist.sorted.bed | wc -l)
    count_nuc_reverse=$(cat ${outputdir}/${sampleid}.reverse_Nucleosome.dist.sorted.bed | wc -l)

    ##### column $10: distance of forward read to the nearest nucleosome
    cat ${outputdir}/${sampleid}.forward_Nucleosome.dist.sorted.bed | cut -f13 > ${outputdir}/forward_nuc.tmp.txt
    paste ${outputdir}/${sampleid}.frag.tsv ${outputdir}/forward_nuc.tmp.txt  > ${outputdir}/${sampleid}.modified1.tsv

    ##### column $11: distance of reverse read to the nearest nucleosome
    cat ${outputdir}/${sampleid}.reverse_Nucleosome.dist.sorted.bed | cut -f13 > ${outputdir}/reverse_nuc.tmp.txt
    paste ${outputdir}/${sampleid}.modified1.tsv ${outputdir}/reverse_nuc.tmp.txt  > ${outputdir}/${sampleid}.modified2.tsv

    ##### column $12: 4bp end motif from forward reads
    cat ${outputdir}/${sampleid}.forward_endmotif4bp.sorted.txt | cut -f2 > ${outputdir}/forward_4bpEM.tmp.txt
    paste ${outputdir}/${sampleid}.modified2.tsv ${outputdir}/forward_4bpEM.tmp.txt  > ${outputdir}/${sampleid}.modified3.tsv

    # ##### column $13: 4bp end motif from reverse reads
    cat ${outputdir}/${sampleid}.reverse_endmotif4bp.sorted.txt | cut -f2 > ${outputdir}/reverse_4bpEM.tmp.txt
    paste ${outputdir}/${sampleid}.modified3.tsv ${outputdir}/reverse_4bpEM.tmp.txt  > ${outputdir}/${sampleid}.modified4.tsv

    mv ${outputdir}/${sampleid}.modified4.tsv ${outputdir}/${sampleid}.final_output.tsv
    rm -rf ${outputdir}/${sampleid}.modified{1,2,3,4}.tsv
fi

