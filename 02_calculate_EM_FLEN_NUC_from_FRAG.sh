while getopts "i:o:n:" opt; do
  case ${opt} in
    i )
      inputfrag=$OPTARG
      ;;
    o )
      outputdir=$OPTARG
      ;;
    \? )
      echo "Usage: cmd [-i] input .frag file [-o] outputdir"
      exit 1
      ;;
  esac
done

sampleid=$(echo ${inputfrag} | xargs -n 1 basename)
sampleid=${sampleid%.frag*}

mkdir -p ${outputdir};

# bash 02_calculate_EM_FLEN_NUC_from_FRAG.sh -i /Users/hieunguyen/src/data/bam_files/WGShg19.frag.tsv  -o ./output/ 
##### check if the FLEN column is already in the file. 
##### our in-house data has FLEN pre-calculated, 
##### for external data in frag.tsv format, FLEN has not been calculated yet. 


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