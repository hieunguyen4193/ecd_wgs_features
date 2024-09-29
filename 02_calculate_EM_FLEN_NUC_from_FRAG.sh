while getopts "i:o:f:r:n:" opt; do
  case ${opt} in
    i )
      inputfrag=$OPTARG
      ;;
    o )
      outputdir=$OPTARG
      ;;
    f )
      path_to_fa=$OPTARG
      ;;
    r )
      nucleosome_ref=$OPTARG
      ;;  
    \? )
      echo "Usage: cmd [-i] input .frag file [-o] outputdir [-f] path_to_fa [-r] nucleosome_ref"
      exit 1
      ;;
  esac
done

sampleid=$(echo ${inputfrag} | xargs -n 1 basename)
sampleid=${sampleid%.frag*}

mkdir -p ${outputdir};
export PATH=/home/hieunguyen/bedtools2/bin/:$PATH

# bash 02_calculate_EM_FLEN_NUC_from_FRAG.sh -i ./output/WGShg19.frag.tsv  -o ./output/ -f /Users/hieunguyen/data/resources/hg19.fa -r /Users/hieunguyen/data/resources/rpr_map_EXP0779.sorted.bed
# bash 02_calculate_EM_FLEN_NUC_from_FRAG.sh -i ./output/WGShg19.frag.tsv  -o ./output/ -f /media/hieunguyen/HNSD01/resources/hg19.fa -r /media/hieunguyen/HNSD01/resources/rpr_map_Budhraja_STM2023.bed
# bash 02_calculate_EM_FLEN_NUC_from_FRAG.sh -i ./output_debug/1-ZLAAO90NB_S7509-S7709.sorted.frag.tsv  -o ./output_debug/ -f /Users/hieunguyen/data/resources/hg19.fa -r /Users/hieunguyen/data/resources/rpr_map_EXP0779.sorted.bed

##### check if the FLEN column is already in the file. 
##### our in-house data has FLEN pre-calculated, 
##### for external data in frag.tsv format, FLEN has not been calculated yet. 
echo -e "Working on the file " ${inputfrag}
count_col=$(awk -F'\t' '{print NF; exit}' ${inputfrag})

if [ $count_col -lt 3 ]; then
    echo -e "column FLEN does not exist in the frag.tsv file, generate new column FLEN ...";
    cat ${inputfrag} | awk -v OFS='\t' '{$4 = $3 - $2}; print $0' > ${inputfrag%.frag*}.addedFLEN.frag.tsv
    inputfrag=${inputfrag%.frag*}.addedFLEN.frag.tsv
fi

#####----------------------------------------------------------------------#####
##### 4bp END MOTIF
#####----------------------------------------------------------------------#####
if [ ! -f "${outputdir}/${sampleid}.finished_4bpEM.txt" ]; then
    echo -e "getting 4bp end motif"
    cat ${inputfrag} | \
      awk '{start=$2 - 1; end= $2 - 1 + 4; name= $5; strand = "+"; print $1 "\t" start "\t" end "\t" name "\t" "1" "\t" strand}' \
    > ${outputdir}/${sampleid}.forward_endcoord4bp.bed;

    cat ${inputfrag} | \
      awk '{start=$3 - 1 - 4; end= $3 - 1; name= $5; strand = "-"; print $1 "\t" start "\t" end "\t" name "\t" "1" "\t" strand}' \
    > ${outputdir}/${sampleid}.reverse_endcoord4bp.bed;

    bedtools getfasta -s -name -tab -fi ${path_to_fa} -bed ${outputdir}/${sampleid}.forward_endcoord4bp.bed | \
      awk -v OFS='\t' '{split($0, a, "::"); $1=a[1]; print $0}'  > ${outputdir}/${sampleid}.forward_endmotif4bp.txt
    bedtools getfasta -s -name -tab -fi ${path_to_fa} -bed ${outputdir}/${sampleid}.reverse_endcoord4bp.bed | \
    awk -v OFS='\t' '{split($0, a, "::"); $1=a[1]; print $0}'  > ${outputdir}/${sampleid}.reverse_endmotif4bp.txt

    # rm -rf ${outputdir}/${sampleid}.forward_endcoord4bp.bed
    # rm -rf ${outputdir}/${sampleid}.reverse_endcoord4bp.bed

    sort -k1,1 ${outputdir}/${sampleid}.forward_endmotif4bp.txt > ${outputdir}/${sampleid}.forward_endmotif4bp.sorted.txt
    sort -k1,1 ${outputdir}/${sampleid}.reverse_endmotif4bp.txt > ${outputdir}/${sampleid}.reverse_endmotif4bp.sorted.txt

    ##### Note: temporaily use this scripts to get the End MOTIF features. This is not the BEST way to get the nucleosome features, but we will use it for now.
    cat ${inputfrag} | \
      awk '{if($4 > 0){print $0}}' |\
      awk '{if($6 >= 30){print $0}}' |\
      awk '{start=$2 - 1; end= $2 - 1 + 4; name= $5; strand = "+"; print $1 "\t" start "\t" end "\t" name "\t" "1" "\t" strand}' \
    > ${outputdir}/${sampleid}.full_endcoord4bp.bed;

    cat ${inputfrag} | \
      awk '{if($4 < 0){print $0}}' |\
      awk '{if($6 >= 30){print $0}}' |\
      awk '{start=$3 - 1 - 4; end= $3 - 1; name= $5; strand = "-"; print $1 "\t" start "\t" end "\t" name "\t" "1" "\t" strand}' \
    >> ${outputdir}/${sampleid}.full_endcoord4bp.bed;

    bedtools getfasta -s -name -tab -fi ${path_to_fa} -bed ${outputdir}/${sampleid}.full_endcoord4bp.bed > ${outputdir}/${sampleid}.full_endmotif4bp.sorted.txt

    touch ${outputdir}/${sampleid}.finished_4bpEM.txt


fi

count_4bpEM_forward=$(cat ${outputdir}/${sampleid}.forward_endmotif4bp.sorted.txt | wc -l)
count_4bpEM_reverse=$(cat ${outputdir}/${sampleid}.reverse_endmotif4bp.sorted.txt | wc -l)

#####----------------------------------------------------------------------#####
##### NUCLEOSOME FOOTPRINT
#####----------------------------------------------------------------------#####
if [ ! -f "${outputdir}/${sampleid}.finished_Nucleosome.txt" ]; then
  echo -e "generating nucleosome features ..."
  cat ${inputfrag} | cut -f1,2,4,5 | \
    awk -v OFS='\t' '{$5=$2 + 1; print $1 "\t" $2 "\t" $5 "\t" $4 "\t" $3}' \
    > ${outputdir}/${sampleid}.forward_Nucleosome.bed
  cat ${inputfrag} | cut -f1,3,4,5 | \
    awk -v OFS='\t' '{$5=$2 + 1; print $1 "\t" $2 "\t" $5 "\t" $4 "\t" $3}' \
    > ${outputdir}/${sampleid}.reverse_Nucleosome.bed

  # Sort your generated BED files
  sort -k 1V,1 -k 2n,2 ${outputdir}/${sampleid}.forward_Nucleosome.bed -o ${outputdir}/${sampleid}.sortedNuc.forward_Nucleosome.bed
  sort -k 1V,1 -k 2n,2 ${outputdir}/${sampleid}.reverse_Nucleosome.bed -o ${outputdir}/${sampleid}.sortedNuc.reverse_Nucleosome.bed

  ##### sort with -k 1V,1 to get the correct order of chromosome, add -t first to get first nucleosome only, match row numbers.
  bedtools closest -a ${outputdir}/${sampleid}.sortedNuc.forward_Nucleosome.bed -b ${nucleosome_ref} -t first | awk -v OFS='\t' '{$13=$11 - $2;print $0}' > ${outputdir}/${sampleid}.forward_Nucleosome.dist.bed
  bedtools closest -a ${outputdir}/${sampleid}.sortedNuc.reverse_Nucleosome.bed -b ${nucleosome_ref} -t first | awk -v OFS='\t' '{$13=$11 - $2;print $0}' > ${outputdir}/${sampleid}.reverse_Nucleosome.dist.bed

  ##### Note: temporaily use this scripts to get the nucleosome features. This is not the BEST way to get the nucleosome features, but we will use it for now.
  ##### to ensure that the features are reproducible between the exploratory phase and the deployment in commercial.
  ##### nguyen nhan tao nen su khac biet giua 2 lan chay la vi option "-t first" trong bedtools closest, no se lay nucleosome dau tien ma no gap duoc,
  ##### co the la nucleosome gan nhat hoac la nucleosome xa nhat, tuy thuoc vao option "-t first" hoac "-t last". Khi default la "-t all", no se lay tat ca. 
  cat ${outputdir}/${sampleid}.forward_Nucleosome.bed | awk '{if($5 > 0) print $0}' > ${outputdir}/${sampleid}.full_Nucleosome.bed
  cat ${outputdir}/${sampleid}.reverse_Nucleosome.bed | awk '{if($5 <= 0) print $0}' >> ${outputdir}/${sampleid}.full_Nucleosome.bed
  python convert_full_bed_nucleosome.py ${outputdir}/${sampleid}.full_Nucleosome.bed ${outputdir}/${sampleid}.full_Nucleosome.sorted.bed
  bedtools closest -a ${outputdir}/${sampleid}.full_Nucleosome.sorted.bed -b ${nucleosome_ref} | cut -f1,2,11 > ${outputdir}/${sampleid}.full_Nucleosome.dist.bed
  awk -v OFS='\t' '{$4=$3-$2; print $0}' ${outputdir}/${sampleid}.full_Nucleosome.dist.bed > ${outputdir}/${sampleid}.
  
  #####

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
    paste ${inputfrag} ${outputdir}/forward_nuc.tmp.txt  > ${outputdir}/${sampleid}.modified1.tsv

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
    # rm -rf ${outputdir}/${sampleid}.modified{1,2,3,4}.tsv
fi
